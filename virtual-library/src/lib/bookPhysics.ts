// Deterministic "physical" properties of a book, shared by the browser
// (visitor recommendations) and the Goodreads import script.
import type { Book } from "../data/books";

export const INK_DARK = "#241f19";
export const INK_LIGHT = "#faf7f0";

export type Palette = { spine: string; band?: string; ink: string };

export type Physical = Pick<
  Book,
  "binding" | "finish" | "face" | "caps" | "width" | "height" | "lean" | "depth" | "wear"
>;

/** FNV-1a — stable 32-bit hash of a string. */
export function hashString(input: string): number {
  let h = 0x811c9dc5;
  for (let i = 0; i < input.length; i++) {
    h ^= input.charCodeAt(i);
    h = Math.imul(h, 0x01000193);
  }
  return h >>> 0;
}

/** mulberry32 — a tiny seeded PRNG returning floats in [0, 1). */
export function seeded(seed: number): () => number {
  let a = seed >>> 0;
  return () => {
    a = (a + 0x6d2b79f5) >>> 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

const clamp = (v: number, lo: number, hi: number) => Math.min(hi, Math.max(lo, v));
const round1 = (v: number) => Math.round(v * 10) / 10;
const round2 = (v: number) => Math.round(v * 100) / 100;

/** Binding, finish, size, lean, depth, wear and lettering, all from one key. */
export function derivePhysical(key: string, pages: number | undefined): Physical {
  const rand = seeded(hashString(key));
  const p = pages && pages > 0 ? pages : 300;

  const binding: Book["binding"] = p > 420 ? "hardcover" : p < 260 ? "mass" : "paperback";
  const finish: Book["finish"] = binding === "hardcover" ? "cloth" : rand() < 0.5 ? "gloss" : "matte";
  const [lo, hi] = binding === "hardcover" ? [236, 254] : binding === "mass" ? [196, 210] : [214, 230];
  const height = Math.round(lo + rand() * (hi - lo));
  const width = Math.round(clamp(p * 0.055 + (rand() * 4 - 2), 16, 58));
  // Most books stand straight; a few rest against their neighbour.
  const lean = rand() < 0.72 ? 0 : round1(-(1 + rand() * 4));
  const depth = round1(rand() * 14 - 7);
  const wear = round2(rand() * 0.35);
  const faceRoll = rand();
  const face: Book["face"] = faceRoll < 0.55 ? "serif" : faceRoll < 0.85 ? "sans" : "mono";
  const caps = rand() < 0.35;

  return { binding, finish, face, caps, width, height, lean, depth, wear };
}

// ── colour helpers ──────────────────────────────────────────────────────────

const toHex = (n: number) =>
  Math.round(clamp(n, 0, 255))
    .toString(16)
    .padStart(2, "0");
export const rgbToHex = (r: number, g: number, b: number) => `#${toHex(r)}${toHex(g)}${toHex(b)}`;

export function hexToRgb(hex: string): [number, number, number] {
  const h = hex.replace("#", "");
  const full =
    h.length === 3
      ? h
          .split("")
          .map((c) => c + c)
          .join("")
      : h;
  const n = parseInt(full, 16);
  return [(n >> 16) & 255, (n >> 8) & 255, n & 255];
}

/** Relative luminance, 0 (black) – 1 (white). */
export function luminance(r: number, g: number, b: number): number {
  const lin = (c: number) => {
    const s = c / 255;
    return s <= 0.03928 ? s / 12.92 : ((s + 0.055) / 1.055) ** 2.4;
  };
  return 0.2126 * lin(r) + 0.7152 * lin(g) + 0.0722 * lin(b);
}

export const inkFor = (spine: string) => (luminance(...hexToRgb(spine)) > 0.55 ? INK_DARK : INK_LIGHT);

function hslToHex(h: number, s: number, l: number): string {
  const a = s * Math.min(l, 1 - l);
  const f = (n: number) => {
    const k = (n + h / 30) % 12;
    return l - a * Math.max(-1, Math.min(k - 3, 9 - k, 1));
  };
  return rgbToHex(f(0) * 255, f(8) * 255, f(4) * 255);
}

/** A muted cloth-like colour when there is no cover to sample. */
export function fallbackPalette(key: string): Palette {
  const rand = seeded(hashString(`palette:${key}`));
  const hue = Math.floor(rand() * 360);
  const spine = hslToHex(hue, 0.22 + rand() * 0.22, 0.22 + rand() * 0.2);
  const band = hslToHex((hue + 25 + rand() * 30) % 360, 0.45, 0.62);
  return { spine, band, ink: inkFor(spine) };
}

/**
 * Sample a palette from RGBA pixels of a cover: the left 6% averaged into
 * `spine`, the most saturated mid-luminance pixel into `band`.
 */
export function paletteFromPixels(data: ArrayLike<number>, w: number, h: number): Palette {
  const edgeW = Math.max(1, Math.round(w * 0.06));
  let r = 0,
    g = 0,
    b = 0,
    n = 0;
  for (let y = 0; y < h; y++) {
    for (let x = 0; x < edgeW; x++) {
      const i = (y * w + x) * 4;
      r += data[i];
      g += data[i + 1];
      b += data[i + 2];
      n++;
    }
  }
  const spine = rgbToHex(r / n, g / n, b / n);

  let best = -1;
  let band: string | undefined;
  for (let i = 0; i < w * h * 4; i += 4 * 3) {
    const pr = data[i],
      pg = data[i + 1],
      pb = data[i + 2];
    const max = Math.max(pr, pg, pb),
      min = Math.min(pr, pg, pb);
    const lum = luminance(pr, pg, pb);
    if (lum < 0.08 || lum > 0.7 || max === 0) continue;
    const sat = (max - min) / max;
    if (sat > best) {
      best = sat;
      band = rgbToHex(pr, pg, pb);
    }
  }
  return { spine, band: best > 0.25 ? band : undefined, ink: inkFor(spine) };
}

export function slugify(text: string): string {
  return text
    .toLowerCase()
    .normalize("NFKD")
    .replace(/[̀-ͯ]/g, "")
    .replace(/[^a-z0-9]+/g, "-")
    .replace(/^-+/, "");
}
