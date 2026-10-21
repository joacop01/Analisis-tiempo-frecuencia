import type { Book } from "../data/books";
import { derivePhysical, fallbackPalette, paletteFromPixels, slugify, type Palette } from "./bookPhysics";

export type OpenLibraryResult = {
  key: string;
  title: string;
  author: string;
  year: number;
  cover: string;
  thumb: string;
  pages?: number;
  publisher: string;
};

type SearchDoc = {
  key: string;
  title?: string;
  author_name?: string[];
  first_publish_year?: number;
  cover_i?: number;
  number_of_pages_median?: number;
  publisher?: string[];
};

export const coverUrl = (id: number, size: "S" | "M" | "L" = "L") =>
  `https://covers.openlibrary.org/b/id/${id}-${size}.jpg`;

export async function searchBooks(q: string, signal?: AbortSignal): Promise<OpenLibraryResult[]> {
  const params = new URLSearchParams({
    q,
    limit: "12",
    fields: "key,title,author_name,first_publish_year,cover_i,number_of_pages_median,publisher",
  });
  const res = await fetch(`https://openlibrary.org/search.json?${params}`, { signal });
  if (!res.ok) throw new Error(`Open Library search failed (${res.status})`);
  const json = (await res.json()) as { docs?: SearchDoc[] };

  return (json.docs ?? [])
    .filter((d) => d.title && d.cover_i)
    .map((d) => ({
      key: d.key,
      title: d.title!,
      author: d.author_name?.[0] ?? "Unknown author",
      year: d.first_publish_year ?? 0,
      cover: coverUrl(d.cover_i!, "L"),
      thumb: coverUrl(d.cover_i!, "S"),
      pages: d.number_of_pages_median,
      publisher: d.publisher?.[0] ?? "",
    }));
}

export function readCoverPalette(src: string): Promise<Palette> {
  return new Promise((resolve, reject) => {
    const img = new Image();
    img.crossOrigin = "anonymous";
    img.onload = () => {
      try {
        const w = 80;
        const h = Math.max(1, Math.round((img.naturalHeight / img.naturalWidth) * w));
        const canvas = document.createElement("canvas");
        canvas.width = w;
        canvas.height = h;
        const ctx = canvas.getContext("2d", { willReadFrequently: true });
        if (!ctx) throw new Error("no 2d context");
        ctx.drawImage(img, 0, 0, w, h);
        resolve(paletteFromPixels(ctx.getImageData(0, 0, w, h).data, w, h));
      } catch (err) {
        reject(err);
      }
    };
    img.onerror = () => reject(new Error("cover failed to load"));
    img.src = src;
  });
}

export async function buildBook(result: OpenLibraryResult): Promise<Book> {
  const palette = await readCoverPalette(result.cover).catch(() => fallbackPalette(result.key));
  return {
    id: `${slugify(result.title).slice(0, 28)}-${slugify(result.key)}`,
    title: result.title,
    author: result.author,
    genres: [],
    cover: result.cover,
    year: result.year,
    blurb: "",
    rating: 0,
    finished: "",
    publisher: result.publisher,
    ...palette,
    ...derivePhysical(result.key, result.pages),
  };
}
