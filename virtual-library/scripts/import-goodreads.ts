// Regenerates src/data/books.ts from a Goodreads library export.
//
//   npm run import:goodreads -- [path/to/goodreads_library_export.csv]
//
// Reads only "read" / "currently-reading" rows, looks each book up on Open
// Library (cover, blurb, subjects), samples the cover's colours for the spine,
// derives the physical props deterministically, and overwrites the `books`
// array (the `Book` type above it is kept as is).
import { existsSync, readFileSync, writeFileSync } from "node:fs";
import { resolve } from "node:path";
import { fileURLToPath } from "node:url";
import jpeg from "jpeg-js";
import type { Book } from "../src/data/books.ts";
import { derivePhysical, fallbackPalette, paletteFromPixels, slugify, type Palette } from "../src/lib/bookPhysics.ts";

const BOOKS_FILE = fileURLToPath(new URL("../src/data/books.ts", import.meta.url));
const OL = "https://openlibrary.org";
const HEADERS = { "User-Agent": "VirtualLibraryImport/1.0 (personal bookshelf generator)" };
const CONCURRENCY = 3;
const MONTHS = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"];
const GENRES = ["Fiction", "Nonfiction", "Sci-Fi", "Mystery & Thriller", "Fantasy", "Romance"] as const;

type Row = Record<string, string>;

// ── CSV ─────────────────────────────────────────────────────────────────────

function parseCsv(text: string): Row[] {
  const rows: string[][] = [];
  let row: string[] = [];
  let field = "";
  let quoted = false;
  const src = text.replace(/^﻿/, "");
  for (let i = 0; i < src.length; i++) {
    const c = src[i];
    if (quoted) {
      if (c === '"' && src[i + 1] === '"') {
        field += '"';
        i++;
      } else if (c === '"') quoted = false;
      else field += c;
    } else if (c === '"') quoted = true;
    else if (c === ",") {
      row.push(field);
      field = "";
    } else if (c === "\n" || c === "\r") {
      if (c === "\r" && src[i + 1] === "\n") i++;
      row.push(field);
      rows.push(row);
      row = [];
      field = "";
    } else field += c;
  }
  if (field || row.length) {
    row.push(field);
    rows.push(row);
  }
  const [header, ...body] = rows;
  return body
    .filter((r) => r.some((v) => v.trim()))
    .map((r) => Object.fromEntries(header.map((h, i) => [h.trim(), (r[i] ?? "").trim()])));
}

// Goodreads writes ISBNs as ="0446310786".
const cleanIsbn = (v = "") => v.replace(/[="\s-]/g, "");

function formatMonth(goodreadsDate: string): string {
  const m = goodreadsDate.match(/^(\d{4})\/(\d{1,2})/);
  return m ? `${MONTHS[Number(m[2]) - 1]} ${m[1]}` : "";
}

// "The Name of the Wind (The Kingkiller Chronicle, #1)" → "The Name of the Wind"
const stripSeries = (title: string) => title.replace(/\s*\([^()]*#\d+(\.\d+)?\)\s*$/, "").trim();

// ── Open Library ────────────────────────────────────────────────────────────

async function getJson<T>(url: string): Promise<T | null> {
  for (let attempt = 0; attempt < 4; attempt++) {
    const res = await fetch(url, { headers: HEADERS }).catch(() => null);
    if (res?.ok) return (await res.json()) as T;
    if (res && res.status !== 429 && res.status < 500) return null;
    await new Promise((r) => setTimeout(r, 1000 * 2 ** attempt));
  }
  return null;
}

type Edition = { covers?: number[]; works?: { key: string }[]; number_of_pages?: number; publishers?: string[] };
type Work = { description?: string | { value: string }; subjects?: string[]; covers?: number[] };
type SearchDoc = {
  key: string;
  cover_i?: number;
  first_publish_year?: number;
  number_of_pages_median?: number;
  publisher?: string[];
  subject?: string[];
};

type Found = {
  workKey?: string;
  coverId?: number;
  pages?: number;
  publisher?: string;
  firstYear?: number;
  description?: string;
  subjects: string[];
};

async function lookUp(row: Row): Promise<Found> {
  const found: Found = { subjects: [] };
  const isbn = cleanIsbn(row["ISBN13"]) || cleanIsbn(row["ISBN"]);

  if (isbn) {
    const edition = await getJson<Edition>(`${OL}/isbn/${isbn}.json`);
    if (edition) {
      found.workKey = edition.works?.[0]?.key;
      found.coverId = edition.covers?.find((c) => c > 0);
      found.pages = edition.number_of_pages;
      found.publisher = edition.publishers?.[0];
    }
  }

  const params = new URLSearchParams({
    title: stripSeries(row["Title"]),
    author: row["Author"],
    limit: "1",
    fields: "key,cover_i,first_publish_year,number_of_pages_median,publisher,subject",
  });
  const search = await getJson<{ docs?: SearchDoc[] }>(`${OL}/search.json?${params}`);
  const doc = search?.docs?.[0];
  if (doc) {
    found.workKey ??= doc.key;
    found.coverId ??= doc.cover_i;
    found.pages ??= doc.number_of_pages_median;
    found.publisher ??= doc.publisher?.[0];
    found.firstYear = doc.first_publish_year;
    found.subjects.push(...(doc.subject ?? []));
  }

  if (found.workKey) {
    const work = await getJson<Work>(`${OL}${found.workKey}.json`);
    if (work) {
      const d = work.description;
      found.description = typeof d === "string" ? d : d?.value;
      found.subjects.push(...(work.subjects ?? []));
      found.coverId ??= work.covers?.find((c) => c > 0);
    }
  }
  return found;
}

async function samplePalette(coverId: number): Promise<Palette | null> {
  try {
    const res = await fetch(`https://covers.openlibrary.org/b/id/${coverId}-M.jpg`, { headers: HEADERS });
    if (!res.ok) return null;
    const img = jpeg.decode(new Uint8Array(await res.arrayBuffer()), { useTArray: true, formatAsRGBA: true });
    return paletteFromPixels(img.data, img.width, img.height);
  } catch {
    return null;
  }
}

// ── Shaping ─────────────────────────────────────────────────────────────────

function cleanBlurb(raw = ""): string {
  const text = raw
    .split(/\n-{3,}/)[0]
    .replace(/\(\[source\]\[\d+\]\)/gi, "")
    .replace(/\[([^\]]+)\]\([^)]*\)/g, "$1")
    .replace(/\[([^\]]+)\]\[\d+\]/g, "$1")
    .replace(/^\s*\[\d+\]:.*$/gm, "")
    .replace(/\s+/g, " ")
    .trim();
  if (text.length <= 420) return text;
  const cut = text.slice(0, 420);
  const end = cut.lastIndexOf(". ");
  return end > 200 ? cut.slice(0, end + 1) : `${cut.replace(/\s+\S*$/, "")}…`;
}

function assignGenres(shelves: string, subjects: string[]): string[] {
  const text = [shelves, ...subjects].join(" | ").toLowerCase();
  const genres: (typeof GENRES)[number][] = [];
  if (/mystery|thriller|detective|crime|suspense|murder/.test(text)) genres.push("Mystery & Thriller");
  if (/science fiction|sci-fi|scifi|space opera|dystopia|cyberpunk/.test(text)) genres.push("Sci-Fi");
  if (/fantasy|magic|dragons|fae\b|faerie|fairies|wizards|witches/.test(text)) genres.push("Fantasy");
  if (/romance|love stories|romantic/.test(text)) genres.push("Romance");
  if (genres.length) return genres.slice(0, 3);
  if (/non-?fiction|biography|memoir|history|essays|self-help|philosophy|science\b/.test(text)) return ["Nonfiction"];
  // "Literary Fiction" and friends → Fiction; anything untagged defaults to Nonfiction.
  return /(?<!non-?)fiction|novel|literary/.test(text) ? ["Fiction"] : ["Nonfiction"];
}

function bookLiteral(book: Book): string {
  const lines = Object.entries(book)
    .filter(([, v]) => v !== undefined)
    .map(([k, v]) => `    ${k}: ${JSON.stringify(v)},`);
  return `  {\n${lines.join("\n")}\n  }`;
}

// ── Main ────────────────────────────────────────────────────────────────────

async function main() {
  const csvPath = resolve(process.argv[2] ?? "goodreads_library_export.csv");
  if (!existsSync(csvPath)) {
    console.error(
      `Can't find ${csvPath}.\nExport it from Goodreads (My Books → Import and export → Export Library) and pass its path.`,
    );
    process.exit(1);
  }

  const rows = parseCsv(readFileSync(csvPath, "utf8")).filter((r) =>
    ["read", "currently-reading"].includes(r["Exclusive Shelf"]),
  );
  // Newest finished first; currently-reading leads.
  const sortKey = (r: Row) =>
    r["Exclusive Shelf"] === "currently-reading" ? "9999" : r["Date Read"] || `0000 ${r["Date Added"]}`;
  rows.sort((a, b) => sortKey(b).localeCompare(sortKey(a)));
  console.log(`Importing ${rows.length} books from ${csvPath}…`);

  const out: Book[] = new Array(rows.length);
  let next = 0;
  let done = 0;
  const worker = async () => {
    while (next < rows.length) {
      const i = next++;
      const row = rows[i];
      const found = await lookUp(row);
      const key = found.workKey ?? (cleanIsbn(row["ISBN13"]) || `goodreads:${row["Book Id"]}`);
      const palette = (found.coverId && (await samplePalette(found.coverId))) || fallbackPalette(key);
      const pages = Number(row["Number of Pages"]) || found.pages;
      const title = stripSeries(row["Title"]);

      out[i] = {
        id: `${slugify(title).slice(0, 28)}-${i}`,
        title,
        author: row["Author"],
        genres: assignGenres(row["Bookshelves"] ?? "", found.subjects),
        cover: found.coverId ? `https://covers.openlibrary.org/b/id/${found.coverId}-L.jpg` : "",
        year: Number(row["Original Publication Year"]) || Number(row["Year Published"]) || found.firstYear || 0,
        blurb: cleanBlurb(found.description),
        rating: Number(row["My Rating"]) || 0,
        finished:
          formatMonth(row["Date Read"]) || (row["Exclusive Shelf"] === "currently-reading" ? "Currently reading" : ""),
        publisher: row["Publisher"] || found.publisher || "",
        ...palette,
        ...derivePhysical(key, pages),
      };
      done++;
      console.log(`[${done}/${rows.length}] ${title}${found.coverId ? "" : "  (no cover found)"}`);
    }
  };
  await Promise.all(Array.from({ length: CONCURRENCY }, worker));

  const source = readFileSync(BOOKS_FILE, "utf8");
  const marker = source.indexOf("export const books");
  if (marker < 0) throw new Error(`Couldn't find "export const books" in ${BOOKS_FILE}`);
  const body = out.length ? `[\n${out.map(bookLiteral).join(",\n")},\n]` : "[]";
  writeFileSync(BOOKS_FILE, `${source.slice(0, marker)}export const books: Book[] = ${body};\n`);
  console.log(`\nWrote ${out.length} books to src/data/books.ts`);
}

main().catch((err) => {
  console.error(err);
  process.exit(1);
});
