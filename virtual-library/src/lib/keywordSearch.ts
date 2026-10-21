import type { Book } from "../data/books";

// Used when no AI key is configured: same UI, simpler ranking.
const STOPWORDS = new Set(
  "a an and the of to in on for with about by from or me my i something some book books novel novels read want looking like that is are".split(
    " ",
  ),
);

export function keywordSearch(query: string, books: Book[]): string[] {
  const terms = query
    .toLowerCase()
    .split(/[^\p{L}\p{N}&]+/u)
    .filter((t) => t.length > 1 && !STOPWORDS.has(t));
  if (terms.length === 0) return [];

  const scored = books.map((book) => {
    const title = book.title.toLowerCase();
    const author = book.author.toLowerCase();
    const genres = (book.genres ?? []).join(" ").toLowerCase();
    const blurb = book.blurb.toLowerCase();
    let score = 0;
    let hits = 0;
    for (const t of terms) {
      const s =
        (title.includes(t) ? 3 : 0) +
        (author.includes(t) ? 3 : 0) +
        (genres.includes(t) ? 2 : 0) +
        (blurb.includes(t) ? 1 : 0);
      if (s > 0) hits++;
      score += s;
    }
    return { id: book.id, score, hits };
  });

  // Prefer books matching every term; otherwise anything that matched at all.
  const all = scored.filter((s) => s.hits === terms.length);
  const pool = all.length > 0 ? all : scored.filter((s) => s.score > 0);
  return pool
    .sort((a, b) => b.score - a.score)
    .slice(0, 20)
    .map((s) => s.id);
}
