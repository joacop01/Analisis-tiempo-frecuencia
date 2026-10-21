export type Book = {
  id: string;
  title: string;
  author: string;
  genres?: string[];          // genre tags shown as filter pills
  cover: string;               // real cover art (Open Library URL)
  year: number;
  blurb: string;
  rating: number;              // 0 means unrated
  finished: string;            // e.g. "Jul 2026"
  recommender?: string;        // set only when a visitor recommended it
  publisher: string;
  binding: "hardcover" | "paperback" | "mass";
  finish: "cloth" | "gloss" | "matte";   // spine surface material
  spine: string;               // base color, sampled from the cover's left edge
  band?: string;               // accent pulled from the cover art
  ink: string;                 // lettering color
  face: "serif" | "sans" | "mono";       // lettering style
  caps?: boolean;
  width: number;               // spine width in px
  height: number;              // spine height in px
  lean: number;                // degrees of lean on the shelf
  depth: number;               // how far forward/back the book sits, px
  wear: number;                // 0–1 edge wear and ink fade
  spineImage?: string;
};

// Generated from goodreads_library_export.csv — see README ("Import your
// Goodreads library"). Regenerate with `npm run import:goodreads`; everything
// below this line is overwritten.
export const books: Book[] = [];
