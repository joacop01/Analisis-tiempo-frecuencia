import type { Book } from "../data/books";

/** Front-cover depth (how far the book reaches back), in un-scaled spine space. */
export const COVER_W = 178;

export const faceFont: Record<Book["face"], string> = {
  serif: "font-display",
  sans: "font-sans",
  mono: "font-mono",
};
