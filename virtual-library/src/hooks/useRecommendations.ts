import { useCallback, useEffect, useRef, useState } from "react";
import type { Book } from "../data/books";
import { supabase } from "../lib/supabase";

type RecommendationRow = {
  id: string;
  created_at: string;
  recommender: string;
  note: string | null;
  title: string;
  author: string;
  cover: string;
  year: number;
  publisher: string;
  binding: Book["binding"];
  finish: Book["finish"];
  spine: string;
  band: string | null;
  ink: string;
  face: Book["face"];
  caps: boolean;
  width: number;
  height: number;
  lean: number;
  depth: number;
  wear: number;
};

export type RecommendInput = { book: Book; recommender: string; note: string };

function toBook(row: RecommendationRow): Book {
  return {
    id: row.id,
    title: row.title,
    author: row.author,
    genres: [],
    cover: row.cover,
    year: row.year,
    blurb: row.note ?? "",
    rating: 0,
    finished: `Recommended by ${row.recommender}`,
    recommender: row.recommender,
    publisher: row.publisher,
    binding: row.binding,
    finish: row.finish,
    spine: row.spine,
    band: row.band ?? undefined,
    ink: row.ink,
    face: row.face,
    caps: row.caps,
    width: row.width,
    height: row.height,
    lean: row.lean,
    depth: row.depth,
    wear: row.wear,
  };
}

export function useRecommendations() {
  const [recommendations, setRecommendations] = useState<Book[]>([]);
  const [justAdded, setJustAdded] = useState<string | null>(null);
  const timer = useRef<ReturnType<typeof setTimeout>>(undefined);

  useEffect(() => {
    if (!supabase) return;
    let cancelled = false;
    supabase
      .from("recommendations")
      .select("*")
      .order("created_at", { ascending: false })
      .limit(200)
      .then(({ data, error }) => {
        if (cancelled) return;
        if (error) console.error("Could not load recommendations:", error.message);
        else setRecommendations((data as RecommendationRow[]).map(toBook));
      });
    return () => {
      cancelled = true;
      clearTimeout(timer.current);
    };
  }, []);

  const recommend = useCallback(async ({ book, recommender, note }: RecommendInput) => {
    if (!supabase) throw new Error("Recommendations aren't enabled on this library.");
    const { data, error } = await supabase
      .from("recommendations")
      .insert({
        recommender: recommender.trim(),
        note: note.trim() || null,
        title: book.title,
        author: book.author,
        cover: book.cover,
        year: book.year,
        publisher: book.publisher,
        binding: book.binding,
        finish: book.finish,
        spine: book.spine,
        band: book.band ?? null,
        ink: book.ink,
        face: book.face,
        caps: book.caps ?? false,
        width: book.width,
        height: book.height,
        lean: book.lean,
        depth: book.depth,
        wear: book.wear,
      })
      .select()
      .single();
    if (error) throw new Error(error.message);

    const added = toBook(data as RecommendationRow);
    setRecommendations((prev) => [added, ...prev]);
    setJustAdded(added.id);
    clearTimeout(timer.current);
    timer.current = setTimeout(() => setJustAdded(null), 1400);
  }, []);

  return { enabled: supabase !== null, recommendations, justAdded, recommend };
}
