import { useEffect, useMemo, useRef, useState } from "react";
import { useServerFn } from "@tanstack/react-start";
import type { Book } from "../data/books";
import { keywordSearch } from "../lib/keywordSearch";
import { smartSearch } from "../lib/librarySearch.functions";

const DEBOUNCE_MS = 600;
const MIN_CHARS = 2;

type Status =
  | { kind: "idle" }
  | { kind: "searching" }
  | { kind: "done"; count: number }
  | { kind: "error"; message: string };

type LibraryFilterProps = {
  books: Book[];
  /** The visible list, or `null` to show everything. */
  onChange: (visible: Book[] | null) => void;
};

export function LibraryFilter({ books, onChange }: LibraryFilterProps) {
  const search = useServerFn(smartSearch);
  const [query, setQuery] = useState("");
  const [genre, setGenre] = useState<string | null>(null);
  const [aiIds, setAiIds] = useState<string[] | null>(null);
  const [status, setStatus] = useState<Status>({ kind: "idle" });
  const reqId = useRef(0);

  // Pills come from whatever genres exist in the data, most common first.
  const genres = useMemo(() => {
    const counts = new Map<string, number>();
    for (const b of books) for (const g of b.genres ?? []) counts.set(g, (counts.get(g) ?? 0) + 1);
    return [...counts.entries()].sort((a, b) => b[1] - a[1] || a[0].localeCompare(b[0])).map(([g]) => g);
  }, [books]);

  useEffect(() => {
    const q = query.trim();
    const id = ++reqId.current;
    if (q.length < MIN_CHARS) {
      setAiIds(null);
      setStatus({ kind: "idle" });
      return;
    }
    setStatus({ kind: "searching" });
    const t = setTimeout(async () => {
      try {
        const result = await search({ data: { query: q } });
        if (id !== reqId.current) return;
        if (result.status === "error") {
          setStatus({ kind: "error", message: result.message });
          return;
        }
        const ids = result.status === "ok" ? result.ids : keywordSearch(q, books);
        setAiIds(ids);
        setStatus({ kind: "done", count: ids.length });
      } catch {
        if (id !== reqId.current) return;
        // Server unreachable: fall back to plain matching rather than failing.
        const ids = keywordSearch(q, books);
        setAiIds(ids);
        setStatus({ kind: "done", count: ids.length });
      }
    }, DEBOUNCE_MS);
    return () => clearTimeout(t);
  }, [query, books, search]);

  useEffect(() => {
    const matchesGenre = (b: Book) => !genre || (b.genres ?? []).includes(genre);
    if (aiIds) {
      const byId = new Map(books.map((b) => [b.id, b]));
      // Keep the AI's ranking, narrowed to the selected genre.
      onChange(aiIds.map((id) => byId.get(id)).filter((b): b is Book => !!b && matchesGenre(b)));
    } else {
      onChange(genre ? books.filter(matchesGenre) : null);
    }
  }, [aiIds, genre, books, onChange]);

  const pill = (active: boolean) =>
    `shrink-0 rounded-full border px-3.5 py-1.5 text-[13px] transition-colors duration-300 focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-ring/60 ${
      active
        ? "border-primary/60 bg-primary text-primary-foreground"
        : "border-border/80 bg-card/40 text-foreground/75 hover:border-primary/40 hover:text-foreground"
    }`;

  return (
    <div className="w-full max-w-2xl">
      <label className="relative block">
        <span className="sr-only">Search the library</span>
        <svg
          aria-hidden
          viewBox="0 0 24 24"
          className="pointer-events-none absolute left-4 top-1/2 h-4 w-4 -translate-y-1/2 text-muted-foreground"
          fill="none"
          stroke="currentColor"
          strokeWidth="1.6"
        >
          <circle cx="11" cy="11" r="6.5" />
          <path d="m16 16 4 4" strokeLinecap="round" />
        </svg>
        <input
          type="search"
          value={query}
          onChange={(e) => setQuery(e.target.value)}
          placeholder="What are you looking for?"
          className="w-full rounded-full border border-border/80 bg-card/60 py-3 pl-11 pr-5 text-[15px] text-foreground shadow-[0_10px_40px_-24px_var(--shade)] backdrop-blur-sm placeholder:text-muted-foreground/80 [&::-webkit-search-cancel-button]:appearance-none focus:border-primary/50 focus:outline-none focus:ring-2 focus:ring-ring/30"
        />
      </label>
      <p
        aria-live="polite"
        className={`mt-2 h-4 font-mono text-[10px] uppercase tracking-[0.16em] ${
          status.kind === "error" ? "text-destructive" : "text-muted-foreground"
        }`}
      >
        {status.kind === "searching" && "Reading the shelves…"}
        {status.kind === "done" && `${status.count} found`}
        {status.kind === "error" && status.message}
      </p>

      {genres.length > 0 && (
        <div className="no-scrollbar mt-2 flex flex-nowrap gap-2 overflow-x-auto pb-1">
          <button type="button" className={pill(genre === null)} onClick={() => setGenre(null)}>
            All
          </button>
          {genres.map((g) => (
            <button
              key={g}
              type="button"
              aria-pressed={genre === g}
              className={pill(genre === g)}
              onClick={() => setGenre((cur) => (cur === g ? null : g))}
            >
              {g}
            </button>
          ))}
        </div>
      )}
    </div>
  );
}
