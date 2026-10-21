import { useEffect, useRef, useState, type FormEvent } from "react";
import { createPortal } from "react-dom";
import type { RecommendInput } from "../hooks/useRecommendations";
import { buildBook, searchBooks, type OpenLibraryResult } from "../lib/openLibrary";

const DEBOUNCE_MS = 280;
const MIN_CHARS = 2;

type RecommendBookDialogProps = {
  open: boolean;
  onClose: () => void;
  onRecommend: (input: RecommendInput) => Promise<void>;
};

export function RecommendBookDialog({ open, onClose, onRecommend }: RecommendBookDialogProps) {
  if (!open) return null;
  return createPortal(<DialogBody onClose={onClose} onRecommend={onRecommend} />, document.body);
}

function DialogBody({ onClose, onRecommend }: Omit<RecommendBookDialogProps, "open">) {
  const [query, setQuery] = useState("");
  const [results, setResults] = useState<OpenLibraryResult[]>([]);
  const [searching, setSearching] = useState(false);
  const [searchError, setSearchError] = useState<string | null>(null);
  const [picked, setPicked] = useState<OpenLibraryResult | null>(null);
  const [name, setName] = useState("");
  const [note, setNote] = useState("");
  const [submitting, setSubmitting] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const inputRef = useRef<HTMLInputElement>(null);

  useEffect(() => {
    inputRef.current?.focus();
    const onKey = (e: KeyboardEvent) => e.key === "Escape" && onClose();
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, [onClose]);

  // Live, debounced, abortable Open Library search.
  useEffect(() => {
    const q = query.trim();
    if (q.length < MIN_CHARS) {
      setResults([]);
      setSearching(false);
      setSearchError(null);
      return;
    }
    const controller = new AbortController();
    setSearching(true);
    const t = setTimeout(() => {
      searchBooks(q, controller.signal)
        .then((r) => {
          setResults(r);
          setSearchError(null);
        })
        .catch((err: unknown) => {
          if ((err as Error).name !== "AbortError") setSearchError("Couldn't reach Open Library. Try again?");
        })
        .finally(() => !controller.signal.aborted && setSearching(false));
    }, DEBOUNCE_MS);
    return () => {
      clearTimeout(t);
      controller.abort();
    };
  }, [query]);

  const nameOk = name.trim().length >= 1 && name.trim().length <= 60;
  const noteOk = note.length <= 500;

  const submit = async (e: FormEvent) => {
    e.preventDefault();
    if (!picked || !nameOk || !noteOk || submitting) return;
    setSubmitting(true);
    setError(null);
    try {
      const book = await buildBook(picked);
      await onRecommend({ book, recommender: name, note });
      onClose();
    } catch (err) {
      setError(err instanceof Error ? err.message : "Something went wrong — please try again.");
      setSubmitting(false);
    }
  };

  const field =
    "w-full rounded-xl border border-border/80 bg-background/60 px-4 py-2.5 text-[15px] text-foreground placeholder:text-muted-foreground/80 focus:border-primary/50 focus:outline-none focus:ring-2 focus:ring-ring/30";

  return (
    <div className="fixed inset-0 z-[95] flex items-end justify-center p-4 sm:items-center">
      <div className="absolute inset-0 animate-fade-in bg-background/60 backdrop-blur-md" onClick={onClose} />
      <div
        role="dialog"
        aria-modal="true"
        aria-labelledby="recommend-title"
        className="relative flex max-h-[88vh] w-full max-w-lg animate-rise flex-col overflow-hidden rounded-2xl border border-border/70 bg-card shadow-[0_40px_100px_-40px_var(--shade)]"
      >
        <div className="flex items-start justify-between gap-4 px-6 pb-2 pt-6">
          <div>
            <p className="font-mono text-[10px] uppercase tracking-[0.22em] text-muted-foreground">For the shelf</p>
            <h2 id="recommend-title" className="mt-1 font-display text-3xl font-light italic text-card-foreground">
              Recommend a book
            </h2>
          </div>
          <button
            type="button"
            onClick={onClose}
            aria-label="Close"
            className="rounded-full p-2 text-muted-foreground transition-colors hover:text-foreground focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-ring/60"
          >
            <svg viewBox="0 0 24 24" className="h-4 w-4" fill="none" stroke="currentColor" strokeWidth="1.6">
              <path d="M6 6l12 12M18 6 6 18" strokeLinecap="round" />
            </svg>
          </button>
        </div>

        {!picked ? (
          <div className="flex min-h-0 flex-col px-6 pb-6">
            <input
              ref={inputRef}
              type="search"
              value={query}
              onChange={(e) => setQuery(e.target.value)}
              placeholder="Search by title or author…"
              className={`${field} mt-3`}
            />
            <p className="mt-2 h-4 font-mono text-[10px] uppercase tracking-[0.16em] text-muted-foreground">
              {searching ? (
                "Searching Open Library…"
              ) : searchError ? (
                <span className="text-destructive">{searchError}</span>
              ) : (
                ""
              )}
            </p>
            <ul className="mt-2 min-h-0 flex-1 space-y-1 overflow-y-auto">
              {results.map((r) => (
                <li key={r.key}>
                  <button
                    type="button"
                    onClick={() => setPicked(r)}
                    className="group flex w-full items-center gap-4 rounded-xl px-2 py-2 text-left transition-colors hover:bg-foreground/5 focus-visible:bg-foreground/5 focus-visible:outline-none"
                  >
                    <img src={r.thumb} alt="" className="h-14 w-10 shrink-0 rounded-[2px] bg-muted object-cover" />
                    <span className="min-w-0 flex-1">
                      <span className="block truncate font-display text-lg leading-tight text-card-foreground">
                        {r.title}
                      </span>
                      <span className="block truncate text-[13px] text-muted-foreground">
                        {r.author}
                        {r.year ? ` · ${r.year}` : ""}
                      </span>
                    </span>
                    <span className="font-mono text-[10px] uppercase tracking-[0.18em] text-primary opacity-0 transition-opacity group-hover:opacity-100 group-focus-visible:opacity-100">
                      Pick
                    </span>
                  </button>
                </li>
              ))}
            </ul>
          </div>
        ) : (
          <form onSubmit={submit} className="overflow-y-auto px-6 pb-6">
            <div className="mt-3 flex gap-4">
              <img
                src={picked.cover}
                alt=""
                className="h-36 w-24 shrink-0 rounded-[3px] bg-muted object-cover shadow-[0_14px_30px_-16px_var(--shade)]"
              />
              <div className="min-w-0">
                <p className="font-display text-2xl leading-tight text-card-foreground">{picked.title}</p>
                <p className="mt-1 text-[14px] text-muted-foreground">
                  {picked.author}
                  {picked.year ? ` · ${picked.year}` : ""}
                </p>
                <button
                  type="button"
                  onClick={() => setPicked(null)}
                  className="mt-3 font-mono text-[10px] uppercase tracking-[0.18em] text-primary hover:underline"
                >
                  ← Choose another
                </button>
              </div>
            </div>

            <label className="mt-6 block">
              <span className="font-mono text-[10px] uppercase tracking-[0.18em] text-muted-foreground">Your name</span>
              <input
                value={name}
                onChange={(e) => setName(e.target.value)}
                maxLength={60}
                required
                className={`${field} mt-1.5`}
              />
            </label>
            <label className="mt-4 block">
              <span className="font-mono text-[10px] uppercase tracking-[0.18em] text-muted-foreground">
                Why should I read it?
              </span>
              <textarea
                value={note}
                onChange={(e) => setNote(e.target.value)}
                maxLength={500}
                rows={3}
                className={`${field} mt-1.5 resize-none`}
              />
              <span className="mt-1 block text-right font-mono text-[10px] text-muted-foreground">
                {note.length}/500
              </span>
            </label>

            {error && <p className="mt-2 text-[14px] text-destructive">{error}</p>}

            <button
              type="submit"
              disabled={!nameOk || !noteOk || submitting}
              className="mt-4 w-full rounded-full bg-primary px-5 py-3 font-mono text-[11px] uppercase tracking-[0.2em] text-primary-foreground transition-opacity disabled:opacity-50"
            >
              {submitting ? "Shelving…" : "Put it on the shelf"}
            </button>
          </form>
        )}
      </div>
    </div>
  );
}
