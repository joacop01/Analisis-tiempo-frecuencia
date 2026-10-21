import { useCallback, useEffect, useRef, useState } from "react";
import { createFileRoute } from "@tanstack/react-router";
import { LibraryFilter } from "../components/LibraryFilter";
import { ENTRANCE_STEP_MS } from "../components/motion";
import { RecommendBookDialog } from "../components/RecommendBookDialog";
import { Shelf } from "../components/Shelf";
import { TypedTitle } from "../components/TypedTitle";
import { books, type Book } from "../data/books";
import { useRecommendations } from "../hooks/useRecommendations";

const TITLE = "Your virtual library";
const DESCRIPTION = "A personal archive of the books I've read — pull one off the shelf.";

// The haze drifts along this path (eased in and out, back and forth) every 32s.
const DRIFT_MS = 32000;
const DRIFT_PATH = [
  { x: -2, y: 0, s: 1 },
  { x: 2, y: -2, s: 1.06 },
  { x: -1, y: 1.5, s: 1.02 },
];
// It moves well under a pixel per frame, so it is redrawn only a couple of times a second: a full-screen layer
// animating every frame would make the compositor redraw the whole 3D shelf every frame, even at rest. It is
// also kept off its own compositor layer (a 2D transform, no will-change): beneath the overlay-blended grain, a
// separate layer there turns every small change elsewhere, like the caret's blink, into a full-screen redraw.
const DRIFT_STEP_MS = 500;

export const Route = createFileRoute("/")({
  head: () => ({
    meta: [
      { title: TITLE },
      { name: "description", content: DESCRIPTION },
      { property: "og:title", content: TITLE },
      { property: "og:description", content: DESCRIPTION },
      { property: "og:type", content: "website" },
      { name: "twitter:card", content: "summary_large_image" },
      { name: "twitter:title", content: TITLE },
      { name: "twitter:description", content: DESCRIPTION },
    ],
  }),
  component: Home,
});

function Home() {
  const [filtered, setFiltered] = useState<Book[] | null>(null);
  const [recommending, setRecommending] = useState(false);
  const { enabled, recommendations, justAdded, recommend } = useRecommendations();
  const shown = filtered ?? books;
  const closeRecommend = useCallback(() => setRecommending(false), []);

  return (
    <div className="grain min-h-screen overflow-x-hidden">
      <Haze />

      <div className="relative z-10">
        <header className="mx-auto flex max-w-6xl flex-col items-center px-4 pt-16 text-center sm:pt-24">
          <p className="animate-rise font-mono text-[10px] uppercase tracking-[0.32em] text-muted-foreground">
            A personal archive
          </p>
          <div className="mt-5">
            <TypedTitle />
          </div>
          <div
            className="mt-6 flex animate-rise items-center gap-4"
            style={{ animationDelay: `${ENTRANCE_STEP_MS}ms` }}
          >
            <span className="font-mono text-[11px] uppercase tracking-[0.2em] text-muted-foreground">
              {books.length} {books.length === 1 ? "volume" : "volumes"}
            </span>
            {enabled && (
              <>
                <span aria-hidden className="h-1 w-1 rounded-full bg-muted-foreground/50" />
                <button
                  type="button"
                  onClick={() => setRecommending(true)}
                  className="rounded-full border border-primary/40 px-4 py-1.5 font-mono text-[11px] uppercase tracking-[0.2em] text-primary transition-colors hover:bg-primary hover:text-primary-foreground focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-ring/60"
                >
                  Recommend a book
                </button>
              </>
            )}
          </div>
          {books.length > 0 && (
            <div
              className="mt-9 flex w-full animate-rise justify-center"
              style={{ animationDelay: `${ENTRANCE_STEP_MS * 2}ms` }}
            >
              <LibraryFilter books={books} onChange={setFiltered} />
            </div>
          )}
        </header>

        <main className="mt-4 animate-rise pb-2" style={{ animationDelay: `${ENTRANCE_STEP_MS * 3}ms` }}>
          {books.length === 0 ? (
            <EmptyLibrary />
          ) : shown.length === 0 ? (
            <p className="py-28 text-center font-display text-2xl italic text-muted-foreground">No books match</p>
          ) : (
            <Shelf books={shown} arrowKeys />
          )}
        </main>

        {recommendations.length > 0 && (
          <section className="pb-16 pt-2">
            <div className="mx-auto flex max-w-6xl items-baseline justify-center gap-3 px-4">
              <h2 className="font-display text-2xl font-light italic text-foreground/85">Recommended to me</h2>
              <span className="font-mono text-[10px] uppercase tracking-[0.2em] text-muted-foreground">
                {recommendations.length}
              </span>
            </div>
            <Shelf books={recommendations} justAdded={justAdded} />
          </section>
        )}
      </div>

      <RecommendBookDialog open={recommending} onClose={closeRecommend} onRecommend={recommend} />
    </div>
  );
}

/** Ease-in-out (cubic-bezier(0.42, 0, 0.58, 1)), as the CSS keyword. */
function easeInOut(p: number) {
  let lo = 0;
  let hi = 1;
  for (let i = 0; i < 24; i++) {
    const t = (lo + hi) / 2;
    const x = 3 * (1 - t) ** 2 * t * 0.42 + 3 * (1 - t) * t * t * 0.58 + t ** 3;
    if (x < p) lo = t;
    else hi = t;
  }
  const t = (lo + hi) / 2;
  return 3 * t * t - 2 * t ** 3;
}

function driftAt(ms: number) {
  const cycle = Math.floor(ms / DRIFT_MS);
  let p = (ms % DRIFT_MS) / DRIFT_MS;
  if (cycle % 2) p = 1 - p;
  const half = p < 0.5 ? 0 : 1;
  const e = easeInOut(half ? (p - 0.5) * 2 : p * 2);
  const [a, b] = [DRIFT_PATH[half]!, DRIFT_PATH[half + 1]!];
  const at = (k: "x" | "y" | "s") => (a[k] + (b[k] - a[k]) * e).toFixed(4);
  return `translate(${at("x")}%, ${at("y")}%) scale(${at("s")})`;
}

function Haze() {
  const ref = useRef<HTMLDivElement>(null);

  useEffect(() => {
    const el = ref.current;
    if (!el) return;
    // Reduced motion holds the haze still from the first paint on (see its class); this stops the clock too, and
    // follows the setting when it changes.
    const reduce = window.matchMedia("(prefers-reduced-motion: reduce)");
    // The drift's own clock, which stands still while a dialog covers the page.
    let clock = performance.now();
    let last = clock;
    let id: ReturnType<typeof setInterval> | undefined;
    const tick = () => {
      const now = performance.now();
      if (!document.querySelector('[aria-modal="true"]')) clock += now - last;
      last = now;
      el.style.transform = driftAt(clock);
    };
    const start = () => {
      clearInterval(id);
      if (reduce.matches) return;
      last = performance.now();
      tick();
      id = setInterval(tick, DRIFT_STEP_MS);
    };
    start();
    reduce.addEventListener("change", start);
    return () => {
      clearInterval(id);
      reduce.removeEventListener("change", start);
    };
  }, []);

  return (
    <div aria-hidden className="pointer-events-none fixed inset-0 z-0 overflow-hidden">
      <div
        ref={ref}
        className="absolute -inset-[25%] motion-reduce:transform-none! bg-[radial-gradient(ellipse_at_32%_18%,var(--glow)_0%,transparent_55%),radial-gradient(ellipse_at_78%_70%,color-mix(in_oklch,var(--primary)_14%,transparent)_0%,transparent_50%)]"
        style={{ transform: driftAt(0) }}
      />
    </div>
  );
}

function EmptyLibrary() {
  return (
    <div className="mx-auto max-w-md px-4 py-24 text-center">
      <div
        aria-hidden
        className="mx-auto h-px w-64 bg-gradient-to-r from-transparent via-foreground/25 to-transparent"
      />
      <p className="mt-8 font-display text-3xl font-light italic text-foreground/80">The shelves are waiting.</p>
      <p className="mt-3 text-[15px] leading-relaxed text-muted-foreground">
        Export your library from Goodreads and run the import to fill them with your books.
      </p>
    </div>
  );
}
