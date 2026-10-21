import { useState } from "react";
import { createFileRoute } from "@tanstack/react-router";
import { LibraryFilter } from "../components/LibraryFilter";
import { RecommendBookDialog } from "../components/RecommendBookDialog";
import { Shelf } from "../components/Shelf";
import { TypedTitle } from "../components/TypedTitle";
import { books, type Book } from "../data/books";
import { useRecommendations } from "../hooks/useRecommendations";

const TITLE = "Your virtual library";
const DESCRIPTION = "A personal archive of the books I've read — pull one off the shelf.";

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

  return (
    <div className="grain min-h-screen overflow-x-hidden">
      <div aria-hidden className="pointer-events-none fixed inset-0 z-0 overflow-hidden">
        <div className="absolute -inset-[25%] animate-drift bg-[radial-gradient(ellipse_at_32%_18%,var(--glow)_0%,transparent_55%),radial-gradient(ellipse_at_78%_70%,color-mix(in_oklch,var(--primary)_14%,transparent)_0%,transparent_50%)]" />
      </div>

      <div className="relative z-10">
        <header className="mx-auto flex max-w-6xl flex-col items-center px-4 pt-16 text-center sm:pt-24">
          <p className="animate-rise font-mono text-[10px] uppercase tracking-[0.32em] text-muted-foreground">
            A personal archive
          </p>
          <div className="mt-5">
            <TypedTitle />
          </div>
          <div className="mt-6 flex animate-rise items-center gap-4 [animation-delay:400ms]">
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
            <div className="mt-9 flex w-full animate-rise justify-center [animation-delay:600ms]">
              <LibraryFilter books={books} onChange={setFiltered} />
            </div>
          )}
        </header>

        <main className="mt-4 animate-rise pb-2 [animation-delay:800ms]">
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

      <RecommendBookDialog open={recommending} onClose={() => setRecommending(false)} onRecommend={recommend} />
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
