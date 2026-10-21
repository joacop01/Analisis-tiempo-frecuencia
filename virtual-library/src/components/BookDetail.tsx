import { useCallback, useEffect, useRef, useState } from "react";
import { createPortal } from "react-dom";
import type { Book } from "../data/books";
import { COVER_W } from "./bookFaces";
import { BookBody } from "./BookSpine";

export type Rect = { left: number; top: number; width: number; height: number };

const PULL_TRANSITION = "transform 900ms cubic-bezier(0.16, 1, 0.3, 1)";
const RETRACT_CLOSE_MS = 620;
const SHELF_POSE = "translate3d(0, 0, 0) scale(1) rotateY(-26deg)";
const MONTH_YEAR = /^[A-Z][a-z]{2} \d{4}$/;

type Layout = {
  narrow: boolean;
  scale: number;
  coverLeft: number;
  coverTop: number;
  coverW: number;
  coverH: number;
};

function computeLayout(vp: { w: number; h: number }, book: Book): Layout {
  const narrow = vp.w < 720;
  const coverH = narrow ? vp.h * 0.42 : Math.min(vp.h * 0.6, 480);
  const scale = coverH / book.height;
  const coverW = COVER_W * scale;
  if (narrow) {
    return { narrow, scale, coverW, coverH, coverLeft: (vp.w - coverW) / 2, coverTop: Math.max(20, vp.h * 0.06) };
  }
  const panelW = Math.min(440, vp.w * 0.38);
  const total = coverW + 56 + panelW;
  return { narrow, scale, coverW, coverH, coverLeft: Math.max(24, (vp.w - total) / 2), coverTop: (vp.h - coverH) / 2 };
}

function useViewport() {
  const [vp, setVp] = useState(() => ({ w: window.innerWidth, h: window.innerHeight }));
  useEffect(() => {
    const onResize = () => setVp({ w: window.innerWidth, h: window.innerHeight });
    window.addEventListener("resize", onResize);
    return () => window.removeEventListener("resize", onResize);
  }, []);
  return vp;
}

/** Flips to true two frames after mount, so the start pose paints before transitioning. */
function useArrived() {
  const [arrived, setArrived] = useState(false);
  useEffect(() => {
    let inner = 0;
    const outer = requestAnimationFrame(() => (inner = requestAnimationFrame(() => setArrived(true))));
    return () => {
      cancelAnimationFrame(outer);
      cancelAnimationFrame(inner);
    };
  }, []);
  return arrived;
}

function PulledBook({ book, rect, layout, closing }: { book: Book; rect: Rect; layout: Layout; closing: boolean }) {
  const out = useArrived() && !closing;
  // Hinge (the spine's right edge) is the transform origin: the front cover swings to face you.
  const dx = layout.coverLeft - (rect.left + book.width);
  const dy = layout.coverTop - rect.top;
  const transform = out ? `translate3d(${dx}px, ${dy}px, 0) scale(${layout.scale}) rotateY(-90deg)` : SHELF_POSE;

  return (
    <div
      className="absolute [transform-style:preserve-3d]"
      style={{
        left: rect.left,
        top: rect.top,
        width: book.width,
        height: book.height,
        transformOrigin: `${book.width}px 0px`,
        transform,
        transition: PULL_TRANSITION,
      }}
    >
      <BookBody book={book} showPages={false} />
    </div>
  );
}

function Stars({ rating }: { rating: number }) {
  if (!rating)
    return <span className="font-mono text-[11px] uppercase tracking-[0.2em] text-muted-foreground">Unrated</span>;
  return (
    <span aria-label={`${rating} of 5 stars`} className="text-lg tracking-[0.2em] text-primary">
      {"★".repeat(rating)}
      <span className="text-muted-foreground/40">{"★".repeat(Math.max(0, 5 - rating))}</span>
    </span>
  );
}

type BookDetailProps = {
  book: Book;
  rect: Rect;
  canStep: boolean;
  onPrev: () => void;
  onNext: () => void;
  onClose: () => void;
};

export function BookDetail({ book, rect, canStep, onPrev, onNext, onClose }: BookDetailProps) {
  const vp = useViewport();
  const shown = useArrived();
  const [closing, setClosing] = useState(false);
  const closeTimer = useRef<ReturnType<typeof setTimeout>>(undefined);
  const dialogRef = useRef<HTMLDivElement>(null);
  const layout = computeLayout(vp, book);

  const retract = useCallback(() => {
    if (closeTimer.current) return;
    setClosing(true);
    closeTimer.current = setTimeout(onClose, RETRACT_CLOSE_MS);
  }, [onClose]);

  useEffect(() => {
    const body = document.body;
    const prevOverflow = body.style.overflow;
    body.dataset.bookOpen = "true";
    body.style.overflow = "hidden";
    dialogRef.current?.focus({ preventScroll: true });
    return () => {
      delete body.dataset.bookOpen;
      body.style.overflow = prevOverflow;
      clearTimeout(closeTimer.current);
    };
  }, []);

  useEffect(() => {
    const onKey = (e: KeyboardEvent) => {
      if (closing) return;
      if (e.key === "Escape") retract();
      else if (canStep && e.key === "ArrowLeft") onPrev();
      else if (canStep && e.key === "ArrowRight") onNext();
      else return;
      e.preventDefault();
    };
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, [closing, canStep, retract, onPrev, onNext]);

  const label = book.recommender
    ? `Recommended by ${book.recommender}`
    : MONTH_YEAR.test(book.finished)
      ? `Finished ${book.finished}`
      : book.finished;

  const panelStyle = layout.narrow
    ? { left: 20, right: 20, top: layout.coverTop + layout.coverH + 28, bottom: 20 }
    : {
        left: layout.coverLeft + layout.coverW + 56,
        width: Math.min(440, vp.w * 0.38),
        top: layout.coverTop,
        height: layout.coverH,
      };

  const control =
    "rounded-full border border-border/80 px-4 py-2 font-mono text-[11px] uppercase tracking-[0.18em] text-foreground/80 transition-colors hover:border-primary/60 hover:text-primary focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-ring/60";

  return createPortal(
    <div
      ref={dialogRef}
      role="dialog"
      aria-modal="true"
      aria-label={`${book.title} by ${book.author}`}
      tabIndex={-1}
      className="fixed inset-0 z-[90] outline-none"
    >
      <div
        className="absolute inset-0 bg-background/70 backdrop-blur-xl transition-opacity duration-700 ease-out"
        style={{ opacity: shown && !closing ? 1 : 0 }}
        onClick={retract}
      />

      {/* Soft ground shadow under the pulled-out cover (kept outside the 3D context). */}
      <div
        aria-hidden
        className="pointer-events-none absolute rounded-sm bg-shade/30 blur-2xl transition-opacity duration-700"
        style={{
          left: layout.coverLeft + 10,
          top: layout.coverTop + 22,
          width: layout.coverW,
          height: layout.coverH,
          opacity: shown && !closing ? 1 : 0,
        }}
      />

      <div
        className="pointer-events-none absolute inset-0"
        style={{
          perspective: "1400px",
          perspectiveOrigin: `${layout.coverLeft + layout.coverW / 2}px ${layout.coverTop + layout.coverH / 2}px`,
        }}
      >
        <PulledBook key={book.id} book={book} rect={rect} layout={layout} closing={closing} />
      </div>

      <div
        className="absolute flex flex-col overflow-y-auto transition-opacity duration-300"
        style={{ ...panelStyle, opacity: closing ? 0 : 1 }}
      >
        <div key={book.id} className="my-auto animate-rise [animation-delay:260ms]">
          {label && <p className="font-mono text-[11px] uppercase tracking-[0.22em] text-muted-foreground">{label}</p>}
          <h2 className="mt-3 font-display text-4xl font-light leading-[1.05] text-foreground md:text-5xl">
            {book.title}
          </h2>
          <p className="mt-2 font-display text-xl italic text-foreground/75">{book.author}</p>
          <p className="mt-3 font-mono text-[11px] uppercase tracking-[0.18em] text-muted-foreground">
            {[book.year || null, book.publisher || null, book.binding === "mass" ? "mass market" : book.binding]
              .filter(Boolean)
              .join(" · ")}
          </p>
          {book.blurb && (
            <p className="mt-5 max-h-[34vh] overflow-y-auto pr-2 text-[15px] leading-relaxed text-foreground/80">
              {book.blurb}
            </p>
          )}
          <div className="mt-5 flex flex-wrap items-center gap-3">
            <Stars rating={book.rating} />
            {book.genres?.map((g) => (
              <span key={g} className="rounded-full bg-foreground/5 px-3 py-1 text-[12px] text-foreground/70">
                {g}
              </span>
            ))}
          </div>
          <div className="mt-7 flex flex-wrap gap-2">
            {canStep && (
              <>
                <button type="button" className={control} onClick={onPrev}>
                  ← Previous
                </button>
                <button type="button" className={control} onClick={onNext}>
                  Next →
                </button>
              </>
            )}
            <button type="button" className={`${control} border-primary/50 text-primary`} onClick={retract}>
              Shelve it
            </button>
          </div>
        </div>
      </div>
    </div>,
    document.body,
  );
}
