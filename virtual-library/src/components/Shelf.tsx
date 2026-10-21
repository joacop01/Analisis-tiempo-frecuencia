import { useCallback, useEffect, useLayoutEffect, useMemo, useRef, useState } from "react";
import type { Book } from "../data/books";
import { BookDetail, type Rect } from "./BookDetail";
import { BookSpine } from "./BookSpine";

const GAP = 2;
const LOOP_MIN_WIDTH = 2600;
const MAX_RY = 34;
const ARROW_STEP = 320;
const SHELVE_MS = 1100;

type ShelfProps = {
  books: Book[];
  /** A freshly shelved book id to glide in. */
  justAdded?: string | null;
  /** Let ← / → scroll this shelf when no book is open. */
  arrowKeys?: boolean;
};

type Slot = { el: HTMLElement; mid: number; ry: number };

// Edge fades as a mask, so they dissolve into whatever lies behind the shelf.
const FADE =
  "linear-gradient(90deg, transparent 0, #000 clamp(40px, 7vw, 110px), #000 calc(100% - clamp(40px, 7vw, 110px)), transparent 100%)";
const EDGE_FADE = { maskImage: FADE, WebkitMaskImage: FADE };

export function Shelf({ books, justAdded = null, arrowKeys = false }: ShelfProps) {
  const scrollerRef = useRef<HTMLDivElement>(null);
  const rowRef = useRef<HTMLDivElement>(null);
  const slots = useRef<Slot[]>([]);
  const frame = useRef(0);
  const drag = useRef({ active: false, moved: false, startX: 0, startScroll: 0 });
  const overflowingRef = useRef(true);
  const [overflowing, setOverflowing] = useState(true);
  const [shelving, setShelving] = useState<string | null>(null);
  const [open, setOpen] = useState<{ pos: number; rect: Rect } | null>(null);

  // One copy's width; a long row renders three copies and lives in the middle one.
  const segment = useMemo(() => books.reduce((sum, b) => sum + b.width + GAP, 0), [books]);
  const copies = segment > LOOP_MIN_WIDTH ? 3 : 1;
  const rendered = useMemo(
    () =>
      Array.from({ length: copies }, (_, c) => books.map((book, i) => ({ book, i, key: `${c}:${book.id}` }))).flat(),
    [books, copies],
  );

  // Each spine turns with its distance from the viewport centre; the middle stays flat.
  const updateCurve = useCallback(() => {
    const scroller = scrollerRef.current;
    if (!scroller) return;
    const half = scroller.clientWidth / 2 || 1;
    const center = scroller.scrollLeft + half;
    for (const slot of slots.current) {
      const t = Math.max(-1, Math.min(1, (slot.mid - center) / half));
      const ry = Math.round(-Math.sign(t) * MAX_RY * Math.abs(t) ** 1.35 * 100) / 100;
      if (ry !== slot.ry) {
        slot.ry = ry;
        slot.el.style.setProperty("--ry", `${ry}deg`);
      }
    }
  }, []);

  const scheduleCurve = useCallback(() => {
    cancelAnimationFrame(frame.current);
    frame.current = requestAnimationFrame(updateCurve);
  }, [updateCurve]);

  const measure = useCallback(() => {
    const scroller = scrollerRef.current;
    const row = rowRef.current;
    if (!scroller || !row) return;
    // Layout positions (unaffected by transforms) in the scroller's content coordinates.
    // The preserve-3d row is the spines' offsetParent in some engines, so add its offset.
    slots.current = Array.from(row.children as HTMLCollectionOf<HTMLElement>).map((el) => {
      const left = el.offsetParent === row ? row.offsetLeft + el.offsetLeft : el.offsetLeft;
      return { el, mid: left + el.offsetWidth / 2, ry: Number.NaN };
    });
    const over = row.scrollWidth > scroller.clientWidth + 1;
    overflowingRef.current = over;
    setOverflowing(over);
    updateCurve();
  }, [updateCurve]);

  useLayoutEffect(() => {
    const scroller = scrollerRef.current;
    if (scroller) scroller.scrollLeft = copies === 3 ? segment : 0;
    measure();
  }, [rendered, copies, segment, measure]);

  useEffect(() => {
    const scroller = scrollerRef.current;
    const row = rowRef.current;
    if (!scroller || !row) return;
    const ro = new ResizeObserver(() => measure());
    ro.observe(scroller);
    ro.observe(row);

    const onWheel = (e: WheelEvent) => {
      if (!overflowingRef.current || Math.abs(e.deltaY) <= Math.abs(e.deltaX)) return;
      scroller.scrollLeft += e.deltaY;
      e.preventDefault();
    };
    scroller.addEventListener("wheel", onWheel, { passive: false });
    return () => {
      ro.disconnect();
      scroller.removeEventListener("wheel", onWheel);
      cancelAnimationFrame(frame.current);
    };
  }, [measure]);

  const onScroll = () => {
    const scroller = scrollerRef.current;
    if (!scroller) return;
    if (copies === 3) {
      // Seamless loop: jump by one copy at the segment boundaries.
      const sl = scroller.scrollLeft;
      const shift = sl < segment * 0.5 ? segment : sl > segment * 1.5 ? -segment : 0;
      if (shift) {
        scroller.scrollLeft = sl + shift;
        drag.current.startScroll += shift;
      }
    }
    scheduleCurve();
  };

  useEffect(() => {
    if (!arrowKeys) return;
    const onKey = (e: KeyboardEvent) => {
      if (e.key !== "ArrowLeft" && e.key !== "ArrowRight") return;
      if (document.body.dataset.bookOpen) return;
      const target = e.target as HTMLElement | null;
      if (target?.closest("input, textarea, select, [contenteditable='true'], [role='dialog']")) return;
      e.preventDefault();
      scrollerRef.current?.scrollBy({ left: e.key === "ArrowRight" ? ARROW_STEP : -ARROW_STEP, behavior: "smooth" });
    };
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, [arrowKeys]);

  useEffect(() => {
    const row = rowRef.current;
    if (!justAdded || !row) return;
    const matches = row.querySelectorAll<HTMLElement>(`[data-book-id="${CSS.escape(justAdded)}"]`);
    const el = matches[Math.floor(matches.length / 2)];
    if (!el) return;
    setShelving(justAdded);
    el.scrollIntoView({ inline: "center", block: "nearest", behavior: "smooth" });
    const t = setTimeout(() => setShelving(null), SHELVE_MS);
    return () => clearTimeout(t);
  }, [justAdded]);

  const rectOf = (el: Element): Rect => {
    const r = el.getBoundingClientRect();
    return { left: r.left, top: r.top, width: r.width, height: r.height };
  };

  // Capture where the spine sits so the detail view can slide out from it.
  const openAt = (pos: number, button: HTMLElement) => setOpen({ pos, rect: rectOf(button) });

  const step = (dir: 1 | -1) => {
    if (!open) return;
    const copy = Math.floor(open.pos / books.length);
    const i = ((open.pos % books.length) + dir + books.length) % books.length;
    const pos = copy * books.length + i;
    const el = rowRef.current?.children[pos];
    if (!el) return;
    el.scrollIntoView({ inline: "center", block: "nearest", behavior: "instant" });
    setOpen({ pos, rect: rectOf(el) });
  };

  const openBook = open ? rendered[open.pos]?.book : undefined;

  return (
    <div className="relative">
      <div
        ref={scrollerRef}
        className={`no-scrollbar relative flex select-none overflow-x-auto ${overflowing ? "justify-start" : "justify-center"}`}
        style={{ perspective: "1400px", perspectiveOrigin: "50% 65%", ...(overflowing ? EDGE_FADE : null) }}
        onScroll={onScroll}
        onPointerDown={(e) => {
          if (e.pointerType !== "mouse" || e.button !== 0) return;
          drag.current = { active: true, moved: false, startX: e.clientX, startScroll: e.currentTarget.scrollLeft };
        }}
        onPointerMove={(e) => {
          const d = drag.current;
          if (!d.active) return;
          const dx = e.clientX - d.startX;
          if (Math.abs(dx) > 4) d.moved = true;
          if (d.moved) e.currentTarget.scrollLeft = d.startScroll - dx;
        }}
        onPointerUp={() => {
          drag.current.active = false;
          // Let this gesture's click see `moved`, then forget it.
          setTimeout(() => (drag.current.moved = false), 0);
        }}
        onPointerLeave={() => (drag.current.active = false)}
      >
        <div
          ref={rowRef}
          className={`flex w-max items-end gap-[2px] pb-6 pt-16 [transform-style:preserve-3d] ${copies === 1 ? "px-6 sm:px-12" : ""}`}
        >
          {rendered.map(({ book, key }, pos) => (
            <BookSpine
              key={key}
              book={book}
              pulled={open?.pos === pos}
              className={shelving === book.id ? "animate-shelve-in" : ""}
              onOpen={(button) => openAt(pos, button)}
              shouldIgnoreClick={() => drag.current.moved}
            />
          ))}
        </div>
      </div>

      <div
        aria-hidden
        className="pointer-events-none relative mx-auto -mt-6 h-px w-[94%] bg-gradient-to-r from-transparent via-foreground/25 to-transparent"
      />
      <div
        aria-hidden
        className="pointer-events-none mx-auto h-7 w-[90%] bg-gradient-to-b from-foreground/8 to-transparent blur-[3px]"
      />

      {open && openBook && (
        <BookDetail
          book={openBook}
          rect={open.rect}
          canStep={books.length > 1}
          onPrev={() => step(-1)}
          onNext={() => step(1)}
          onClose={() => setOpen(null)}
        />
      )}
    </div>
  );
}
