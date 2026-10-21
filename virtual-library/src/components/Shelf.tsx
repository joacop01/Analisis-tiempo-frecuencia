import {
  memo,
  startTransition,
  useCallback,
  useEffect,
  useLayoutEffect,
  useMemo,
  useRef,
  useState,
  type CSSProperties,
} from "react";
import { flushSync } from "react-dom";
import type { Book } from "../data/books";
import { BookDetail, type Rect } from "./BookDetail";
import { BookSpine, shelfCover } from "./BookSpine";
import { SHELVE_MS } from "./motion";

const GAP = 2;
const LOOP_MIN_WIDTH = 2600;
const MAX_RY = 34;
const ARROW_STEP = 320;
// Spines further than this outside the visible stretch are hidden: every live 3D book costs compositing work on
// each frame, on-screen or not.
const CULL_MARGIN = 300;
// A looping row builds only the books within this much more of the visible stretch, in every copy, so the jump
// between copies always lands on built books; the others are empty slots of the same size.
const BUILD_MARGIN = 500;
// Before hydration the viewport width is unknown: build a copy's start out to the right edge of a wide (QHD) screen,
// and the cull margin left of it.
const FIRST_VIEW = 2560;

type ShelfProps = {
  books: Book[];
  /** A freshly shelved book id to glide in. */
  justAdded?: string | null;
  /** Let ← / → scroll this shelf when no book is open. */
  arrowKeys?: boolean;
};

// Shelf-size covers already requested in this page (shared by every shelf).
const warmed = new Set<string>();

/** A spine's slot, in the scroller's content coordinates; `turn` is the span that carries its rotateY. */
type Slot = { el: HTMLElement; turn: HTMLElement; left: number; right: number; ry: number; hidden: boolean };
/**
 * The built stretch of a looping row, in one copy's coordinates (it may run past either end and wraps); `row` is the
 * rendered row it belongs to.
 */
type Built = { row: unknown; from: number; to: number };

// Edge fades as a mask, so they dissolve into whatever lies behind the shelf.
const FADE =
  "linear-gradient(90deg, transparent 0, #000 clamp(40px, 7vw, 110px), #000 calc(100% - clamp(40px, 7vw, 110px)), transparent 100%)";
const EDGE_FADE = { maskImage: FADE, WebkitMaskImage: FADE };

export const Shelf = memo(function Shelf({ books, justAdded = null, arrowKeys = false }: ShelfProps) {
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

  const firstBuilt = useMemo<Built>(() => ({ row: rendered, from: -CULL_MARGIN, to: FIRST_VIEW }), [rendered]);
  const [builtState, setBuilt] = useState(firstBuilt);
  // A new row (even the same books again, after a filter is cleared) starts over at its start, like its scroll.
  const built = builtState.row === rendered ? builtState : firstBuilt;
  // What is built on screen, and what is built or on its way (a deferred update may still be rendering).
  const shownRef = useRef(built);
  const builtRef = useRef(built);
  const loop = useRef({ segment, copies });

  // Which books (by index within a copy) are built; null when the row is short enough to build whole.
  const near = useMemo(() => {
    if (copies === 1 || built.to - built.from >= segment) return null;
    let x = 0;
    return books.map((b) => {
      const at = x;
      x += b.width + GAP;
      return [-segment, 0, segment].some((k) => at + k <= built.to && at + b.width + k >= built.from);
    });
  }, [books, copies, segment, built]);

  /**
   * Makes sure the books over the visible stretch [left, right] (content coordinates) and the cull margin around
   * it are built, in every copy. The books in view are built before the next paint, or the shelf would show empty
   * slots after a long jump or a fast fling; the margin can wait until React has time.
   */
  const keepBuilt = useCallback((left: number, right: number) => {
    const { segment: seg, copies: n } = loop.current;
    const cur = builtRef.current;
    if (n === 1 || cur.to - cur.from >= seg || !slots.current.length) return;
    const len = right - left;
    const at = (((left - slots.current[0]!.left) % seg) + seg) % seg;
    const covers = (w: Built, from: number, to: number) =>
      [-seg, 0, seg].some((k) => w.from <= at + from + k && at + to + k <= w.to);
    const inView = covers(shownRef.current, 0, len);
    const ahead = covers(cur, -CULL_MARGIN, len + CULL_MARGIN);
    if (inView && ahead) return;
    const next = ahead
      ? cur
      : { row: cur.row, from: at - CULL_MARGIN - BUILD_MARGIN, to: at + len + CULL_MARGIN + BUILD_MARGIN };
    builtRef.current = next;
    const later = () => startTransition(() => setBuilt(next));
    if (inView) return later();
    // Just the view first. flushSync is not allowed while React commits (measure() also runs from a layout effect);
    // a microtask still runs before the paint. The margin is queued after it, so it is the update that stays.
    const now = { row: cur.row, from: at, to: at + len };
    queueMicrotask(() => {
      flushSync(() => setBuilt(now));
      later();
    });
  }, []);

  useLayoutEffect(() => {
    shownRef.current = built;
    builtRef.current = built;
    loop.current = { segment, copies };
    // A deferred window can commit after the view has already moved on (an older update landing late): check
    // the view against what is now on screen.
    const scroller = scrollerRef.current;
    if (scroller) keepBuilt(scroller.scrollLeft, scroller.scrollLeft + scroller.clientWidth);
  }, [built, segment, copies, keepBuilt]);

  // Each spine turns with its distance from the viewport centre; the middle stays flat.
  const updateCurve = useCallback(() => {
    const scroller = scrollerRef.current;
    if (!scroller) return;
    const left = scroller.scrollLeft;
    const width = scroller.clientWidth;
    const half = width / 2 || 1;
    const center = left + half;
    const from = left - CULL_MARGIN;
    const to = left + width + CULL_MARGIN;
    const cull = overflowingRef.current;
    // Never the focused spine (hiding it would drop keyboard focus to the page) or its neighbours (Tab and
    // Shift+Tab skip hidden ones).
    const focused = document.activeElement;
    const keep =
      focused && focused.parentElement === rowRef.current
        ? [focused, focused.previousElementSibling, focused.nextElementSibling]
        : [];
    for (const slot of slots.current) {
      const t = Math.max(-1, Math.min(1, ((slot.left + slot.right) / 2 - center) / half));
      const ry = Math.round(-Math.sign(t) * MAX_RY * Math.abs(t) ** 1.35 * 100) / 100;
      // Written on the turning span itself: a custom property on the button would restyle the whole book.
      if (ry !== slot.ry) {
        slot.ry = ry;
        slot.turn.style.transform = `rotateY(${ry}deg)`;
      }
      // visibility, not display: the slot keeps its place in the row. It has to be the button itself: a visible
      // one still takes part in the 3D scene, and all of them are re-layerized on every frame.
      const hidden = cull && (slot.right < from || slot.left > to) && !keep.includes(slot.el);
      if (hidden !== slot.hidden) {
        slot.hidden = hidden;
        slot.el.style.visibility = hidden ? "hidden" : "";
      }
    }
    keepBuilt(left, left + width);
  }, [keepBuilt]);

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
    // Slots that survive keep what was last written to them, so only real changes are written again.
    const prev = new Map(slots.current.map((slot) => [slot.el, slot]));
    slots.current = Array.from(row.children as HTMLCollectionOf<HTMLElement>).map((el) => {
      const left = el.offsetParent === row ? row.offsetLeft + el.offsetLeft : el.offsetLeft;
      const turn = el.firstElementChild as HTMLElement;
      const old = prev.get(el);
      const kept = old?.turn === turn;
      return {
        el,
        turn,
        left,
        right: left + el.offsetWidth,
        ry: kept ? old.ry : Number.NaN,
        hidden: kept ? old.hidden : false,
      };
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

  // A looping row builds books only as they come near the view; fetch their covers in the background once the page
  // is idle, so they are in the cache by then instead of popping in after a long scroll.
  useEffect(() => {
    if (copies === 1) return;
    const urls = books.flatMap((b) => (b.cover && !warmed.has(shelfCover(b.cover)) ? [shelfCover(b.cover)] : []));
    if (!urls.length) return;
    const warm = () => {
      for (const url of urls) {
        if (warmed.has(url)) continue;
        warmed.add(url);
        const img = new Image();
        img.decoding = "async";
        img.fetchPriority = "low";
        img.src = url;
      }
    };
    if ("requestIdleCallback" in window) {
      const id = requestIdleCallback(warm, { timeout: 4000 });
      return () => cancelIdleCallback(id);
    }
    const t = setTimeout(warm, 1500);
    return () => clearTimeout(t);
  }, [books, copies]);

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
        // A focused spine is now a copy away: hand focus to its twin in view.
        const row = rowRef.current;
        const focused = document.activeElement;
        if (row && focused instanceof HTMLElement && focused.parentElement === row) {
          updateCurve();
          const pos = Array.prototype.indexOf.call(row.children, focused) + Math.sign(shift) * books.length;
          (row.children[pos] as HTMLElement | undefined)?.focus({ preventScroll: true });
        }
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

  // Capture where the spine sits so the detail view can slide out from it.
  const openAt = useCallback((pos: number, button: HTMLElement) => setOpen({ pos, rect: rectOf(button) }), []);
  const dragged = useCallback(() => drag.current.moved, []);

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
        className="no-scrollbar relative flex select-none overflow-x-auto"
        style={{ perspective: "1400px", perspectiveOrigin: "50% 65%", ...(overflowing ? EDGE_FADE : null) }}
        onScroll={onScroll}
        // The focused spine (and its neighbours) must not stay hidden.
        onFocus={scheduleCurve}
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
        {/* Auto margins centre a short row and give way once it overflows, so it is centred from the server render
            on, and spine positions never depend on whether it overflows. */}
        <div
          ref={rowRef}
          className={`mx-auto flex w-max items-end gap-[2px] pb-6 pt-16 [transform-style:preserve-3d] ${copies === 1 ? "px-6 sm:px-12" : ""}`}
          style={{ "--shelve-ms": `${SHELVE_MS}ms` } as CSSProperties}
        >
          {rendered.map(({ book, i, key }, pos) => (
            <BookSpine
              key={key}
              book={book}
              pos={pos}
              built={near ? near[i] : true}
              pulled={open?.pos === pos}
              className={shelving === book.id ? "animate-shelve-in" : ""}
              onOpen={openAt}
              shouldIgnoreClick={dragged}
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
});

function rectOf(el: Element): Rect {
  const r = el.getBoundingClientRect();
  return { left: r.left, top: r.top, width: r.width, height: r.height };
}
