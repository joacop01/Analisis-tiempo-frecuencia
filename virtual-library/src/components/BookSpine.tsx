import { memo, useEffect, useRef, useState, type CSSProperties } from "react";
import { createPortal } from "react-dom";
import type { Book } from "../data/books";
import { COVER_W, faceFont } from "./bookFaces";
import { HOVER_CARD_MS, HOVER_MS } from "./motion";

const HOVER_PULL = 96;
const HOVER_LIFT = -26;
const LEAVE_GRACE_MS = 90;

/** The shelf only ever shows a cover as a thin strip or a foreshortened plane: use the medium size there. */
export const shelfCover = (url: string) => url.replace(/-L\.jpg$/, "-M.jpg");

// The constant surface layers (material texture, cylindrical sheen, wear, edges) are classes in styles.css;
// only per-book values stay inline.

function publisherMark(publisher: string): string {
  return publisher
    .replace(/[^\p{L}\s]/gu, " ")
    .split(/\s+/)
    .filter((w) => w.length > 2 && !/^(the|and|books|press|publishing|publishers|inc|ltd)$/i.test(w))
    .slice(0, 2)
    .map((w) => w[0]!.toUpperCase())
    .join("");
}

/** The front cover: real art, or a typeset card when there is none. */
export function CoverFace({
  book,
  src = book.cover,
  className = "",
}: {
  book: Book;
  src?: string;
  className?: string;
}) {
  if (book.cover) {
    return (
      <img
        src={src}
        alt=""
        draggable={false}
        decoding="async"
        className={`h-full w-full object-cover ${className}`}
        style={{ backgroundColor: book.spine }}
      />
    );
  }
  const accent = book.band ?? book.ink;
  return (
    <span
      className={`flex h-full w-full flex-col items-center justify-center gap-[8%] px-[12%] text-center ${faceFont[book.face]} ${className}`}
      style={{ backgroundColor: book.spine, color: book.ink }}
    >
      <span className="h-px w-1/3" style={{ backgroundColor: accent }} />
      <span className={`text-[22px] leading-tight ${book.caps ? "uppercase tracking-[0.12em] text-[17px]" : ""}`}>
        {book.title}
      </span>
      <span className="text-[13px] italic opacity-80">{book.author}</span>
      <span className="h-px w-1/3" style={{ backgroundColor: accent }} />
    </span>
  );
}

type BookBodyProps = {
  book: Book;
  showPages?: boolean;
  /** Front-cover art; defaults to the shelf size (the detail view passes the large one). */
  coverSrc?: string;
};

/** The book as a 3D object: spine face, hinged covers and the page block. */
export const BookBody = memo(function BookBody({ book, showPages = true, coverSrc }: BookBodyProps) {
  const { width: w, height: h } = book;
  const showAuthor = w >= 44;
  const showMark = w >= 30;
  const titleSize = Math.max(9, Math.min(showAuthor ? w * 0.34 : w * 0.46, 18));
  const rule = book.band ?? book.ink;
  const lettering: CSSProperties = { color: book.ink, opacity: 1 - book.wear * 0.45 };
  const mark = showMark ? publisherMark(book.publisher) : "";
  const art = book.cover ? shelfCover(book.cover) : "";

  return (
    <span className="absolute inset-0 block [transform-style:preserve-3d]">
      {/* Spine face */}
      <span className="absolute inset-0 block overflow-hidden rounded-[1.5px]" style={{ backgroundColor: book.spine }}>
        {art && !book.spineImage && (
          <span
            className="absolute inset-0"
            style={{
              backgroundImage: `url("${art}")`,
              backgroundSize: "auto 100%",
              backgroundPosition: "left center",
            }}
          />
        )}
        <span className="absolute inset-0" style={{ backgroundColor: book.spine, opacity: book.cover ? 0.78 : 1 }} />
        {book.spineImage ? (
          <img src={book.spineImage} alt="" draggable={false} className="absolute inset-0 h-full w-full object-cover" />
        ) : (
          <>
            <span className="absolute inset-x-0 top-[6%] h-[2px]" style={{ backgroundColor: rule, opacity: 0.8 }} />
            <span
              className="absolute inset-x-0 top-[calc(6%+5px)] h-px"
              style={{ backgroundColor: rule, opacity: 0.55 }}
            />
            <span className="absolute inset-x-0 bottom-[6%] h-[2px]" style={{ backgroundColor: rule, opacity: 0.8 }} />

            <span
              className="absolute inset-x-0 top-[11%] flex justify-center"
              style={{ bottom: showMark ? "17%" : "11%" }}
            >
              <span
                className={`flex h-full flex-col items-center justify-center ${faceFont[book.face]} ${
                  book.caps ? "uppercase tracking-[0.12em]" : ""
                }`}
                style={{ ...lettering, writingMode: "vertical-rl" }}
              >
                <span
                  className="max-h-full overflow-hidden text-ellipsis whitespace-nowrap leading-none"
                  style={{ fontSize: titleSize, fontWeight: book.face === "serif" ? 500 : 400 }}
                >
                  {book.title}
                </span>
                {showAuthor && (
                  <span
                    className="mr-[3px] max-h-full overflow-hidden text-ellipsis whitespace-nowrap leading-none opacity-80"
                    style={{ fontSize: Math.max(8, titleSize * 0.62) }}
                  >
                    {book.author}
                  </span>
                )}
              </span>
            </span>

            {mark && (
              <span
                className="absolute inset-x-0 bottom-[9%] text-center font-mono text-[7px] font-bold leading-none tracking-[0.08em]"
                style={lettering}
              >
                {mark}
              </span>
            )}
          </>
        )}

        {/* Material grain laid over the lettering, so the type reads as printed into it. */}
        <span className={`absolute inset-0 spine-texture-${book.finish}`} />
        <span className={`absolute inset-0 spine-sheen-${book.finish}`} />
        <span className="absolute inset-0 spine-wear" style={{ opacity: book.wear }} />
        <span className="absolute inset-0 rounded-[1.5px] spine-edge" />
      </span>

      {/* Front cover, hinged at the spine's right edge */}
      <span
        className="absolute left-full top-0 block overflow-hidden rounded-r-[2px] book-front"
        style={{ width: COVER_W, height: h }}
      >
        <CoverFace book={book} src={coverSrc ?? art} />
        <span className="absolute inset-0 book-front-shade" />
      </span>

      {/* Back cover */}
      <span
        className="absolute right-full top-0 block book-back"
        style={{ width: COVER_W, height: h, backgroundColor: book.spine }}
      />

      {/* Page block on top; hardcovers get a headband */}
      {showPages && (
        <span
          className="absolute bottom-full left-[1px] block book-pages"
          style={{ width: w - 2, height: COVER_W - 6 }}
        >
          {book.binding === "hardcover" && (
            <span
              className="absolute inset-x-0 bottom-0 h-[3px] book-headband"
              style={{ "--band": rule } as CSSProperties}
            />
          )}
        </span>
      )}
    </span>
  );
});

function Stars({ rating }: { rating: number }) {
  if (!rating) return <span>Unrated</span>;
  return (
    <span aria-label={`${rating} of 5 stars`} className="tracking-[0.15em] text-primary">
      {"★".repeat(rating)}
      <span className="text-muted-foreground/50">{"★".repeat(Math.max(0, 5 - rating))}</span>
    </span>
  );
}

/** Whether a spine's centre is within the window (the edge fades hide anything less). */
const onScreen = (r: DOMRect) => r.left + r.width / 2 > 0 && r.left + r.width / 2 < window.innerWidth;

function HoverCard({ book, anchor }: { book: Book; anchor: DOMRect }) {
  const half = 124;
  const left = Math.min(Math.max(anchor.left + anchor.width / 2, half + 8), window.innerWidth - half - 8);
  return createPortal(
    <div
      role="tooltip"
      className="pointer-events-none fixed z-[100] w-[248px] animate-fade-in rounded-xl border border-border/70 bg-card/90 px-4 py-3 shadow-[0_24px_60px_-28px_var(--shade)] backdrop-blur-md"
      style={{
        animationDuration: `${HOVER_CARD_MS}ms`,
        left,
        top: anchor.top - 14,
        transform: `translate(-50%, calc(-100% + ${HOVER_LIFT}px))`,
      }}
    >
      <p className="font-display text-[20px] leading-tight text-card-foreground">{book.title}</p>
      <p className="mt-1 text-[15px] text-muted-foreground">{book.author}</p>
      <p className="mt-2 text-[14px] capitalize text-muted-foreground">
        {[book.year || null, book.binding === "mass" ? "mass market" : book.binding].filter(Boolean).join(" · ")}
      </p>
      <p className="text-[14px] text-muted-foreground">
        <Stars rating={book.rating} />
      </p>
      {book.genres && book.genres.length > 0 && (
        <p className="mt-1 text-[14px] text-foreground/70">{book.genres.join(" · ")}</p>
      )}
    </div>,
    document.body,
  );
}

type BookSpineProps = {
  book: Book;
  /** Position in the shelf row, handed back to `onOpen`. */
  pos: number;
  /** False keeps only the sized, empty slot: the shelf builds the books near the viewport. */
  built?: boolean;
  pulled?: boolean;
  className?: string;
  onOpen: (pos: number, button: HTMLButtonElement) => void;
  shouldIgnoreClick?: () => boolean;
};

export const BookSpine = memo(function BookSpine({
  book,
  pos,
  built = true,
  pulled = false,
  className = "",
  onOpen,
  shouldIgnoreClick,
}: BookSpineProps) {
  const ref = useRef<HTMLButtonElement>(null);
  const leaveTimer = useRef<ReturnType<typeof setTimeout>>(undefined);
  // Whether the card is up for keyboard focus rather than for the pointer.
  const byKeyboard = useRef(false);
  const [anchor, setAnchor] = useState<DOMRect | null>(null);
  const hovered = anchor !== null && !pulled;

  useEffect(() => () => clearTimeout(leaveTimer.current), []);

  // No pointer moves a keyboard-focus card along: it follows its book as the shelf scrolls, and shows only while
  // the book is on screen.
  useEffect(() => {
    if (!hovered || !byKeyboard.current) return;
    const follow = () => ref.current && setAnchor(ref.current.getBoundingClientRect());
    window.addEventListener("scroll", follow, { capture: true, passive: true });
    return () => window.removeEventListener("scroll", follow, { capture: true });
  }, [hovered]);

  const enter = (keyboard: boolean) => {
    clearTimeout(leaveTimer.current);
    byKeyboard.current = keyboard;
    if (ref.current) setAnchor(ref.current.getBoundingClientRect());
  };
  const leave = () => {
    clearTimeout(leaveTimer.current);
    leaveTimer.current = setTimeout(() => setAnchor(null), LEAVE_GRACE_MS);
  };

  const lean = hovered ? 0 : book.lean;
  const pull = hovered ? HOVER_PULL : 0;
  const lift = hovered ? HOVER_LIFT : 0;

  return (
    // The button never transforms: it is the stable hit target while the book moves.
    <button
      ref={ref}
      type="button"
      data-book-id={book.id}
      data-hovered={hovered || undefined}
      aria-label={`${book.title} by ${book.author}`}
      className={`relative shrink-0 cursor-pointer rounded-[2px] outline-none [transform-style:preserve-3d] focus-visible:ring-2 focus-visible:ring-ring/70 ${className}`}
      style={
        {
          width: book.width,
          height: book.height,
          zIndex: hovered ? 40 : undefined,
          "--spine-w": `${book.width}px`,
        } as CSSProperties
      }
      onPointerEnter={(e) => e.pointerType === "mouse" && enter(false)}
      onPointerLeave={leave}
      // Only keyboard focus stands in for the pointer: a press does not mean the pointer is over this book, and no
      // pointerleave would ever take the card down. When the shelf hands focus to this book's twin as it loops,
      // the card carries over only if the twin had it up.
      onFocus={(e) => {
        const from = e.relatedTarget;
        const twin = from instanceof HTMLElement && from !== e.currentTarget && from.dataset.bookId === book.id;
        if (twin ? from.dataset.hovered : e.currentTarget.matches(":focus-visible")) enter(true);
      }}
      onBlur={leave}
      onClick={() => {
        if (shouldIgnoreClick?.()) return;
        setAnchor(null);
        if (ref.current) onOpen(pos, ref.current);
      }}
    >
      {/* Fixed size (not inset-0) so the book keeps its shape while the shelve-in gap opens.
          The shelf writes this span's rotateY (and hides it far off-screen) directly as it scrolls. */}
      <span
        className="absolute left-0 top-0 block [transform-style:preserve-3d]"
        style={{ width: book.width, height: book.height, transform: "rotateY(0deg)" }}
      >
        {built && (
          <span
            className="absolute inset-0 block [transform-style:preserve-3d]"
            style={{
              transform: `rotateZ(${lean}deg) translateZ(${pull + book.depth}px) translateY(${lift}px)`,
              transformOrigin: "0% 100%",
              transition: `transform ${HOVER_MS}ms cubic-bezier(0.22, 1, 0.36, 1)`,
              // visibility, not opacity: opacity < 1 would flatten the 3D book.
              visibility: pulled ? "hidden" : undefined,
            }}
          >
            <BookBody book={book} />
          </span>
        )}
      </span>
      {hovered && onScreen(anchor) && <HoverCard book={book} anchor={anchor} />}
    </button>
  );
});
