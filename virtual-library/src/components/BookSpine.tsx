import { useEffect, useRef, useState, type CSSProperties } from "react";
import { createPortal } from "react-dom";
import type { Book } from "../data/books";
import { COVER_W, faceFont } from "./bookFaces";

const HOVER_PULL = 96;
const HOVER_LIFT = -26;
const LEAVE_GRACE_MS = 90;

const mix = (token: string, pct: number) => `color-mix(in oklch, var(${token}) ${pct}%, transparent)`;

// Cylindrical highlight across the spine, by surface material.
const finishSheen: Record<Book["finish"], string> = {
  gloss: `linear-gradient(90deg, ${mix("--shade", 42)} 0%, transparent 16%, ${mix("--light", 42)} 30%, ${mix("--light", 8)} 40%, transparent 58%, ${mix("--shade", 34)} 100%)`,
  matte: `linear-gradient(90deg, ${mix("--shade", 36)} 0%, transparent 20%, ${mix("--light", 14)} 34%, transparent 60%, ${mix("--shade", 30)} 100%)`,
  cloth: `linear-gradient(90deg, ${mix("--shade", 40)} 0%, transparent 22%, ${mix("--light", 9)} 38%, transparent 62%, ${mix("--shade", 36)} 100%)`,
};

const NOISE =
  "url(\"data:image/svg+xml;utf8,<svg xmlns='http://www.w3.org/2000/svg' width='120' height='120'><filter id='n'><feTurbulence type='fractalNoise' baseFrequency='1.1' numOctaves='2' stitchTiles='stitch'/><feColorMatrix type='saturate' values='0'/></filter><rect width='100%' height='100%' filter='url(%23n)'/></svg>\")";

// Material grain laid over the lettering, so the type reads as printed into it.
const finishTexture: Record<Book["finish"], CSSProperties> = {
  cloth: {
    backgroundImage: `repeating-linear-gradient(0deg, ${mix("--shade", 22)} 0 1px, transparent 1px 3px), repeating-linear-gradient(90deg, ${mix("--light", 14)} 0 1px, transparent 1px 3px)`,
    opacity: 0.55,
    mixBlendMode: "multiply",
  },
  matte: { backgroundImage: NOISE, opacity: 0.22, mixBlendMode: "multiply" },
  gloss: {
    backgroundImage: `linear-gradient(180deg, ${mix("--light", 10)}, transparent 30%, transparent 70%, ${mix("--light", 6)})`,
    opacity: 1,
  },
};

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
export function CoverFace({ book, className = "" }: { book: Book; className?: string }) {
  if (book.cover) {
    return (
      <img
        src={book.cover}
        alt=""
        draggable={false}
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

/** The book as a 3D object: spine face, hinged covers and the page block. */
export function BookBody({ book, showPages = true }: { book: Book; showPages?: boolean }) {
  const { width: w, height: h } = book;
  const showAuthor = w >= 44;
  const showMark = w >= 30;
  const titleSize = Math.max(9, Math.min(showAuthor ? w * 0.34 : w * 0.46, 18));
  const rule = book.band ?? book.ink;
  const lettering: CSSProperties = { color: book.ink, opacity: 1 - book.wear * 0.45 };
  const mark = showMark ? publisherMark(book.publisher) : "";

  return (
    <span className="absolute inset-0 block [transform-style:preserve-3d]">
      {/* Spine face */}
      <span className="absolute inset-0 block overflow-hidden rounded-[1.5px]" style={{ backgroundColor: book.spine }}>
        {book.cover && !book.spineImage && (
          <span
            className="absolute inset-0"
            style={{
              backgroundImage: `url("${book.cover}")`,
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

        <span className="absolute inset-0" style={finishTexture[book.finish]} />
        <span className="absolute inset-0" style={{ backgroundImage: finishSheen[book.finish] }} />
        <span
          className="absolute inset-0"
          style={{
            opacity: book.wear,
            backgroundImage: `linear-gradient(180deg, ${mix("--light", 55)} 0%, transparent 10%, transparent 90%, ${mix("--light", 45)} 100%), linear-gradient(90deg, ${mix("--light", 35)} 0%, transparent 7%, transparent 93%, ${mix("--light", 35)} 100%), linear-gradient(180deg, ${mix("--glow", 30)}, transparent 55%)`,
          }}
        />
        <span
          className="absolute inset-0 rounded-[1.5px]"
          style={{ boxShadow: `inset 1px 0 0 ${mix("--light", 18)}, inset -1px 0 0 ${mix("--shade", 30)}` }}
        />
      </span>

      {/* Front cover, hinged at the spine's right edge */}
      <span
        className="absolute left-full top-0 block overflow-hidden rounded-r-[2px]"
        style={{ width: COVER_W, height: h, transformOrigin: "left center", transform: "rotateY(90deg)" }}
      >
        <CoverFace book={book} />
        <span
          className="absolute inset-0"
          style={{
            backgroundImage: `linear-gradient(90deg, ${mix("--shade", 38)} 0%, ${mix("--light", 12)} 3%, transparent 9%, transparent 85%, ${mix("--shade", 12)} 100%)`,
          }}
        />
      </span>

      {/* Back cover */}
      <span
        className="absolute right-full top-0 block"
        style={{
          width: COVER_W,
          height: h,
          transformOrigin: "right center",
          transform: "rotateY(-90deg)",
          backgroundColor: book.spine,
          backgroundImage: `linear-gradient(270deg, ${mix("--shade", 45)} 0%, transparent 12%, ${mix("--shade", 25)} 100%)`,
        }}
      />

      {/* Page block on top; hardcovers get a headband */}
      {showPages && (
        <span
          className="absolute bottom-full left-[1px] block"
          style={{
            width: w - 2,
            height: COVER_W - 6,
            transformOrigin: "bottom center",
            transform: "rotateX(78deg)",
            backgroundColor: "var(--paper)",
            backgroundImage: `repeating-linear-gradient(90deg, ${mix("--shade", 7)} 0 1px, transparent 1px 2px), linear-gradient(0deg, ${mix("--shade", 18)}, transparent 40%)`,
          }}
        >
          {book.binding === "hardcover" && (
            <span
              className="absolute inset-x-0 bottom-0 h-[3px]"
              style={{
                backgroundImage: `repeating-linear-gradient(90deg, ${rule} 0 2px, var(--paper) 2px 3px)`,
              }}
            />
          )}
        </span>
      )}
    </span>
  );
}

function Stars({ rating }: { rating: number }) {
  if (!rating) return <span>Unrated</span>;
  return (
    <span aria-label={`${rating} of 5 stars`} className="tracking-[0.15em] text-primary">
      {"★".repeat(rating)}
      <span className="text-muted-foreground/50">{"★".repeat(Math.max(0, 5 - rating))}</span>
    </span>
  );
}

function HoverCard({ book, anchor }: { book: Book; anchor: DOMRect }) {
  const half = 124;
  const left = Math.min(Math.max(anchor.left + anchor.width / 2, half + 8), window.innerWidth - half - 8);
  return createPortal(
    <div
      role="tooltip"
      className="pointer-events-none fixed z-[100] w-[248px] animate-fade-in rounded-xl border border-border/70 bg-card/90 px-4 py-3 shadow-[0_24px_60px_-28px_var(--shade)] backdrop-blur-md [animation-duration:240ms]"
      style={{ left, top: anchor.top - 14, transform: `translate(-50%, calc(-100% + ${HOVER_LIFT}px))` }}
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
  pulled?: boolean;
  className?: string;
  onOpen: (button: HTMLButtonElement) => void;
  shouldIgnoreClick?: () => boolean;
};

export function BookSpine({ book, pulled = false, className = "", onOpen, shouldIgnoreClick }: BookSpineProps) {
  const ref = useRef<HTMLButtonElement>(null);
  const leaveTimer = useRef<ReturnType<typeof setTimeout>>(undefined);
  const [anchor, setAnchor] = useState<DOMRect | null>(null);
  const hovered = anchor !== null && !pulled;

  useEffect(() => () => clearTimeout(leaveTimer.current), []);

  const enter = () => {
    clearTimeout(leaveTimer.current);
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
      onPointerEnter={(e) => e.pointerType === "mouse" && enter()}
      onPointerLeave={leave}
      onFocus={enter}
      onBlur={leave}
      onClick={() => {
        if (shouldIgnoreClick?.()) return;
        setAnchor(null);
        if (ref.current) onOpen(ref.current);
      }}
    >
      {/* Fixed size (not inset-0) so the book keeps its shape while the shelve-in gap opens. */}
      <span
        className="absolute left-0 top-0 block [transform-style:preserve-3d]"
        style={{ width: book.width, height: book.height, transform: "rotateY(var(--ry, 0deg))" }}
      >
        <span
          className="absolute inset-0 block [transform-style:preserve-3d]"
          style={{
            transform: `rotateZ(${lean}deg) translateZ(${pull + book.depth}px) translateY(${lift}px)`,
            transformOrigin: "0% 100%",
            transition: "transform 700ms cubic-bezier(0.22, 1, 0.36, 1)",
            // visibility, not opacity: opacity < 1 would flatten the 3D book.
            visibility: pulled ? "hidden" : undefined,
          }}
        >
          <BookBody book={book} />
        </span>
      </span>
      {hovered && <HoverCard book={book} anchor={anchor} />}
    </button>
  );
}
