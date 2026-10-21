// Animation timings, all in one place (milliseconds). Lower = snappier.

/** A book sliding toward you when the pointer rests on it. */
export const HOVER_MS = 350;
/** The floating card that appears above a hovered book. */
export const HOVER_CARD_MS = 160;
/** A clicked book leaving the shelf and turning to show its cover (and going back). */
export const PULL_MS = 600;
/** After Escape / "Shelve it": when the detail view closes, while the book is nearly home (~70% of PULL_MS). */
export const RETRACT_CLOSE_MS = 420;
/** The blurred backdrop behind a pulled-out book. */
export const BACKDROP_MS = 450;
/** The details panel next to a pulled-out book: its delay and rise. */
export const PANEL_DELAY_MS = 150;
export const PANEL_RISE_MS = 600;
/** A recommended book gliding onto its shelf. */
export const SHELVE_MS = 800;
/** Gap between the page's parts fading in on load (counter, search, shelf). */
export const ENTRANCE_STEP_MS = 120;
