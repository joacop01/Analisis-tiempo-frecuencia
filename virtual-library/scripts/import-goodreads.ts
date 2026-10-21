// Regenerates src/data/books.ts from a Goodreads library export.
//
//   npm run import:goodreads -- [path/to/goodreads_library_export.csv] [--contact you@example.com] [--refresh]
//
// Reads only "read" / "currently-reading" rows, looks each book up on Open
// Library (cover, blurb, subjects), samples the cover's colours for the spine,
// derives the physical props deterministically, and overwrites the `books`
// array (the `Book` type above it is kept as is).
//
// Open Library allows 1 request/second from anonymous clients and 3/second from
// clients that put a contact email in their User-Agent, so `--contact` (or
// OPENLIBRARY_CONTACT=you@example.com) makes the import ~3x faster. Lookups and
// sampled cover colours are cached in .cache/openlibrary.json, so an interrupted
// or repeated run only fetches what it hasn't seen yet; `--refresh` ignores it.
import { existsSync, mkdirSync, readFileSync, renameSync, writeFileSync } from "node:fs";
import { dirname, resolve } from "node:path";
import { fileURLToPath } from "node:url";
import jpeg from "jpeg-js";
import type { Book } from "../src/data/books.ts";
import { derivePhysical, fallbackPalette, paletteFromPixels, slugify, type Palette } from "../src/lib/bookPhysics.ts";

const BOOKS_FILE = fileURLToPath(new URL("../src/data/books.ts", import.meta.url));
const CACHE_FILE = fileURLToPath(new URL("../.cache/openlibrary.json", import.meta.url));
const OL = "https://openlibrary.org";
const COVERS = "https://covers.openlibrary.org";
const MONTHS = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"];
const GENRES = ["Fiction", "Nonfiction", "Sci-Fi", "Mystery & Thriller", "Fantasy", "Romance"] as const;
const USAGE = "node scripts/import-goodreads.ts [goodreads_library_export.csv] [--contact you@example.com] [--refresh]";

const BATCH_SIZE = 50; // ISBNs per Books API request (~1000-char URLs)
const COVER_CONCURRENCY = 6; // covers fetched by id aren't rate limited
const MAX_ATTEMPTS = 6;
const TIMEOUT_MS = 30_000;

type Row = Record<string, string>;

// ── Options ─────────────────────────────────────────────────────────────────

function fail(message: string): never {
  console.error(message);
  process.exit(1);
}

function readOptions(argv: string[]) {
  let csv: string | undefined;
  let contact = process.env.OPENLIBRARY_CONTACT?.trim() ?? "";
  let refresh = false;
  for (let i = 0; i < argv.length; i++) {
    const arg = argv[i];
    if (arg === "--help" || arg === "-h") {
      console.log(`Usage: ${USAGE}`);
      process.exit(0);
    } else if (arg === "--refresh") refresh = true;
    else if (arg === "--contact" || arg.startsWith("--contact=")) {
      contact = (arg === "--contact" ? (argv[++i] ?? "") : arg.slice("--contact=".length)).trim();
      if (!contact) fail("--contact needs an email address, e.g. --contact you@example.com");
    } else if (arg.startsWith("-")) fail(`Unknown option ${arg}.\nUsage: ${USAGE}`);
    else csv ??= arg;
  }
  // It goes into the User-Agent header, which only takes plain ASCII.
  if (contact && !(/^[^\s@()<>]+@[^\s@()<>]+\.[^\s@()<>]+$/.test(contact) && /^[!-~]+$/.test(contact))) {
    fail(
      `The contact must be a plain email address (no spaces or accents), e.g. --contact you@example.com (got "${contact}").`,
    );
  }
  return { csv: resolve(csv ?? "goodreads_library_export.csv"), contact, refresh };
}

const options = readOptions(process.argv.slice(2));
// Open Library asks clients to identify themselves; with a contact email it allows 3 requests/second.
const HEADERS = { "User-Agent": `VirtualLibraryImport/1.0 (${options.contact || "personal bookshelf generator"})` };
const PER_SECOND = options.contact ? 3 : 1;

// ── CSV ─────────────────────────────────────────────────────────────────────

function parseCsv(text: string): Row[] {
  const rows: string[][] = [];
  let row: string[] = [];
  let field = "";
  let quoted = false;
  const src = text.replace(/^﻿/, "");
  for (let i = 0; i < src.length; i++) {
    const c = src[i];
    if (quoted) {
      if (c === '"' && src[i + 1] === '"') {
        field += '"';
        i++;
      } else if (c === '"') quoted = false;
      else field += c;
    } else if (c === '"') quoted = true;
    else if (c === ",") {
      row.push(field);
      field = "";
    } else if (c === "\n" || c === "\r") {
      if (c === "\r" && src[i + 1] === "\n") i++;
      row.push(field);
      rows.push(row);
      row = [];
      field = "";
    } else field += c;
  }
  if (field || row.length) {
    row.push(field);
    rows.push(row);
  }
  const [header, ...body] = rows;
  return body
    .filter((r) => r.some((v) => v.trim()))
    .map((r) => Object.fromEntries(header.map((h, i) => [h.trim(), (r[i] ?? "").trim()])));
}

// Goodreads writes ISBNs as ="0446310786".
const cleanIsbn = (v = "") => v.replace(/[="\s-]/g, "");
const ISBN = /^(\d{9}[\dXx]|\d{13})$/;
// ISBN13 first, then the ISBN10 if that one isn't on Open Library.
const isbnsOf = (row: Row) =>
  [cleanIsbn(row["ISBN13"]), cleanIsbn(row["ISBN"])].filter((v, i, all) => ISBN.test(v) && all.indexOf(v) === i);

function formatMonth(goodreadsDate: string): string {
  const m = goodreadsDate.match(/^(\d{4})\/(\d{1,2})/);
  return m ? `${MONTHS[Number(m[2]) - 1]} ${m[1]}` : "";
}

// "The Name of the Wind (The Kingkiller Chronicle, #1)" → "The Name of the Wind"
// "Title (Spanish Edition)" / "(Edición Kindle)" → "Title" (repeatedly: series and edition can both be there).
const TRAILING_SERIES_OR_EDITION = /\s*\((?:[^()]*#\d+(?:\.\d+)?|[^()]*\b(?:edition|edici[oó]n)\b[^()]*)\)\s*$/i;
const stripSeries = (title: string) => {
  let t = title.trim();
  for (let prev = ""; prev !== t; ) [prev, t] = [t, t.replace(TRAILING_SERIES_OR_EDITION, "").trim()];
  return t;
};

// ── Cache ───────────────────────────────────────────────────────────────────
// API responses keyed by URL (trimmed to the fields used here; null = not on
// Open Library) and sampled cover palettes keyed by cover id. Failures are never
// cached, so re-running retries exactly what's missing.

type CacheFile = { version: 1; api: Record<string, unknown>; palettes: Record<string, Palette | null> };

const cache = loadCache();
const fetchedThisRun = new Set<string>();
let cacheDirty = false;
let cacheTimer: ReturnType<typeof setTimeout> | undefined;
let cacheWarned = false;

function loadCache(): CacheFile {
  const empty: CacheFile = { version: 1, api: {}, palettes: {} };
  if (!existsSync(CACHE_FILE)) return empty;
  try {
    const parsed: unknown = JSON.parse(readFileSync(CACHE_FILE, "utf8"));
    if (isRecord(parsed) && parsed.version === 1 && isRecord(parsed.api) && isRecord(parsed.palettes)) {
      return parsed as CacheFile;
    }
  } catch {
    // fall through: start over
  }
  console.warn(`Ignoring unreadable cache ${CACHE_FILE}.`);
  return empty;
}

/** The cached response for `url` (null = known to be missing), or undefined if it has to be fetched. */
function cached(url: string): unknown {
  if (options.refresh && !fetchedThisRun.has(url)) return undefined;
  return Object.hasOwn(cache.api, url) ? cache.api[url] : undefined;
}

let keptStale = 0;

/** With --refresh, a lookup that fails keeps what an earlier run cached instead of losing it. */
function orCached<T>(got: Got<T>, store: Record<string, unknown>, key: string): Got<T> {
  if (got.ok || !options.refresh || !Object.hasOwn(store, key)) return got;
  keptStale++;
  return { ok: true, value: store[key] as T | null };
}

function remember(url: string, value: unknown) {
  cache.api[url] = value;
  fetchedThisRun.add(url);
  markDirty();
}

function markDirty() {
  cacheDirty = true;
  cacheTimer ??= setTimeout(() => {
    cacheTimer = undefined;
    saveCache();
  }, 2000).unref();
}

function saveCache() {
  if (!cacheDirty) return;
  cacheDirty = false;
  try {
    mkdirSync(dirname(CACHE_FILE), { recursive: true });
    writeFileSync(`${CACHE_FILE}.tmp`, JSON.stringify(cache));
    renameSync(`${CACHE_FILE}.tmp`, CACHE_FILE);
  } catch (err) {
    cacheDirty = true; // try again on the next save (e.g. the file was briefly locked on Windows)
    if (!cacheWarned) console.warn(`Couldn't save the cache to ${CACHE_FILE} (${describe(err)}).`);
    cacheWarned = true;
  }
}

// ── Requests ────────────────────────────────────────────────────────────────

/** ok + null = Open Library doesn't have it; !ok = it failed even after retries. */
type Got<T> = { ok: true; value: T | null } | { ok: false; reason: string };
type Limiter = ReturnType<typeof makeLimiter>;

// Milliseconds between API requests: 10% under the documented limit.
const spacing = (perSecond: number) => Math.ceil((1000 / perSecond) * 1.1);

/**
 * Hands out evenly spaced turns for every openlibrary.org API request. Retries go
 * first, then the lowest `order` (the book's place on the shelf), so books finish
 * steadily from the top of the shelf down.
 */
function makeLimiter(perSecond: number) {
  const base = spacing(perSecond);
  let interval = base;
  let calm = 0; // successes since the last 429
  let nextAt = 0;
  let timer: ReturnType<typeof setTimeout> | undefined;
  let seq = 0;
  let released = false;
  const waiting: { order: number; seq: number; go: () => void }[] = [];
  const pump = () => {
    if (timer || !waiting.length) return;
    const wait = nextAt - Date.now();
    if (wait > 0) {
      timer = setTimeout(() => {
        timer = undefined;
        pump();
      }, wait);
      return;
    }
    nextAt = Date.now() + interval;
    limiter.started++;
    let next = 0;
    for (let i = 1; i < waiting.length; i++) {
      const a = waiting[i];
      const b = waiting[next];
      if (a.order < b.order || (a.order === b.order && a.seq < b.seq)) next = i;
    }
    waiting.splice(next, 1)[0].go();
    pump();
  };
  const limiter = {
    started: 0,
    get queued() {
      return waiting.length;
    },
    get interval() {
      return interval;
    },
    turn(order: number) {
      if (released) return Promise.resolve();
      return new Promise<void>((go) => {
        waiting.push({ order, seq: seq++, go });
        pump();
      });
    },
    /** Lets everyone through at once (the service is gone, so nobody will actually send anything). */
    release() {
      released = true;
      clearTimeout(timer);
      timer = undefined;
      for (const w of waiting.splice(0)) w.go();
    },
    /** The server asked us to back off: pause everyone, and (on a 429) go a bit slower for a while. */
    holdOff(ms: number, slower: boolean) {
      nextAt = Math.max(nextAt, Date.now() + ms);
      if (slower) {
        interval = Math.min(Math.round(interval * 1.2), base * 2);
        calm = 0;
      }
    },
    succeeded() {
      if (interval > base && ++calm >= 30) {
        interval = Math.max(base, Math.round(interval / 1.2));
        calm = 0;
      }
    },
  };
  return limiter;
}

/**
 * Notices when a service goes away mid-run: once 5 requests in a row have failed after every
 * retry and nothing has answered for a minute, the rest fail straight away (they're listed at the
 * end and a re-run fills them in) instead of each one retrying for half a minute.
 */
function makeBreaker(name: string, what: string, onOpen = () => {}) {
  let lastAnswer = Date.now();
  let failedInARow = 0;
  const breaker = {
    open: "",
    answered() {
      lastAnswer = Date.now();
      failedInARow = 0;
    },
    failed(reason: string) {
      const silent = Date.now() - lastAnswer;
      if (breaker.open || ++failedInARow < 5 || silent < 60_000) return;
      breaker.open = `skipped: ${name} stopped answering`;
      console.log(
        `\n${name} hasn't answered for ${formatDuration(silent)} (${reason}), so the remaining ${what} ` +
          "are skipped. Run the same command again later to fill them in.\n",
      );
      onOpen();
    },
  };
  return breaker;
}

const api = makeLimiter(PER_SECOND);
const apiBreaker = makeBreaker("Open Library", "lookups", () => api.release());
const coverBreaker = makeBreaker("The Open Library cover server", "cover images");
let apiReached = false;
let networkErrors = 0;
let refusals = 0;
let lastNotice = 0;

const sleep = (ms: number) => new Promise((r) => setTimeout(r, ms));
// 1, 2, 4, 8, 16 s (±20%) between attempts
const backoff = (attempt: number) => Math.min(30_000, 500 * 2 ** attempt) * (0.8 + Math.random() * 0.4);

function describe(err: unknown): string {
  if (!(err instanceof Error)) return String(err);
  return err.cause instanceof Error ? `${err.message}: ${err.cause.message}` : err.message;
}

function retryAfterMs(res: Response): number | undefined {
  const header = res.headers.get("retry-after");
  if (!header) return undefined;
  const seconds = Number(header);
  const ms = Number.isFinite(seconds) ? seconds * 1000 : Date.parse(header) - Date.now();
  return Number.isFinite(ms) ? Math.max(0, ms) : undefined;
}

function giveUp(problem: string): never {
  saveCache();
  console.error(
    `\n${problem}\nsrc/data/books.ts was not changed. Lookups that already finished are cached, so the next run picks up from here.`,
  );
  process.exit(1);
}

/**
 * GET with retries. API requests (`limiter` given) wait for their turn on every
 * attempt; 429/503 + Retry-After pause the whole limiter. `target` may be a
 * function so a batch can keep collecting ISBNs until its turn comes.
 */
async function request<T>(
  target: string | (() => string),
  read: (res: Response) => Promise<T | null>,
  limiter?: Limiter,
  order = 0,
): Promise<Got<T>> {
  let url = "";
  let reason = "";
  const breaker = limiter ? apiBreaker : coverBreaker;
  for (let attempt = 1; ; attempt++) {
    if (breaker.open) return { ok: false, reason: breaker.open };
    if (limiter) await limiter.turn(attempt > 1 ? -1 : order);
    if (breaker.open) return { ok: false, reason: breaker.open };
    url ||= typeof target === "string" ? target : target();
    let wait = backoff(attempt);
    let held = false;
    const init = { headers: HEADERS, signal: AbortSignal.timeout(TIMEOUT_MS) };
    const res = await fetch(url, init).catch((err: unknown) => {
      reason = describe(err);
      return null;
    });
    if (!res) {
      if (limiter && !apiReached && ++networkErrors >= 8) {
        giveUp(`Can't reach Open Library (${reason}). Check your internet connection and run the command again.`);
      }
    } else if (res.status === 404 || res.status === 410) {
      await res.body?.cancel().catch(() => {});
      if (limiter) apiReached = true;
      breaker.answered();
      return { ok: true, value: null };
    } else if (res.ok) {
      try {
        const value = await read(res);
        if (limiter) {
          apiReached = true;
          limiter.succeeded();
        }
        breaker.answered();
        return { ok: true, value };
      } catch (err) {
        reason = `unexpected response (${describe(err)})`;
        await res.body?.cancel().catch(() => {}); // free the connection if `read` stopped early
      }
    } else {
      await res.body?.cancel().catch(() => {});
      reason = `HTTP ${res.status}`;
      if (!(res.status === 408 || res.status === 425 || res.status === 429 || res.status >= 500)) break;
      const hinted = retryAfterMs(res);
      wait = Math.min(hinted ?? wait, 120_000);
      // Too many requests (or told when to come back): every API request waits, not just this one.
      if (limiter && (res.status === 429 || hinted !== undefined)) {
        limiter.holdOff(wait, res.status === 429);
        held = true;
      }
    }
    if (attempt >= MAX_ATTEMPTS) break;
    if (wait >= 5000 && Date.now() - lastNotice > 10_000) {
      lastNotice = Date.now();
      console.log(`  … an Open Library request failed (${reason}); retrying in ${Math.round(wait / 1000)} s.`);
    }
    if (!held) await sleep(wait);
  }
  if (limiter && !apiReached && ++refusals >= 5) {
    giveUp(`Open Library keeps refusing requests (${reason}). Wait a while and run the command again.`);
  }
  breaker.failed(reason);
  return { ok: false, reason };
}

function isRecord(v: unknown): v is Record<string, unknown> {
  return typeof v === "object" && v !== null && !Array.isArray(v);
}

/** Parses a JSON body and keeps only what `trim` returns (throwing = malformed, retried). */
const readJson = (trim: (json: unknown) => unknown) => async (res: Response) => trim(await res.json());

const apiCalls = new Map<string, Promise<Got<unknown>>>();

/** A cached, de-duplicated, rate-limited openlibrary.org GET (`order`: see makeLimiter). */
function getApi(url: string, trim: (json: unknown) => unknown, order: number): Promise<Got<unknown>> {
  const hit = cached(url);
  if (hit !== undefined) return Promise.resolve({ ok: true, value: hit });
  let call = apiCalls.get(url);
  if (!call) {
    call = request(url, readJson(trim), api, order).then((got) => {
      if (got.ok) remember(url, got.value);
      return orCached(got, cache.api, url);
    });
    apiCalls.set(url, call);
  }
  return call;
}

// ── Open Library ────────────────────────────────────────────────────────────

type Edition = { workKey?: string; coverId?: number; pages?: number; publisher?: string };
type Work = { description?: string; subjects: string[]; coverId?: number };
type SearchDoc = {
  key?: string;
  cover_i?: number;
  first_publish_year?: number;
  number_of_pages_median?: number;
  publisher?: string;
  subject: string[];
};

type Found = {
  workKey?: string;
  coverId?: number;
  pages?: number;
  publisher?: string;
  firstYear?: number;
  description?: string;
  subjects: string[];
};

const firstCover = (v: unknown) =>
  Array.isArray(v) ? v.find((c): c is number => typeof c === "number" && c > 0) : undefined;
const num = (v: unknown) => (typeof v === "number" && Number.isFinite(v) ? v : undefined);
const strings = (v: unknown) => (Array.isArray(v) ? v.filter((s): s is string => typeof s === "string") : []);
function firstName(v: unknown): string | undefined {
  const first: unknown = Array.isArray(v) ? v[0] : undefined;
  if (typeof first === "string") return first;
  return isRecord(first) && typeof first.name === "string" ? first.name : undefined;
}

// Edition records come from the Books API, many ISBNs per request:
//   /api/books?bibkeys=ISBN:a,ISBN:b&jscmd=details&format=json → { "ISBN:a": { details: <edition> }, … }
// Each ISBN is cached under its own single-ISBN URL so later runs can batch differently.
const booksApiUrl = (isbns: string[]) =>
  `${OL}/api/books?bibkeys=${isbns.map((isbn) => `ISBN:${isbn}`).join(",")}&jscmd=details&format=json`;

function trimEdition(details: unknown) {
  if (!isRecord(details)) return null;
  const { covers, works, number_of_pages, publishers } = details;
  return { covers, works, number_of_pages, publishers };
}

function toEdition(raw: unknown): Edition | null {
  if (!isRecord(raw)) return null;
  const work: unknown = Array.isArray(raw.works) ? raw.works[0] : undefined;
  return {
    workKey: isRecord(work) && typeof work.key === "string" ? work.key : undefined,
    coverId: firstCover(raw.covers),
    pages: num(raw.number_of_pages),
    publisher: firstName(raw.publishers),
  };
}

type Batch = { isbns: string[]; order: number; settle: Map<string, (got: Got<unknown>) => void> };
let openBatch: Batch | undefined;
const editions = new Map<string, Promise<Got<Edition>>>();

function editionFor(isbn: string, order: number): Promise<Got<Edition>> {
  let edition = editions.get(isbn);
  if (!edition) {
    edition = fetchEdition(isbn, order);
    editions.set(isbn, edition);
  }
  return edition;
}

async function fetchEdition(isbn: string, order: number): Promise<Got<Edition>> {
  const hit = cached(booksApiUrl([isbn]));
  if (hit !== undefined) return { ok: true, value: toEdition(hit) };
  const got = await new Promise<Got<unknown>>((settle) => {
    if (!openBatch || openBatch.isbns.length >= BATCH_SIZE) {
      openBatch = { isbns: [], order, settle: new Map() };
      void sendBatch(openBatch);
    }
    openBatch.isbns.push(isbn);
    openBatch.settle.set(isbn, settle);
  });
  return got.ok ? { ok: true, value: toEdition(got.value) } : got;
}

async function sendBatch(batch: Batch) {
  // The batch keeps collecting ISBNs until the limiter gives it a turn.
  const close = () => {
    if (openBatch === batch) openBatch = undefined;
    return booksApiUrl(batch.isbns);
  };
  const got = await request(
    close,
    readJson((json) => {
      if (!isRecord(json)) throw new Error("expected an object of bibkeys");
      return json;
    }),
    api,
    batch.order,
  );
  // A request skipped without sending (Open Library stopped answering) never closed the batch: close it
  // now, or ISBNs added later would wait for an answer that never comes.
  close();
  const byKey = new Map<string, unknown>();
  if (got.ok && got.value) {
    for (const [key, entry] of Object.entries(got.value)) {
      if (isRecord(entry)) byKey.set(key.toUpperCase(), entry.details);
    }
  }
  for (const isbn of batch.isbns) {
    const settle = batch.settle.get(isbn);
    if (!got.ok || !got.value) {
      settle?.(orCached(got.ok ? { ok: false, reason: "Books API not found" } : got, cache.api, booksApiUrl([isbn])));
      continue;
    }
    const edition = trimEdition(byKey.get(`ISBN:${isbn}`.toUpperCase()));
    remember(booksApiUrl([isbn]), edition);
    settle?.({ ok: true, value: edition });
  }
}

const SEARCH_FIELDS = "key,cover_i,first_publish_year,number_of_pages_median,publisher,subject";

function searchUrl(row: Row): string {
  const params = new URLSearchParams({
    title: stripSeries(row["Title"]),
    author: row["Author"],
    limit: "1",
    fields: SEARCH_FIELDS,
  });
  return `${OL}/search.json?${params}`;
}

function trimSearch(json: unknown) {
  if (!isRecord(json) || !Array.isArray(json.docs)) throw new Error("expected search docs");
  return { docs: json.docs.slice(0, 1) };
}

function toDoc(raw: unknown): SearchDoc | null {
  const doc: unknown = isRecord(raw) && Array.isArray(raw.docs) ? raw.docs[0] : undefined;
  if (!isRecord(doc)) return null;
  return {
    key: typeof doc.key === "string" ? doc.key : undefined,
    cover_i: num(doc.cover_i),
    first_publish_year: num(doc.first_publish_year),
    number_of_pages_median: num(doc.number_of_pages_median),
    publisher: firstName(doc.publisher),
    subject: strings(doc.subject),
  };
}

const workUrl = (key: string) => `${OL}${key}.json`;

function trimWork(json: unknown) {
  // Every Open Library record has a key; anything else (e.g. `{}` from a struggling server) isn't one to cache.
  if (!isRecord(json) || typeof json.key !== "string") throw new Error("expected a work record");
  const { description, subjects, covers } = json;
  return { description, subjects, covers };
}

function toWork(raw: unknown): Work | null {
  if (!isRecord(raw)) return null;
  const d = raw.description;
  return {
    description: typeof d === "string" ? d : isRecord(d) && typeof d.value === "string" ? d.value : undefined,
    subjects: strings(raw.subjects),
    coverId: firstCover(raw.covers),
  };
}

function workFor(key: string, order: number): Promise<Got<unknown>> {
  if (!/^\/works\/OL\d+W$/.test(key)) return Promise.resolve({ ok: true, value: null });
  return getApi(workUrl(key), trimWork, order);
}

async function lookUp(row: Row, order: number, failures: string[]): Promise<{ found: Found; matched: boolean }> {
  const value = (got: Got<unknown>, what: string) => {
    if (got.ok) return got.value;
    failures.push(`${what}: ${got.reason}`);
    return null;
  };

  let edition: Edition | null = null;
  let editionFailure: Got<unknown> | undefined;
  for (const isbn of isbnsOf(row)) {
    const got = await editionFor(isbn, order);
    if (!got.ok) editionFailure ??= got;
    else if ((edition = got.value)) break;
  }
  if (!edition && editionFailure) value(editionFailure, "edition");
  // The edition's own cover always wins, so start sampling it straight away.
  if (edition?.coverId) void paletteFor(edition.coverId);

  // A title search used to run for every book; now it only runs when it can still add
  // something (no edition/work, or the CSV and edition lack the year, pages or publisher).
  const needSearch =
    !edition?.workKey ||
    !(Number(row["Original Publication Year"]) || Number(row["Year Published"])) ||
    (!Number(row["Number of Pages"]) && edition.pages === undefined) ||
    (!row["Publisher"] && edition.publisher === undefined);
  const search = needSearch ? getApi(searchUrl(row), trimSearch, order) : undefined;
  const editionWork = edition?.workKey ? workFor(edition.workKey, order) : undefined;

  let doc = search ? toDoc(value(await search, "search")) : null;
  const workKey = edition?.workKey ?? doc?.key;
  const work = workKey ? toWork(value(await (editionWork ?? workFor(workKey, order)), "work")) : null;
  // Same precedence as ever (edition → search → work): search for a cover only if neither has one.
  if (!search && !edition?.coverId && !work?.coverId) {
    doc = toDoc(value(await getApi(searchUrl(row), trimSearch, order), "search"));
  }

  return { found: merge(edition, doc, work), matched: Boolean(edition || doc) };
}

function merge(edition: Edition | null, doc: SearchDoc | null, work: Work | null): Found {
  const found: Found = { subjects: [] };
  if (edition) {
    found.workKey = edition.workKey;
    found.coverId = edition.coverId;
    found.pages = edition.pages;
    found.publisher = edition.publisher;
  }
  if (doc) {
    found.workKey ??= doc.key;
    found.coverId ??= doc.cover_i;
    found.pages ??= doc.number_of_pages_median;
    found.publisher ??= doc.publisher;
    found.firstYear = doc.first_publish_year;
    found.subjects.push(...doc.subject);
  }
  if (work) {
    found.description = work.description;
    found.subjects.push(...work.subjects);
    found.coverId ??= work.coverId;
  }
  return found;
}

// Cover images run outside the API limiter, a few at a time.
function makePool(size: number) {
  let active = 0;
  const queue: (() => void)[] = [];
  return async <T>(task: () => Promise<T>): Promise<T> => {
    if (active < size) active++;
    else await new Promise<void>((go) => queue.push(go)); // the finishing task hands over its slot
    try {
      return await task();
    } finally {
      const next = queue.shift();
      if (next) next();
      else active--;
    }
  };
}

const coverSlot = makePool(COVER_CONCURRENCY);
const palettes = new Map<number, Promise<Got<Palette>>>();

function paletteFor(coverId: number): Promise<Got<Palette>> {
  let palette = palettes.get(coverId);
  if (!palette) {
    palette = samplePalette(coverId);
    palettes.set(coverId, palette);
  }
  return palette;
}

async function samplePalette(coverId: number): Promise<Got<Palette>> {
  const key = String(coverId);
  if (!options.refresh && Object.hasOwn(cache.palettes, key)) return { ok: true, value: cache.palettes[key] };
  // ?default=false → 404 instead of a placeholder image when there's no cover.
  const got = await coverSlot(() =>
    request(`${COVERS}/b/id/${coverId}-M.jpg?default=false`, async (res) => {
      if ((res.headers.get("content-type") ?? "").startsWith("text/")) throw new Error("not an image");
      const bytes = new Uint8Array(await res.arrayBuffer());
      try {
        const img = jpeg.decode(bytes, { useTArray: true, formatAsRGBA: true });
        return paletteFromPixels(img.data, img.width, img.height);
      } catch {
        return null; // unreadable image → plain spine, as if there were no cover
      }
    }),
  );
  if (got.ok) {
    cache.palettes[key] = got.value;
    markDirty();
  }
  return orCached(got, cache.palettes, key);
}

// ── Shaping ─────────────────────────────────────────────────────────────────

function cleanBlurb(raw = ""): string {
  const text = raw
    .split(/\n-{3,}/)[0]
    .replace(/\(\[source\]\[\d+\]\)/gi, "")
    .replace(/\[([^\]]+)\]\([^)]*\)/g, "$1")
    .replace(/\[([^\]]+)\]\[\d+\]/g, "$1")
    .replace(/^\s*\[\d+\]:.*$/gm, "")
    .replace(/\s+/g, " ")
    .trim();
  if (text.length <= 420) return text;
  const cut = text.slice(0, 420);
  const end = cut.lastIndexOf(". ");
  return end > 200 ? cut.slice(0, end + 1) : `${cut.replace(/\s+\S*$/, "")}…`;
}

// Matched against shelves and subjects lower-cased, without accents, and with hyphens and underscores as spaces, so
// Goodreads shelves ("ciencia-ficcion") and English or Spanish Open Library subjects ("Ciencia ficción") read alike.
const GENRE_PATTERNS: [(typeof GENRES)[number], RegExp][] = [
  [
    "Mystery & Thriller",
    /mystery|thriller|detective|crime|suspense|murder|misterio|policia(l|ca|co)|suspenso|novela negra/,
  ],
  ["Sci-Fi", /science fiction|sci ?fi|space opera|dystopia|cyberpunk|ciencia ficcion|ficcion cientifica|distopia/],
  [
    "Fantasy",
    /fantasy|magic|dragons|fae\b|faerie|fairies|wizards|witches|fantasia|literatura fantastica|dragones|hadas/,
  ],
  ["Romance", /romance|love stories|romantic|novela rosa|historias de amor/],
];
const NONFICTION =
  /non ?fiction|biography|memoir|history|essays|self help|philosophy|science\b|no ficcion|biografia|memorias|historia|ensayo|autoayuda|filosofia|divulgacion|ciencia\b/;
const FICTION = /(?<!non ?)fiction|novel|literary|(?<!no )ficcion|novela|narrativa|literatura|cuentos/;

function assignGenres(shelves: string, subjects: string[]): string[] {
  const text = [shelves, ...subjects]
    .join(" | ")
    .toLowerCase()
    .normalize("NFD")
    .replace(/[\u0300-\u036f]/g, "")
    .replace(/[-_]/g, " ");
  const genres = GENRE_PATTERNS.filter(([, re]) => re.test(text)).map(([genre]) => genre);
  if (genres.length) return genres.slice(0, 3);
  if (NONFICTION.test(text)) return ["Nonfiction"];
  // "Literary Fiction" and friends → Fiction; anything untagged defaults to Nonfiction.
  return FICTION.test(text) ? ["Fiction"] : ["Nonfiction"];
}

function bookLiteral(book: Book): string {
  const lines = Object.entries(book)
    .filter(([, v]) => v !== undefined)
    .map(([k, v]) => `    ${k}: ${JSON.stringify(v)},`);
  return `  {\n${lines.join("\n")}\n  }`;
}

type Outcome = { book: Book; matched: boolean; hasCover: boolean; failures: string[] };

async function importRow(row: Row, i: number): Promise<Outcome> {
  const failures: string[] = [];
  const { found, matched } = await lookUp(row, i, failures);
  const key = found.workKey ?? (cleanIsbn(row["ISBN13"]) || `goodreads:${row["Book Id"]}`);
  let sampled: Palette | null = null;
  if (found.coverId) {
    const got = await paletteFor(found.coverId);
    if (got.ok) sampled = got.value;
    else failures.push(`cover image: ${got.reason}`);
  }
  const palette = sampled || fallbackPalette(key);
  const pages = Number(row["Number of Pages"]) || found.pages;
  const title = stripSeries(row["Title"]);

  const book: Book = {
    id: `${slugify(title).slice(0, 28)}-${i}`,
    title,
    author: row["Author"],
    genres: assignGenres(row["Bookshelves"] ?? "", found.subjects),
    cover: found.coverId ? `https://covers.openlibrary.org/b/id/${found.coverId}-L.jpg` : "",
    year: Number(row["Original Publication Year"]) || Number(row["Year Published"]) || found.firstYear || 0,
    blurb: cleanBlurb(found.description),
    rating: Number(row["My Rating"]) || 0,
    finished:
      formatMonth(row["Date Read"]) || (row["Exclusive Shelf"] === "currently-reading" ? "Currently reading" : ""),
    publisher: row["Publisher"] || found.publisher || "",
    ...palette,
    ...derivePhysical(key, pages),
  };
  return { book, matched, hasCover: Boolean(found.coverId), failures };
}

// ── Progress ────────────────────────────────────────────────────────────────

const roughly = (ms: number) =>
  ms < 90_000 ? `${Math.max(10, Math.round(ms / 10_000) * 10)} s` : `${Math.round(ms / 60_000)} min`;

function formatDuration(ms: number): string {
  const s = Math.max(1, Math.round(ms / 1000));
  if (s < 60) return `${s} s`;
  if (s < 600) return `${Math.floor(s / 60)} min ${String(s % 60).padStart(2, "0")} s`;
  return `${Math.round(s / 60)} min`;
}

/** Rough number of API requests a book still needs, judging by what's cached. */
function expectedRequests(row: Row): number {
  const isbn = isbnsOf(row)[0];
  let workKey: string | undefined;
  if (isbn) {
    const edition = cached(booksApiUrl([isbn]));
    if (edition === undefined) return 1 / BATCH_SIZE + 1.15; // batch share + work (+ the odd search)
    workKey = toEdition(edition)?.workKey;
  }
  if (!workKey) {
    const search = cached(searchUrl(row));
    if (search === undefined) return 2; // search + work
    workKey = toDoc(search)?.key;
  }
  return workKey && /^\/works\/OL\d+W$/.test(workKey) && cached(workUrl(workKey)) === undefined ? 1 : 0;
}

const plural = (n: number, word: string) => `${n} ${word}${n === 1 ? "" : "s"}`;

function listBooks(lines: string[]) {
  for (const line of lines.slice(0, 25)) console.log(`  - ${line}`);
  if (lines.length > 25) console.log(`  …and ${lines.length - 25} more`);
}

// ── Main ────────────────────────────────────────────────────────────────────

async function main() {
  const csvPath = options.csv;
  if (!existsSync(csvPath)) {
    console.error(
      `Can't find ${csvPath}.\nExport it from Goodreads (My Books → Import and export → Export Library) and pass its path.`,
    );
    process.exit(1);
  }

  const rows = parseCsv(readFileSync(csvPath, "utf8")).filter((r) =>
    ["read", "currently-reading"].includes(r["Exclusive Shelf"]),
  );
  // Newest finished first; currently-reading leads.
  const sortKey = (r: Row) =>
    r["Exclusive Shelf"] === "currently-reading" ? "9999" : r["Date Read"] || `0000 ${r["Date Added"]}`;
  rows.sort((a, b) => sortKey(b).localeCompare(sortKey(a)));
  if (!rows.length) {
    // Most likely the wrong file: don't empty the shelf over it.
    fail(`${csvPath} has no books on the "read" or "currently-reading" shelves, so src/data/books.ts was not changed.`);
  }
  console.log(`Importing ${rows.length} books from ${csvPath}…`);

  // Keep whatever was fetched if the run is stopped or crashes.
  process.on("exit", saveCache);
  process.on("SIGINT", () => {
    saveCache();
    console.log("\nStopped. Progress is saved; run the same command again to continue where it left off.");
    process.exit(130);
  });
  // Node exits quietly with code 0 if a lookup somehow never settles; say so and fail instead.
  let finished = false;
  process.once("beforeExit", () => {
    if (finished) return;
    console.error("\nThe import stalled before finishing, so src/data/books.ts was not changed. Run it again.");
    process.exit(1);
  });

  const expected = rows.map(expectedRequests);
  const toLookUp = expected.filter((n) => n > 0).length;
  let estimate = Math.ceil(expected.reduce((a, b) => a + b, 0));
  const rate = `${plural(PER_SECOND, "request")}/second`;
  if (!toLookUp) {
    console.log(`All of them are cached in .cache/openlibrary.json, so this is quick (--refresh looks them up again).`);
  } else {
    const fromCache = rows.length - toLookUp;
    console.log(
      `${fromCache ? `${fromCache} are cached; looking up the other ${toLookUp}` : `Looking up ${plural(toLookUp, "book")}`} ` +
        `on Open Library: about ${plural(estimate, "request")} at ${rate}, so roughly ${roughly(estimate * api.interval)}.`,
    );
    if (!options.contact && estimate * api.interval >= 60_000) {
      const faster = roughly(estimate * spacing(3));
      console.log(
        `Tip: add --contact your@email (or set OPENLIBRARY_CONTACT) and it takes about ${faster} instead: ` +
          "Open Library allows 3x more requests from clients that leave a contact address.",
      );
    }
  }

  const startedAt = Date.now();
  const out: Book[] = new Array(rows.length);
  const notFound: string[] = [];
  const failed: string[] = [];
  let withoutCover = 0;
  let done = 0;
  const timeLeft = () => {
    estimate = Math.max(estimate, api.started + api.queued);
    const requests = Math.min(Math.max(estimate - api.started, api.queued), (rows.length - done) * 2.2);
    return requests >= 1 ? `  · ~${formatDuration(requests * api.interval)} left` : "";
  };

  await Promise.all(
    rows.map(async (row, i) => {
      const { book, matched, hasCover, failures } = await importRow(row, i);
      out[i] = book;
      done++;
      const label = `${book.title}, ${book.author}`;
      if (failures.length) failed.push(`${label} (${failures.join("; ")})`);
      else if (!matched) notFound.push(label);
      if (!hasCover) withoutCover++;
      const note = failures.length ? "  (incomplete, see the end)" : hasCover ? "" : "  (no cover found)";
      console.log(`[${done}/${rows.length}] ${book.title}${note}${timeLeft()}`);
    }),
  );

  const source = readFileSync(BOOKS_FILE, "utf8");
  const marker = source.indexOf("export const books");
  if (marker < 0) throw new Error(`Couldn't find "export const books" in ${BOOKS_FILE}`);
  const body = out.length ? `[\n${out.map(bookLiteral).join(",\n")},\n]` : "[]";
  writeFileSync(BOOKS_FILE, `${source.slice(0, marker)}export const books: Book[] = ${body};\n`);
  finished = true;
  saveCache();

  console.log(
    `\nWrote ${out.length} books to src/data/books.ts in ${formatDuration(Date.now() - startedAt)} ` +
      `(${api.started ? plural(api.started, "Open Library request") : "all from the cache"}).`,
  );
  if (withoutCover) console.log(`${withoutCover} without a cover on Open Library get a plain cloth-coloured spine.`);
  if (keptStale) console.log(`Couldn't refresh ${plural(keptStale, "lookup")}, so they keep what was cached before.`);
  if (notFound.length) {
    console.log(`\nNot found on Open Library (${notFound.length}) — kept with your Goodreads details:`);
    listBooks(notFound);
  }
  if (failed.length) {
    console.log(
      `\nFailed after retries (${failed.length}) — run the same command again to fill these in (everything else is cached):`,
    );
    listBooks(failed);
  }
}

main().catch((err) => {
  console.error(err);
  process.exit(1);
});
