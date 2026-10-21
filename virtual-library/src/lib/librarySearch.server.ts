import Anthropic from "@anthropic-ai/sdk";
import { zodOutputFormat } from "@anthropic-ai/sdk/helpers/zod";
import { z } from "zod";
import { books } from "../data/books";

export type SmartSearchResult =
  | { status: "ok"; ids: string[] }
  | { status: "unconfigured" }
  | { status: "error"; message: string };

const SYSTEM_PROMPT = `You are the librarian of a small personal library. Match the reader's natural-language request against the catalogue below. Each line is: id :: title :: author :: year :: genres :: blurb.

Return ONLY a JSON object of the form {"ids":["id1","id2"]} — the ids of the books that fit the request, best match first, at most 20. Consider mood, themes, setting, genre and author, not just keywords. If nothing fits, return {"ids":[]}.`;

const IdsSchema = z.object({ ids: z.array(z.string()) });

// Built from the same `books` array the shelf renders, at request time.
function catalogue(): string {
  return books
    .map((b) =>
      [b.id, b.title, b.author, b.year || "", (b.genres ?? []).join(", "), b.blurb.slice(0, 160)].join(" :: "),
    )
    .join("\n");
}

function keepKnownIds(ids: unknown[]): string[] {
  const known = new Set(books.map((b) => b.id));
  const seen = new Set<string>();
  for (const id of ids) {
    if (typeof id === "string" && known.has(id)) seen.add(id);
  }
  return [...seen].slice(0, 20);
}

async function searchWithAnthropic(query: string, apiKey: string): Promise<SmartSearchResult> {
  const client = new Anthropic({ apiKey });
  try {
    const response = await client.messages.parse({
      model: process.env.ANTHROPIC_MODEL || "claude-haiku-5-5",
      max_tokens: 4096,
      // The catalogue is identical across searches, so cache it; the query follows.
      system: [
        {
          type: "text",
          text: `${SYSTEM_PROMPT}\n\nCatalogue:\n${catalogue()}`,
          cache_control: { type: "ephemeral" },
        },
      ],
      messages: [{ role: "user", content: query }],
      output_config: { effort: "low", format: zodOutputFormat(IdsSchema) },
    });
    if (response.stop_reason === "refusal") return { status: "ok", ids: [] };
    return { status: "ok", ids: keepKnownIds(response.parsed_output?.ids ?? []) };
  } catch (error) {
    if (error instanceof Anthropic.RateLimitError) {
      return { status: "error", message: "Too many searches — give it a moment and try again." };
    }
    if (error instanceof Anthropic.AuthenticationError || error instanceof Anthropic.PermissionDeniedError) {
      return { status: "error", message: "Search isn't configured correctly." };
    }
    if (error instanceof Anthropic.APIError) {
      console.error(`AI search failed (${error.status}):`, error.message);
      return { status: "error", message: "Search is unavailable right now." };
    }
    throw error;
  }
}

async function searchWithOpenAICompatible(query: string, apiKey: string): Promise<SmartSearchResult> {
  const base = (process.env.OPENAI_BASE_URL || "https://api.openai.com/v1").replace(/\/$/, "");
  const res = await fetch(`${base}/chat/completions`, {
    method: "POST",
    headers: { "Content-Type": "application/json", Authorization: `Bearer ${apiKey}` },
    body: JSON.stringify({
      model: process.env.OPENAI_MODEL || "gpt-4o-mini",
      response_format: { type: "json_object" },
      messages: [
        { role: "system", content: `${SYSTEM_PROMPT}\n\nCatalogue:\n${catalogue()}` },
        { role: "user", content: query },
      ],
    }),
  });

  if (res.status === 429) {
    return { status: "error", message: "Too many searches — give it a moment and try again." };
  }
  if (res.status === 402) return { status: "error", message: "Search credits exhausted." };
  if (!res.ok) {
    console.error(`AI search failed (${res.status}):`, await res.text().catch(() => ""));
    return { status: "error", message: "Search is unavailable right now." };
  }

  const json = (await res.json()) as { choices?: { message?: { content?: string } }[] };
  try {
    const parsed = IdsSchema.parse(JSON.parse(json.choices?.[0]?.message?.content ?? "{}"));
    return { status: "ok", ids: keepKnownIds(parsed.ids) };
  } catch {
    return { status: "ok", ids: [] };
  }
}

export async function smartSearchHandler(query: string): Promise<SmartSearchResult> {
  // Keys are read per request, inside the server handler — never in browser code.
  const anthropicKey = process.env.ANTHROPIC_API_KEY;
  if (anthropicKey) return searchWithAnthropic(query, anthropicKey);

  const openaiKey = process.env.OPENAI_API_KEY;
  if (openaiKey) return searchWithOpenAICompatible(query, openaiKey);

  return { status: "unconfigured" };
}
