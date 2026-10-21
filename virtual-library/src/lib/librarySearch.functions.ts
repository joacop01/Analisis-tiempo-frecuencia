import { createServerFn } from "@tanstack/react-start";
import { smartSearchHandler } from "./librarySearch.server";

export type { SmartSearchResult } from "./librarySearch.server";

export const smartSearch = createServerFn({ method: "POST" })
  .validator((data: unknown) => {
    const query =
      typeof data === "object" && data !== null && "query" in data && typeof data.query === "string"
        ? data.query.trim()
        : "";
    if (query.length < 2 || query.length > 300) throw new Error("Query must be 2–300 characters.");
    return { query };
  })
  .handler(({ data }) => smartSearchHandler(data.query));
