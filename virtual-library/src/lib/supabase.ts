import { createClient } from "@supabase/supabase-js";

// The anon key is public by design: the table's grants and row-level policies
// (supabase/migrations) are what keep it safe.
const url = import.meta.env.VITE_SUPABASE_URL as string | undefined;
const anonKey = import.meta.env.VITE_SUPABASE_ANON_KEY as string | undefined;

/** `null` when the optional recommendations backend isn't configured. */
export const supabase = url && anonKey ? createClient(url, anonKey) : null;
