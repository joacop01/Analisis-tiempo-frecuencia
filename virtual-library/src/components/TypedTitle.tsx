import { useEffect, useState } from "react";

const FULL = "Welcome to my library";
const STEP_MS = 95;

export function TypedTitle() {
  const [count, setCount] = useState(0);
  const done = count >= FULL.length;

  useEffect(() => {
    if (done) return;
    const t = setTimeout(() => setCount((c) => c + 1), STEP_MS);
    return () => clearTimeout(t);
  }, [count, done]);

  return (
    <h1
      aria-label={FULL}
      className="font-display text-5xl font-light italic leading-[1.05] tracking-tight text-foreground sm:text-6xl md:text-7xl"
    >
      <span aria-hidden="true">
        {FULL.slice(0, count)}
        <span
          className={`ml-1 inline-block h-[0.82em] w-[2px] translate-y-[0.08em] bg-primary/80 align-baseline ${
            done ? "animate-caret" : ""
          }`}
        />
      </span>
    </h1>
  );
}
