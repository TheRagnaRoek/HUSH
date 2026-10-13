import { useCallback, useEffect, useState } from "react";

export type ThemeMode = "light" | "dark";
const KEY = "hush-theme";

function read(): ThemeMode {
  try {
    const v = localStorage.getItem(KEY);
    if (v === "light" || v === "dark") return v;
  } catch {
    // storage blocked (private window etc.) - fall through to default
  }
  return "light"; // light is the default
}

/** Theme is an explicit choice persisted per viewer; it sets data-theme on the
 * document root (tokens.css reads it). Light is the default. */
export function useTheme(): [ThemeMode, (m: ThemeMode) => void] {
  const [mode, setMode] = useState<ThemeMode>(read);

  useEffect(() => {
    document.documentElement.setAttribute("data-theme", mode);
    try {
      localStorage.setItem(KEY, mode);
    } catch {
      // ignore persistence failures
    }
  }, [mode]);

  const set = useCallback((m: ThemeMode) => setMode(m), []);
  return [mode, set];
}
