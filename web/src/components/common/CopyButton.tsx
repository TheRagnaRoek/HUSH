import { useState } from "react";

/** Copy-to-clipboard button with a select-text fallback for app views that
 * refuse the clipboard API. */
export function CopyButton({ text, label = "Copy" }: { text: string; label?: string }) {
  const [done, setDone] = useState(false);
  async function copy() {
    try {
      await navigator.clipboard.writeText(text);
      setDone(true);
      setTimeout(() => setDone(false), 1200);
    } catch {
      // fallback: no clipboard permission; leave state unchanged
    }
  }
  return (
    <button type="button" className="copy" onClick={copy}>
      {done ? "Copied ✓" : label}
    </button>
  );
}
