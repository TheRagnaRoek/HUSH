import { useTheme } from "../../hooks/useTheme";

export function ThemeToggle() {
  const [mode, setMode] = useTheme();
  return (
    <div className="toggle" role="group" aria-label="Theme">
      <button type="button" aria-pressed={mode === "light"} onClick={() => setMode("light")}>
        {"☀"} Light
      </button>
      <button type="button" aria-pressed={mode === "dark"} onClick={() => setMode("dark")}>
        {"☾"} Dark
      </button>
    </div>
  );
}
