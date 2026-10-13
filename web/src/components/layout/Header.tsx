import { useSession } from "../../session/useSession";
import { ThemeToggle } from "./ThemeToggle";

export function Header() {
  const { expert, setExpert } = useSession();
  return (
    <header className="topbar">
      <div className="brand">
        <span className="logo">HUSH</span>
        <span className="sep">/</span>
        <span className="tag">allele-aware siRNA &amp; shRNA designer &middot; Project HIDE</span>
      </div>
      <span className="spacer" />
      <span className="status">
        <span className="dot" />
        engine v0.4
      </span>
      <div className="toggle" role="group" aria-label="Mode">
        <button type="button" aria-pressed={!expert} onClick={() => setExpert(false)}>
          Guided
        </button>
        <button type="button" aria-pressed={expert} onClick={() => setExpert(true)}>
          Expert
        </button>
      </div>
      <ThemeToggle />
    </header>
  );
}
