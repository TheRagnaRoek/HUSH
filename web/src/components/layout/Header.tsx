import { useEffect, useState } from "react";

import { health } from "../../api/hush";
import { useSession } from "../../session/useSession";
import { ThemeToggle } from "./ThemeToggle";

export function Header() {
  const { expert, setExpert } = useSession();
  // The indicator reports what the API actually says. With no backend
  // reachable the app still works on sample data, and says so.
  const [version, setVersion] = useState<string | null>(null);
  useEffect(() => {
    let cancelled = false;
    health()
      .then((h) => {
        if (!cancelled) setVersion(h.version);
      })
      .catch(() => {
        if (!cancelled) setVersion(null);
      });
    return () => {
      cancelled = true;
    };
  }, []);
  return (
    <header className="topbar">
      <div className="brand">
        <img className="logo" src="/hush_logo.svg" alt="HUSH" />
        <span className="sep">/</span>
        <span className="tag">allele-aware siRNA &amp; shRNA designer &middot; Project HIDE</span>
      </div>
      <span className="spacer" />
      <span
        className="status"
        title={
          version
            ? "Connected to the HUSH engine"
            : "No engine reachable - showing sample data"
        }
      >
        <span className={version ? "dot" : "dot off"} />
        {version ? `engine v${version}` : "sample data"}
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
