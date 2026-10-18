import { NavLink } from "react-router-dom";
import { OFFTARGET_SHORT, offtargetStatus } from "../../lib/engine";
import { useSession } from "../../session/useSession";
import { ThemeToggle } from "./ThemeToggle";

export function Header() {
  const { expert, setExpert, engine } = useSession();
  const version = engine?.version ?? null;
  // A reachable engine can still be unable to run the off-target check.
  const offtarget = offtargetStatus(engine);
  const problem = offtarget === "no-bowtie" || offtarget === "no-reference" ? offtarget : null;
  return (
    <header className="topbar">
      <div className="brand">
        <img className="logo" src="/hush_logo.svg" alt="HUSH" />
        <span className="sep">/</span>
        <span className="tag">
          <span>Allele-aware siRNA &amp; shRNA designer</span>
          <span>
            <a href="https://2026.igem.wiki/gu-frankfurt/" target="_blank" rel="noopener noreferrer">
              Designed by Team HIDE
            </a>
          </span>
        </span>
      </div>
      <span className="spacer" />
      <span
        className="status"
        title={
          version
            ? `Connected to the HUSH engine${problem ? "\n" + OFFTARGET_SHORT[problem] : ""}`
            : "No engine reachable - showing sample data"
        }
      >
        <span className={!version ? "dot off" : problem ? "dot warn" : "dot"} />
        {version ? `Engine v${version}` : "Sample data"}
      </span>
      <NavLink to="/reference" className={({ isActive }) => (isActive ? "navlink active" : "navlink")}>
        References
      </NavLink>
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
