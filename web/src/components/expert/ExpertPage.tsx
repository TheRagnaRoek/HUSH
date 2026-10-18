import type { ReactNode } from "react";
import { CandidatesStep } from "../steps/CandidatesStep";
import { CriteriaStep } from "../steps/CriteriaStep";
import { InputStep } from "../steps/InputStep";
import { TargetingStep } from "../steps/TargetingStep";
import { ScaffoldOptions } from "../viewer/ScaffoldOptions";
import { useSession } from "../../session/useSession";

function Section({ title, children }: { title: string; children: ReactNode }) {
  return (
    <details className="xp-sec" open>
      <summary>{title}</summary>
      {children}
    </details>
  );
}

/** Expert mode: every setting in one column, the results beside it. */
export function ExpertPage() {
  const { runDesign, runConvert, busy, design, inspect, offline } = useSession();
  // The targeting choice is seeded from the panel, so it restarts with a new one.
  const panelKey = (inspect?.panel?.members ?? []).map((m) => m.id).join("|");

  return (
    <div className="xp">
      <aside className="xp-side">
        <div className="xp-bar">
          <button
            type="button"
            className="hbtn"
            onClick={() => void runDesign()}
            disabled={busy}
          >
            {busy ? "Running..." : "Run design"}
          </button>
          <button
            type="button"
            className="hbtn ghost"
            onClick={() => void runConvert()}
            disabled={busy || !design}
          >
            Convert to shRNA
          </button>
          {offline && <span className="hint">Sample data</span>}
        </div>
        <Section title="Target & panel">
          <InputStep compact />
        </Section>
        <Section title="Targeting">
          <TargetingStep key={panelKey} compact />
        </Section>
        <Section title="Criteria">
          <CriteriaStep compact />
        </Section>
        <Section title="shRNA scaffold">
          <div className="xp-body">
            <ScaffoldOptions />
          </div>
        </Section>
      </aside>
      <main className="xp-main">
        <CandidatesStep compact />
      </main>
    </div>
  );
}
