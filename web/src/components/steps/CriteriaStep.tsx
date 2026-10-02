import { useState } from "react";
import { GateRow } from "../criteria/GateRow";
import { GATES } from "../criteria/gates";
import { useSession } from "../../session/useSession";

function initialOn(): Record<string, boolean> {
  return Object.fromEntries(GATES.map((g) => [g.id, g.defaultOn]));
}

export function CriteriaStep() {
  const { setStep, runDesign, busy } = useSession();
  const [on, setOn] = useState<Record<string, boolean>>(initialOn);
  const [presetApplied, setPresetApplied] = useState(false);

  function applyRelaxed() {
    setOn((prev) => ({ ...prev, seed: false, gc: true, hairpin: true }));
    setPresetApplied(true);
  }

  return (
    <div className="workspace">
      <div className="whead">
        <h2>Design criteria</h2>
        <p>
          Every gate is independent and optional; a candidate is kept only if it clears each
          enabled one. These are the best-practice defaults - open any <b>?</b> for what a gate
          does and where to deviate.
        </p>
      </div>

      <div>
        {GATES.map((g) => (
          <GateRow
            key={g.id}
            gate={g}
            on={on[g.id]}
            onToggle={(v) => setOn((prev) => ({ ...prev, [g.id]: v }))}
          />
        ))}
      </div>

      <div className="divider" />
      <div className="callout warn">
        <span className="ic">{"⚠"}</span>
        <span>
          <b>This target is GC-rich (64.8%).</b> The seed filter is the siDirect industry standard
          but deliberately strict - on HLA it can reject every candidate. Recommended here: turn
          seed <b>off</b> and add the standalone GC and hairpin gates.
          <br />
          <button type="button" className="mini-btn" onClick={applyRelaxed}>
            {presetApplied ? "Relaxed HLA preset applied ✓" : "Apply relaxed HLA preset"}
          </button>
        </span>
      </div>

      <div className="actions">
        <button type="button" className="hbtn ghost" onClick={() => setStep("targeting")}>
          &larr; Back
        </button>
        <button type="button" className="hbtn" onClick={() => void runDesign()} disabled={busy}>
          {busy ? "Running..." : "Run design →"}
        </button>
      </div>
    </div>
  );
}
