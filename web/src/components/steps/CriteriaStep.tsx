import { useState } from "react";
import { GateRow } from "../criteria/GateRow";
import { GATES } from "../criteria/gates";
import { OFFTARGET_LONG, offtargetStatus } from "../../lib/engine";
import { RELAXED_HLA_PATCH } from "../../lib/params";
import { useSession } from "../../session/useSession";

/** compact: the expert-mode panel body, without heading, divider or navigation. */
export function CriteriaStep({ compact = false }: { compact?: boolean }) {
  const { setStep, runDesign, busy, params, updateParams, inspect, engine } = useSession();
  const offtarget = offtargetStatus(engine);
  const [presetApplied, setPresetApplied] = useState(false);
  const gc = inspect?.target.gc_percent;

  function applyRelaxed() {
    updateParams(RELAXED_HLA_PATCH);
    setPresetApplied(true);
  }

  return (
    <div className={compact ? "xp-body" : "workspace"}>
      {!compact && (
        <div className="whead">
          <h2>Design criteria</h2>
          <p>
            Every gate is independent and optional; a candidate is kept only if it clears each
            enabled one. These are the best-practice defaults - open any <b>i</b> marker for what a
            gate does and where to deviate.
          </p>
        </div>
      )}

      <div>
        {GATES.map((g) => (
          <GateRow
            key={g.id}
            gate={g}
            params={params}
            update={updateParams}
            notice={
              g.id === "offtarget" && (offtarget === "no-bowtie" || offtarget === "no-reference")
                ? OFFTARGET_LONG[offtarget]
                : undefined
            }
          />
        ))}
      </div>

      {compact ? (
        <div className="xp-preset">
          <button type="button" className="linkbtn" onClick={applyRelaxed}>
            {presetApplied ? "Relaxed HLA preset applied ✓" : "Apply relaxed HLA preset"}
          </button>
          {gc != null && <span className="hint">Target GC {gc.toFixed(1)}%</span>}
        </div>
      ) : (
        <>
          <div className="divider" />
          <div className="callout warn">
            <span className="ic">{"⚠"}</span>
            <span>
              <b>This target is GC-rich{gc != null ? ` (${gc.toFixed(1)}%)` : ""}.</b> The seed
              filter (Naito et al. 2009) is deliberately strict - on HLA it can reject every
              candidate. Recommended here: turn seed <b>off</b> and add the standalone GC and
              hairpin gates.
              <br />
              <button type="button" className="mini-btn" onClick={applyRelaxed}>
                {presetApplied ? "Relaxed HLA preset applied ✓" : "Apply relaxed HLA preset"}
              </button>
            </span>
          </div>
        </>
      )}

      {!compact && (
        <div className="actions">
          <button type="button" className="hbtn ghost" onClick={() => setStep("targeting")}>
            &larr; Back
          </button>
          <button type="button" className="hbtn" onClick={() => void runDesign()} disabled={busy}>
            {busy ? "Running..." : "Run design →"}
          </button>
        </div>
      )}
    </div>
  );
}
