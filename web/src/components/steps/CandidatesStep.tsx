import { CandidateList } from "../candidates/CandidateList";
import { Funnel } from "../candidates/Funnel";
import { ConstructViewer } from "../viewer/ConstructViewer";
import { terminatorString } from "../../lib/loops";
import { RELAXED_HLA_PATCH } from "../../lib/params";
import { useSession } from "../../session/useSession";

/** The gate fewest candidates cleared: the one to loosen first.
 *
 * Funnel counts are per-gate and independent, not cumulative, so the
 * smallest is the narrowest gate. 'selected' is the outcome, not a gate.
 */
function narrowestGate(
  funnel: Record<string, number>,
): { name: string; passed: number } | null {
  const gates = Object.entries(funnel).filter(([name]) => name !== "selected");
  if (gates.length === 0) return null;
  const [name, passed] = gates.reduce((a, b) => (b[1] < a[1] ? b : a));
  return { name, passed };
}

export function CandidatesStep() {
  const {
    design,
    inspect,
    convert,
    selectedCandidate,
    selectedConstruct,
    selectedRegion,
    enrichment,
    offline,
    error,
    setStep,
    params,
    runDesign,
    busy,
  } = useSession();

  if (!design) return null;
  const rank = design.selected.findIndex((c) => c.region === selectedRegion) + 1;
  const targetId = inspect?.target.id ?? design.target.id;
  const terminator = terminatorString(convert?.terminator_length);
  // Only for a real, finished run: not for the sample, after an error, or
  // while busy (applying the preset changes params before the re-run returns).
  const emptyRun = !offline && !error && !busy && design.selected.length === 0;
  const narrowest = emptyRun ? narrowestGate(design.funnel) : null;
  const seedIsOn = params.use_seed_filter;

  return (
    <div className="workspace">
      <div className="whead">
        <h2>Candidates</h2>
        <p>
          {design.selected.length} of {design.candidate_count} enumerated windows cleared every
          gate. The funnel shows where the rest fell out; pick a candidate to load it into the
          viewer.
        </p>
      </div>

      {error && (
        <div className="banner err">
          Your design did not run: {error} The results below are representative
          sample data, not your target.
        </div>
      )}
      {offline && !error && (
        <div className="banner">
          Representative sample data. Load a target in step 1 and run a design to
          see your own results.
        </div>
      )}

      {emptyRun && (
        <div className="callout warn">
          <span className="ic">{"\u26A0"}</span>
          <span>
            <b>No candidate cleared every gate.</b>
            {narrowest && (
              <>
                {" "}The narrowest was <b>{narrowest.name}</b>, which{" "}
                {narrowest.passed} of {design.candidate_count} passed.
              </>
            )}
            {seedIsOn ? (
              <>
                {" "}The seed filter is the siDirect standard but deliberately
                strict - on a GC-rich target it can reject everything.
                <br />
                <button
                  type="button"
                  className="mini-btn"
                  disabled={busy}
                  onClick={() => void runDesign(RELAXED_HLA_PATCH)}
                >
                  {busy ? "Running..." : "Apply the relaxed preset and re-run"}
                </button>
              </>
            ) : (
              <>
                {" "}Loosen it in{" "}
                <button
                  type="button"
                  className="linkbtn"
                  onClick={() => setStep("criteria")}
                >
                  step 3
                </button>{" "}
                and run again.
              </>
            )}
          </span>
        </div>
      )}

      <div className="split">
        <div className="leftcol">
          <Funnel funnel={design.funnel} offtargetSkipped={design.offtarget_skipped} />
          <div className="listhead">
            <span className="lh">Selected - by discrimination</span>
            <button type="button" className="linkbtn">
              Export TSV
            </button>
          </div>
          <CandidateList />
        </div>
        <div className="rightcol">
          {selectedCandidate && (
            <ConstructViewer
              candidate={selectedCandidate}
              construct={selectedConstruct}
              enrichment={enrichment}
              criteria={design.criteria}
              offtargetSkipped={design.offtarget_skipped}
              rank={rank}
              targetId={targetId}
              terminator={terminator}
              startTab="align"
              onOpenConstruct={() => setStep("construct")}
            />
          )}
        </div>
      </div>
    </div>
  );
}
