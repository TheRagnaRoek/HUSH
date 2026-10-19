import { CandidateList } from "../candidates/CandidateList";
import { Funnel } from "../candidates/Funnel";
import { ConstructViewer } from "../viewer/ConstructViewer";
import { WithheldSummary } from "../viewer/WithheldSummary";
import { terminatorString } from "../../lib/loops";
import { RELAXED_HLA_PATCH } from "../../lib/params";
import { OFFTARGET_REASON, offtargetStatus } from "../../lib/engine";
import { POSITION_ORDER, rankMode } from "../../lib/targeting";
import { acronyms } from "../../lib/text";
import { exportConstruct, exportDesign, type ExportFormat } from "../../lib/tsv";
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

/** compact: the expert-mode results pane (funnel, list and viewer side by side). */
export function CandidatesStep({ compact = false }: { compact?: boolean }) {
  const {
    design,
    inspect,
    convert,
    selectedCandidate,
    selectedConstruct,
    selectedWithheld,
    selectedRegion,
    enrichment,
    offline,
    error,
    setStep,
    params,
    runDesign,
    runConvert,
    engine,
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
  const mode = rankMode(design.targeting);
  const offtarget = offtargetStatus(engine);

  const summary = (
    <>
      {design.selected.length} of {design.candidate_count} enumerated windows cleared every gate.
    </>
  );

  const notices = (
    <>
      {error && (
        <div className="banner err">
          Your design did not run: {error} The results below are representative
          sample data, not your target.
        </div>
      )}
      {offline && !error && (
        <div className="banner">
          {compact
            ? "Representative sample data. Load a target and run a design to see your own results."
            : "Representative sample data. Load a target in step 1 and run a design to see your own results."}
        </div>
      )}
      {design.offtarget_skipped && !offline && !error && (
        <div className="banner warn">
          Off-target screening was skipped, so these guides have not been checked against the
          human transcriptome.{" "}
          {offtarget === "no-bowtie" || offtarget === "no-reference"
            ? OFFTARGET_REASON[offtarget]
            : "No reference index was available to the engine."}
        </div>
      )}
      {mode !== "discriminating" && !error && (
        <div className="banner">{POSITION_ORDER[mode]}</div>
      )}

      {emptyRun && (
        <div className="callout warn">
          <span className="ic">{"⚠"}</span>
          <span>
            <b>No candidate cleared every gate.</b>
            {narrowest && (
              <>
                {" "}The narrowest was <b>{acronyms(narrowest.name)}</b>, which{" "}
                {narrowest.passed} of {design.candidate_count} passed.
              </>
            )}
            {seedIsOn ? (
              <>
                {" "}The seed filter (Naito et al. 2009) is deliberately
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
                {" "}Loosen the criteria
                {compact ? (
                  " "
                ) : (
                  <>
                    {" "}in{" "}
                    <button
                      type="button"
                      className="linkbtn"
                      onClick={() => setStep("criteria")}
                    >
                      step 3
                    </button>{" "}
                  </>
                )}
                and run again.
              </>
            )}
          </span>
        </div>
      )}
    </>
  );

  const listHead = (
    <div className="listhead">
      <span className="lh">
        {mode === "discriminating" ? "Selected - by discrimination" : "Selected - by target position"}
      </span>
      <span className="exports">
        <button
          type="button"
          className="linkbtn"
          disabled={design.selected.length === 0}
          onClick={() => exportDesign(design, offline, "tsv")}
        >
          Export TSV
        </button>
        <button
          type="button"
          className="linkbtn"
          disabled={design.selected.length === 0}
          onClick={() => exportDesign(design, offline, "txt")}
        >
          Export TXT
        </button>
      </span>
    </div>
  );

  const onExport = selectedConstruct
    ? (format: ExportFormat) => exportConstruct(selectedConstruct, targetId, offline, format)
    : undefined;

  if (compact) {
    return (
      <div className="xp-body">
        <p className="xp-sum">{summary}</p>
        {notices}
        <WithheldSummary />
        <div className="xp-top">
          <div>
            <Funnel funnel={design.funnel} offtargetSkipped={design.offtarget_skipped} />
          </div>
          <div>
            {listHead}
            <CandidateList />
          </div>
        </div>
        {selectedCandidate && (
          <ConstructViewer
            candidate={selectedCandidate}
            construct={selectedConstruct}
            withheld={selectedWithheld}
            enrichment={enrichment}
            criteria={design.criteria}
            offtargetSkipped={design.offtarget_skipped}
            rank={rank}
            mode={mode}
            targetId={targetId}
            terminator={terminator}
            startTab="align"
            keepTab
            onExport={onExport}
          />
        )}
      </div>
    );
  }

  return (
    <div className="workspace">
      <div className="whead">
        <h2>Candidates</h2>
        <p>
          {summary} The funnel shows where the rest fell out; pick a candidate to load it into the
          viewer.
        </p>
      </div>

      {notices}

      <div className="split">
        <div className="leftcol">
          <Funnel funnel={design.funnel} offtargetSkipped={design.offtarget_skipped} />
          {listHead}
          <CandidateList />
        </div>
        <div className="rightcol">
          {selectedCandidate && (
            <ConstructViewer
              candidate={selectedCandidate}
              construct={selectedConstruct}
              withheld={selectedWithheld}
              enrichment={enrichment}
              criteria={design.criteria}
              offtargetSkipped={design.offtarget_skipped}
              rank={rank}
              mode={mode}
              targetId={targetId}
              terminator={terminator}
              startTab="align"
              buildView="overview"
              onOpenConstruct={
                selectedConstruct || selectedWithheld
                  ? () => setStep("construct")
                  : () => void runConvert()
              }
              onExport={onExport}
            />
          )}
        </div>
      </div>

      <div className="actions">
        <button type="button" className="hbtn ghost" onClick={() => setStep("criteria")}>
          &larr; Back
        </button>
        <button
          type="button"
          className="hbtn"
          disabled={busy || design.selected.length === 0}
          onClick={() => (convert ? setStep("construct") : void runConvert())}
        >
          {busy ? "Converting..." : convert ? "Continue to shRNA construct →" : "Convert to shRNA →"}
        </button>
        <span className="hint">
          Optional. The siRNA duplex of each candidate is ready to use as it is; convert to get a
          hairpin to express from a vector.
        </span>
      </div>
    </div>
  );
}
