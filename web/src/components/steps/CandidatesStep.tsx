import { CandidateList } from "../candidates/CandidateList";
import { Funnel } from "../candidates/Funnel";
import { ConstructViewer } from "../viewer/ConstructViewer";
import { terminatorString } from "../../lib/loops";
import { useSession } from "../../session/useSession";

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
  } = useSession();

  if (!design) return null;
  const rank = design.selected.findIndex((c) => c.region === selectedRegion) + 1;
  const targetId = inspect?.target.id ?? design.target.id;
  const terminator = terminatorString(convert?.terminator_length);

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

      {error && <div className="banner err">API error: {error}. Showing sample data.</div>}
      {offline && !error && (
        <div className="banner">Showing representative sample data - no backend connected.</div>
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
