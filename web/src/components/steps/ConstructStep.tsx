import { ConstructViewer } from "../viewer/ConstructViewer";
import { useSession } from "../../session/useSession";

export function ConstructStep() {
  const {
    design,
    inspect,
    selectedCandidate,
    selectedConstruct,
    selectedRegion,
    enrichment,
    runConvert,
    busy,
  } = useSession();

  if (!design) return null;
  const rank = design.selected.findIndex((c) => c.region === selectedRegion) + 1;
  const targetId = inspect?.target.id ?? design.target.id;

  return (
    <div className="workspace">
      <div className="whead">
        <h2>Construct viewer</h2>
        <p>
          The selected candidate, grafted onto a stem-loop scaffold - the same viewer as the
          candidates step, opened straight to the shRNA build.
        </p>
      </div>

      {selectedCandidate && !selectedConstruct && (
        <div className="actions" style={{ marginTop: 0, marginBottom: 16 }}>
          <button type="button" className="hbtn" onClick={() => void runConvert()} disabled={busy}>
            {busy ? "Converting..." : "Convert selected to shRNA →"}
          </button>
        </div>
      )}

      {selectedCandidate && (
        <ConstructViewer
          candidate={selectedCandidate}
          construct={selectedConstruct}
          enrichment={enrichment}
          criteria={design.criteria}
          offtargetSkipped={design.offtarget_skipped}
          rank={rank}
          targetId={targetId}
          startTab="build"
        />
      )}
    </div>
  );
}
