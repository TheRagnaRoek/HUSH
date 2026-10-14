import { NoSelfExcludedCallout } from "../criteria/SelfTranscripts";
import { ConstructViewer } from "../viewer/ConstructViewer";
import { ScaffoldOptions } from "../viewer/ScaffoldOptions";
import { terminatorString } from "../../lib/loops";
import { useSession } from "../../session/useSession";

export function ConstructStep() {
  const {
    design,
    inspect,
    convert,
    selectedCandidate,
    selectedConstruct,
    selectedRegion,
    enrichment,
    offline,
    runConvert,
    busy,
  } = useSession();

  if (!design) return null;
  const rank = design.selected.findIndex((c) => c.region === selectedRegion) + 1;
  const targetId = inspect?.target.id ?? design.target.id;
  const terminator = terminatorString(convert?.terminator_length);
  const converted = Boolean(selectedConstruct);

  return (
    <div className="workspace">
      <div className="whead">
        <h2>Construct viewer</h2>
        <p>
          The selected candidate, grafted onto a stem-loop scaffold - the same viewer as the
          candidates step, opened straight to the shRNA build.
        </p>
      </div>

      <NoSelfExcludedCallout criteria={design.criteria} />

      <div className="funnel-h">Scaffold options</div>
      <ScaffoldOptions />
      <div className="actions" style={{ marginTop: 14, marginBottom: 16 }}>
        <button type="button" className="hbtn" onClick={() => void runConvert()} disabled={busy}>
          {busy
            ? "Converting..."
            : converted
              ? "Re-convert with these options →"
              : "Convert selected to shRNA →"}
        </button>
        {converted && (
          <span className="hint">Changed an option? Re-convert to apply it.</span>
        )}
        {offline && (
          <span className="hint">Sample mode - options apply against a connected backend.</span>
        )}
      </div>
      <div className="divider" />

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
          startTab="build"
        />
      )}
    </div>
  );
}
