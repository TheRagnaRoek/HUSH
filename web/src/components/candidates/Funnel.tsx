import { gateLabel } from "../../lib/text";

export function Funnel({
  funnel,
  offtargetSkipped,
}: {
  funnel: Record<string, number>;
  offtargetSkipped: boolean;
}) {
  const entries = Object.entries(funnel);
  const max = Math.max(1, ...entries.map(([, v]) => v));
  return (
    <>
      <div className="funnel-h">Selection funnel</div>
      <div className="funnel">
        {entries.map(([name, v]) => (
          <div className="frow" key={name}>
            <span className="fn">{gateLabel(name)}</span>
            <div className="ftrack">
              <div className="fbar" style={{ width: `${Math.max(5, (v / max) * 100)}%` }} />
            </div>
            <span className="fv num">{v}</span>
          </div>
        ))}
        {offtargetSkipped && (
          <div className="frow">
            <span className="fn">Off-target</span>
            <div className="ftrack">
              <div className="fbar skip" style={{ width: "100%" }} />
            </div>
            <span className="fv" style={{ color: "var(--muted)", fontWeight: 400 }}>
              Skip
            </span>
          </div>
        )}
      </div>
    </>
  );
}
