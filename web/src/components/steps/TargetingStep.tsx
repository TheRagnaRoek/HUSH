import { useSession } from "../../session/useSession";

const INTENTS = [
  {
    t: "Allele-specific",
    d: "Silence one named panel member, discriminate against the rest. Use when the allele to silence is in the panel.",
    checked: true,
  },
  { t: "Pan-allele", d: "Conserved sites across the whole family. Hit every member, spare none.", checked: false },
  {
    t: "Target-specific",
    d: "Silence the loaded target, spare the entire reference panel. Use when the target sits apart from its panel.",
    checked: false,
  },
  { t: "Custom hit / spare", d: "Name explicit subsets to hit and to discriminate against.", checked: false },
];

export function TargetingStep() {
  const { inspect, design, setStep } = useSession();
  const members = inspect?.panel?.members ?? [];
  const hit = new Set(design?.targeting?.hit ?? []);

  return (
    <div className="workspace">
      <div className="whead">
        <h2>What should this construct hit and spare?</h2>
        <p>
          The allele gate is ANDed on top of every efficacy filter: a candidate must hit the
          alleles you want silenced and carry enough mismatches against the ones you want spared.
        </p>
      </div>

      <div className="intent">
        {INTENTS.map((it) => (
          <label key={it.t}>
            <input type="radio" name="intent" defaultChecked={it.checked} />
            <span className="t">{it.t}</span>
            <span className="d">{it.d}</span>
          </label>
        ))}
      </div>

      <div className="divider" />
      <div className="funnel-h">Panel members - click to assign roles</div>
      <div className="members">
        {members.map((m) => {
          const isHit = hit.has(m.id);
          return (
            <span className={`allele${isHit ? " hit" : ""}`} key={m.id}>
              {m.id} <span className="role">{isHit ? "hit" : "spare"}</span>
            </span>
          );
        })}
      </div>

      <p className="hint" style={{ marginTop: 14 }}>
        A request that spares something is ranked by discrimination by default - survivors are
        ordered by how well they distinguish the weakest spared allele.
      </p>

      <div className="actions">
        <button type="button" className="hbtn ghost" onClick={() => setStep("input")}>
          &larr; Back
        </button>
        <button type="button" className="hbtn" onClick={() => setStep("criteria")}>
          Continue &rarr;
        </button>
      </div>
    </div>
  );
}
