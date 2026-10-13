import { useEffect, useMemo, useState } from "react";
import type { TargetingSchema } from "../../api/types";
import { useSession } from "../../session/useSession";

type Intent = "pan" | "allele" | "target" | "custom";

const INTENTS: { id: Intent; t: string; d: string }[] = [
  {
    id: "allele",
    t: "Allele-specific",
    d: "Silence one named panel member, discriminate against the rest. Use when the allele to silence is in the panel.",
  },
  { id: "pan", t: "Pan-allele", d: "Conserved sites across the whole family. Hit every member, spare none." },
  {
    id: "target",
    t: "Target-specific",
    d: "Silence the loaded target, spare the entire reference panel. Use when the target sits apart from its panel.",
  },
  { id: "custom", t: "Custom hit / spare", d: "Assign each panel member to hit or spare explicitly." },
];

function initialIntent(t: TargetingSchema | null | undefined, hasPanel: boolean): Intent {
  if (t?.allele_specific) return "allele";
  if (t?.target_specific) return "target";
  if (t?.hit) return "custom";
  if (t?.pan_allele) return "pan";
  return hasPanel ? "allele" : "pan";
}

export function TargetingStep() {
  const { inspect, params, updateParams, setStep } = useSession();
  const members = useMemo(() => inspect?.panel?.members ?? [], [inspect]);
  const memberIds = useMemo(() => members.map((m) => m.id), [members]);
  const hasPanel = memberIds.length > 0;

  const [intent, setIntent] = useState<Intent>(() =>
    initialIntent(params.targeting, memberIds.length > 0),
  );
  const [specificId, setSpecificId] = useState<string>(
    () => params.targeting?.allele_specific ?? memberIds[0] ?? "",
  );
  const [hitSet, setHitSet] = useState<Set<string>>(
    () => new Set(params.targeting?.hit ?? (memberIds[0] ? [memberIds[0]] : [])),
  );

  // The UI selection is the single source of truth for params.targeting.
  const schema = useMemo<TargetingSchema>(() => {
    if (intent === "pan") return { pan_allele: true };
    if (intent === "allele") return { allele_specific: specificId };
    if (intent === "target") return { target_specific: true };
    return { hit: [...hitSet], spare: memberIds.filter((id) => !hitSet.has(id)) };
  }, [intent, specificId, hitSet, memberIds]);

  useEffect(() => {
    updateParams({ targeting: hasPanel ? schema : null });
  }, [schema, hasPanel, updateParams]);

  function roleOf(id: string): "hit" | "spare" {
    if (intent === "pan") return "hit";
    if (intent === "target") return "spare";
    if (intent === "allele") return id === specificId ? "hit" : "spare";
    return hitSet.has(id) ? "hit" : "spare";
  }

  function onMember(id: string) {
    if (intent === "allele") setSpecificId(id);
    else if (intent === "custom") {
      setHitSet((prev) => {
        const next = new Set(prev);
        if (next.has(id)) next.delete(id);
        else next.add(id);
        return next;
      });
    }
  }

  const interactive = intent === "allele" || intent === "custom";
  const discriminates = intent === "allele" || intent === "target" || intent === "custom";

  return (
    <div className="workspace">
      <div className="whead">
        <h2>What should this construct hit and spare?</h2>
        <p>
          The allele gate is ANDed on top of every efficacy filter: a candidate must hit the
          alleles you want silenced and carry enough mismatches against the ones you want spared.
        </p>
      </div>

      {!hasPanel && (
        <div className="banner">
          No panel loaded - add one in step 1 to use allele targeting. A design with no panel runs
          the efficacy filters alone.
        </div>
      )}

      <div className="intent">
        {INTENTS.map((it) => (
          <label key={it.id}>
            <input
              type="radio"
              name="intent"
              checked={intent === it.id}
              onChange={() => setIntent(it.id)}
            />
            <span className="t">{it.t}</span>
            <span className="d">{it.d}</span>
          </label>
        ))}
      </div>

      {hasPanel && (
        <>
          <div className="divider" />
          <div className="funnel-h">
            Panel members
            {interactive && (
              <>
                {" - "}
                {intent === "allele" ? "click the one to silence" : "click to toggle hit / spare"}
              </>
            )}
          </div>
          <div className="members">
            {members.map((m) => {
              const role = roleOf(m.id);
              const cls = `allele${role === "hit" ? " hit" : ""}`;
              return interactive ? (
                <button type="button" className={cls} key={m.id} onClick={() => onMember(m.id)}>
                  {m.id} <span className="role">{role}</span>
                </button>
              ) : (
                <span className={cls} key={m.id}>
                  {m.id} <span className="role">{role}</span>
                </span>
              );
            })}
          </div>
        </>
      )}

      <p className="hint" style={{ marginTop: 14 }}>
        {discriminates
          ? "A request that spares something is ranked by discrimination by default - survivors are ordered by how well they distinguish the weakest spared allele."
          : "Pan-allele targeting spares nothing, so there is no discrimination to rank; survivors come back in target-position order."}
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
