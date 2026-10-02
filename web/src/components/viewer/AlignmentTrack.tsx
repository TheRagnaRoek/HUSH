import { useState } from "react";
import type { AlignmentModel } from "../../lib/alignment";
import { discProse } from "../../lib/disc";

interface TipState {
  show: boolean;
  x: number;
  y: number;
  title: string;
  lines: string[];
}

const EMPTY: TipState = { show: false, x: 0, y: 0, title: "", lines: [] };

export function AlignmentTrack({
  model,
  targetId,
}: {
  model: AlignmentModel;
  targetId: string;
}) {
  const [tip, setTip] = useState<TipState>(EMPTY);
  const target = [...model.target];
  const guideRev = [...model.guide].reverse();
  const L = target.length;
  const cols = Array.from({ length: L }, (_, i) => i);

  const weakest = model.alleles.find((a) => a.weak) ?? model.alleles[0];

  return (
    <div className="vpane">
      <div className="legend">
        <span>
          <span
            className="sw"
            style={{ background: "var(--brand-soft)", boxShadow: "inset 0 0 0 1px var(--brand)" }}
          />
          target / guide
        </span>
        <span>
          <span className="sw" style={{ background: "var(--mismatch)" }} />
          mismatch vs a spared allele
        </span>
        <span>hover a highlighted base for its discrimination criteria</span>
      </div>

      <div className="aln">
        <div className="alnrow">
          <span className="rl">position</span>
          <div className="track">
            {cols.map((i) => (
              <span className="poslab" key={i}>
                {i + 1}
              </span>
            ))}
          </div>
        </div>

        <div className="alnrow target">
          <span className="rl">
            <b>{targetId}</b>target &middot; sense
          </span>
          <div className="track">
            {target.map((b, i) => (
              <span className="base t" key={i}>
                {b}
              </span>
            ))}
          </div>
        </div>

        <div className="alnrow">
          <span className="rl">guide 3&prime;&rarr;5&prime;</span>
          <div className="track">
            {guideRev.map((b, i) => (
              <span className="base dim" key={i}>
                {b}
              </span>
            ))}
          </div>
        </div>

        <div className="alnrow sep">
          <span className="rl" style={{ color: "var(--cyan-ink)" }}>
            spared alleles
          </span>
          <div className="track" />
        </div>

        {model.alleles.map((a) => (
          <div className="alnrow" key={a.name}>
            <span className="rl">
              <b>{a.name}</b>
              {a.diff} diff{a.weak ? " · weakest" : ""}
            </span>
            <div className="track">
              {cols.map((i) => {
                const m = a.mism.find((x) => x.col === i);
                if (!m) {
                  return (
                    <span className="base dim" key={i}>
                      &middot;
                    </span>
                  );
                }
                return (
                  <span
                    className="base mm"
                    key={i}
                    onMouseEnter={(e) =>
                      setTip({
                        show: true,
                        x: e.clientX + 12,
                        y: e.clientY + 12,
                        title: `${m.from}→${m.to} @ pos P${m.disc.guide_position}`,
                        lines: discProse(m.disc),
                      })
                    }
                    onMouseMove={(e) =>
                      setTip((t) => ({ ...t, x: e.clientX + 12, y: e.clientY + 12 }))
                    }
                    onMouseLeave={() => setTip(EMPTY)}
                  >
                    {m.to}
                  </span>
                );
              })}
            </div>
          </div>
        ))}
      </div>

      <p className="hint" style={{ marginTop: 14 }}>
        The guide matches the target across the window; highlighted columns are where each spared
        allele differs. HUSH scores a window by its <b>weakest</b> spared allele
        {weakest ? <> &mdash; here {weakest.name}</> : null}.
        {!model.complete && (
          <>
            {" "}
            Showing each allele&apos;s deciding mismatch; the full per-base window alignment needs
            the panel surfaced by the API (planned).
          </>
        )}
      </p>

      {tip.show && (
        <div className="tip show" role="tooltip" style={{ left: tip.x, top: tip.y }}>
          <div className="tt">{tip.title}</div>
          {tip.lines.map((l, i) => (
            <div key={i}>{l}</div>
          ))}
        </div>
      )}
    </div>
  );
}
