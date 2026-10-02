import type { CandidateResult } from "../../api/types";
import { discChips } from "../../lib/disc";

export function DiscriminationBreakdown({ candidate }: { candidate: CandidateResult }) {
  const perAllele = candidate.allele?.per_allele ?? [];
  return (
    <div className="vpane">
      <p className="hint" style={{ marginBottom: 4 }}>
        Every spared allele, weakest first. Each mismatch is graded by where it sits on the guide
        and what kind of base change it is - that grading sets the discrimination rank.
      </p>
      {perAllele.map((a, i) => {
        const weak = i === 0;
        return (
          <div className={`pa${weak ? " weak" : ""}`} key={a.allele}>
            <div className="pah">
              <span className="an">{a.allele}</span>
              <span className="df">
                {a.diff_count} difference{a.diff_count > 1 ? "s" : ""} in window
              </span>
              {weak && <span className="weakflag">deciding allele</span>}
            </div>
            {a.discrimination && (
              <div className="discstr">
                {discChips(a.discrimination).map((c) => (
                  <span className={`dtok${c.emphasis ? " pos" : ""}`} key={c.label}>
                    <b>{c.label}</b>
                    {c.value}
                  </span>
                ))}
              </div>
            )}
          </div>
        );
      })}
    </div>
  );
}
