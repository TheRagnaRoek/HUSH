import type { CandidateResult } from "../../api/types";
import { discChips } from "../../lib/disc";
import { NO_ALLELES, type RankMode } from "../../lib/targeting";
import { RefNote } from "../common/RefLink";
import { cap } from "../../lib/text";

export function DiscriminationBreakdown({
  candidate,
  mode = "discriminating",
}: {
  candidate: CandidateResult;
  mode?: RankMode;
}) {
  const perAllele = candidate.allele?.per_allele ?? [];
  if (perAllele.length === 0) {
    return (
      <div className="vpane">
        <p className="hint">
          {mode === "discriminating" ? NO_ALLELES.none : NO_ALLELES[mode]} This tab grades how well
          a candidate tells a spared allele apart, by where each mismatch sits on the guide and
          what base change it is, so it stays empty here.
          {mode === "none" &&
            " Load a panel and choose something to spare to rank candidates by discrimination."}
          {mode === "pan" &&
            " Choose an allele-specific, target-specific or custom request to rank candidates by discrimination."}
        </p>
      </div>
    );
  }
  return (
    <div className="vpane">
      <p className="hint" style={{ marginBottom: 4 }}>
        Every spared allele, weakest first. Each mismatch is graded by where it sits on the guide
        and what kind of base change it is - that grading sets the discrimination rank. Hover a
        qualifier to see what it means.
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
              {weak && <span className="weakflag">Deciding allele</span>}
            </div>
            {a.discrimination && (
              <div className="discstr">
                {discChips(a.discrimination).map((c) => (
                  <span className={`dtok has-tip${c.emphasis ? " pos" : ""}`} key={c.label} tabIndex={0}>
                    <b>{c.label}</b>
                    {c.label === "Region" ? cap(c.value) : c.value}
                    <span className="dtip" role="tooltip">
                      <span className="dtt">{c.label}</span>
                      {c.info}
                      <RefNote cites={c.cites} to="discrimination" />
                    </span>
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
