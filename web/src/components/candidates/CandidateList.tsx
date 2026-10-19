import type { CandidateResult } from "../../api/types";
import { discToken } from "../../lib/disc";
import { rankMode } from "../../lib/targeting";
import { useSession } from "../../session/useSession";

function discLabel(c: CandidateResult): string {
  const d = c.allele?.discrimination;
  return d ? discToken(d).split(" ")[0] : "-";
}

export function CandidateList() {
  const { design, convert, selectedRegion, setSelectedRegion } = useSession();
  const selected = design?.selected ?? [];
  const withheld = new Set(convert?.withheld.map((c) => c.region));
  const discriminating = rankMode(design?.targeting) === "discriminating";
  return (
    <div className="clist">
      {selected.map((c, i) => (
        <button
          key={c.region}
          type="button"
          className="cand"
          aria-selected={c.region === selectedRegion}
          onClick={() => setSelectedRegion(c.region)}
        >
          <span className="rank">{i + 1}</span>
          <span className="seq">{c.guide_5to3}</span>
          <span className="pos">
            Region {c.region} &middot; Guide 5&prime;&rarr;3&prime;
          </span>
          <span className="badges">
            {withheld.has(c.region) && <span className="mini warn">Withheld</span>}
            {discriminating ? (
              <>
                <span className="mini disc">Disc {discLabel(c)}</span>
                {c.allele?.spare_min_mismatch != null && (
                  <span className="mini gc">Min diff {c.allele.spare_min_mismatch}</span>
                )}
              </>
            ) : (
              <span className="mini gc">GC {c.gc_percent.toFixed(1)}%</span>
            )}
          </span>
        </button>
      ))}
    </div>
  );
}
