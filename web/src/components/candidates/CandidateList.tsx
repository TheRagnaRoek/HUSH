import type { CandidateResult } from "../../api/types";
import { discToken } from "../../lib/disc";
import { useSession } from "../../session/useSession";

function discLabel(c: CandidateResult): string {
  const d = c.allele?.discrimination;
  return d ? discToken(d).split(" ")[0] : "-";
}

export function CandidateList() {
  const { design, selectedRegion, setSelectedRegion } = useSession();
  const selected = design?.selected ?? [];
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
            region {c.region} &middot; guide 5&prime;&rarr;3&prime;
          </span>
          <span className="badges">
            <span className="mini disc">disc {discLabel(c)}</span>
            {c.allele?.spare_min_mismatch != null && (
              <span className="mini gc">minDiff {c.allele.spare_min_mismatch}</span>
            )}
          </span>
        </button>
      ))}
    </div>
  );
}
