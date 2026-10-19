import type { CandidateResult } from "../../api/types";
import { OligoRow } from "./OligoRow";

/** The synthetic siRNA duplex of a candidate; needs no shRNA conversion. */
export function SirnaDuplex({ candidate }: { candidate: CandidateResult }) {
  return (
    <>
      <div className="subhead">
        siRNA duplex <span className="note">- Synthetic strands, 5′→3′, 3′ overhangs included</span>
      </div>
      <OligoRow label="Guide (antisense)" seq={candidate.guide_5to3} />
      <OligoRow label="Passenger (sense)" seq={candidate.passenger_5to3} />
    </>
  );
}
