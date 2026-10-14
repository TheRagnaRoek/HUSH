import { useEffect, useState } from "react";
import type { CandidateResult, ConstructResult, CriteriaSummary } from "../../api/types";
import { buildAlignment, type AlleleEnrichment } from "../../lib/alignment";
import { AlignmentTrack } from "./AlignmentTrack";
import { ConstructSchematic } from "./ConstructSchematic";
import { DiscriminationBreakdown } from "./DiscriminationBreakdown";
import { HairpinFold } from "./HairpinFold";
import { OligoRow } from "./OligoRow";
import { QcChips } from "./QcChips";
import { ScoreReport } from "./ScoreReport";

type Tab = "align" | "disc" | "score" | "build";

const TABS: { id: Tab; label: string }[] = [
  { id: "align", label: "Target & alignment" },
  { id: "disc", label: "Discrimination" },
  { id: "score", label: "Why this one" },
  { id: "build", label: "shRNA construct" },
];

export function ConstructViewer({
  candidate,
  construct,
  enrichment,
  criteria,
  offtargetSkipped,
  rank,
  targetId,
  terminator = "UUUUUU",
  startTab = "align",
  onOpenConstruct,
}: {
  candidate: CandidateResult;
  construct?: ConstructResult;
  enrichment?: AlleleEnrichment[];
  criteria: CriteriaSummary;
  offtargetSkipped: boolean;
  rank: number;
  targetId: string;
  terminator?: string;
  startTab?: Tab;
  onOpenConstruct?: () => void;
}) {
  const [tab, setTab] = useState<Tab>(startTab);
  useEffect(() => setTab(startTab), [startTab, candidate.region]);

  const model = buildAlignment(candidate, enrichment);
  // The hairpin encodes both strands in full, so the schematic shows the
  // cores and their overhangs as separate segments rather than hiding them.
  const sense = candidate.passenger_5to3.slice(
    0,
    candidate.passenger_5to3.length - candidate.passenger_overhang.length,
  );
  const antisense = candidate.guide_5to3.slice(
    0,
    candidate.guide_5to3.length - candidate.guide_overhang.length,
  );

  return (
    <div>
      <div className="vh">
        <span className="title">{candidate.guide_5to3}</span>
        <span className="region">
          region {candidate.region} &middot; rank {rank}
        </span>
        <span className="spacer" />
        <span className="mini disc">discrimination rank #{rank}</span>
      </div>

      <div className="segtabs" role="tablist">
        {TABS.map((t) => (
          <button
            key={t.id}
            type="button"
            className="tab"
            aria-selected={t.id === tab}
            onClick={() => setTab(t.id)}
          >
            {t.label}
          </button>
        ))}
      </div>

      {tab === "align" && <AlignmentTrack model={model} targetId={targetId} />}
      {tab === "disc" && <DiscriminationBreakdown candidate={candidate} />}
      {tab === "score" && (
        <ScoreReport
          candidate={candidate}
          criteria={criteria}
          offtargetSkipped={offtargetSkipped}
          rank={rank}
        />
      )}
      {tab === "build" && (
        <div className="vpane">
          {construct ? (
            <>
              <ConstructSchematic
                sense={sense}
                senseOverhang={construct.sense_overhang}
                loop={construct.loop_seq}
                antisense={antisense}
                guideOverhang={construct.guide_overhang}
                terminator={terminator}
              />
              <div className="buildgrid">
                <HairpinFold construct={construct} />
                <div>
                  <div className="subhead">Fold QC</div>
                  <QcChips construct={construct} />
                  {!construct.starts_with_purine && (
                    <p className="hint" style={{ marginTop: 10 }}>
                      Not initiating on A or G - a human U6 or 7SK promoter starts most
                      efficiently on a purine at position 1, so prefer a candidate
                      starting with A or G.
                    </p>
                  )}
                </div>
              </div>
              <div className="divider" />
              <div className="subhead">
                Cloning oligos{" "}
                <span className="note">- vector-agnostic; add your restriction sticky-end tails</span>
              </div>
              <OligoRow label={"insert top 5′→3′"} seq={construct.insert_top_dna_5to3} />
              <OligoRow
                label={"insert bottom 5′→3′"}
                seq={construct.insert_bottom_dna_5to3}
              />
              <div className="actions">
                {onOpenConstruct && (
                  <button type="button" className="hbtn" onClick={onOpenConstruct}>
                    Open full construct view
                  </button>
                )}
                <button type="button" className="hbtn ghost">
                  Export construct (TSV)
                </button>
              </div>
            </>
          ) : (
            <p className="hint">
              No construct yet for this candidate. Convert the selected design to generate the
              shRNA scaffold.
            </p>
          )}
        </div>
      )}
    </div>
  );
}
