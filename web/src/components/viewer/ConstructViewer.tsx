import { useEffect, useState } from "react";
import type { CandidateResult, ConstructResult, CriteriaSummary } from "../../api/types";
import { buildAlignment, type AlleleEnrichment } from "../../lib/alignment";
import type { RankMode } from "../../lib/targeting";
import type { ExportFormat } from "../../lib/tsv";
import { cap } from "../../lib/text";
import { AlignmentTrack } from "./AlignmentTrack";
import { ConstructSchematic } from "./ConstructSchematic";
import { DiscriminationBreakdown } from "./DiscriminationBreakdown";
import { OligoRow } from "./OligoRow";
import { QcChips } from "./QcChips";
import { ScoreReport } from "./ScoreReport";
import { SirnaDuplex } from "./SirnaDuplex";
import { WithheldNotice } from "./WithheldNotice";

type Tab = "align" | "disc" | "score" | "build";

/** "overview" is the short export summary on the candidates step; "full" is
 * the whole construct view on the construct step. */
type BuildView = "overview" | "full";

const TABS: { id: Tab; label: string }[] = [
  { id: "align", label: "Target & alignment" },
  { id: "disc", label: "Discrimination" },
  { id: "score", label: "Reasoning" },
  { id: "build", label: "shRNA construct" },
];

function qcIssues(c: ConstructResult): number {
  return [
    c.folds_as_designed,
    c.starts_with_purine,
    c.premature_terminator_free,
    c.max_homopolymer_run <= 3,
  ].filter((ok) => !ok).length;
}

export function ConstructViewer({
  candidate,
  construct,
  withheld,
  enrichment,
  criteria,
  offtargetSkipped,
  rank,
  mode = "discriminating",
  targetId,
  terminator = "UUUUUU",
  startTab = "align",
  buildView = "full",
  keepTab = false,
  onOpenConstruct,
  onExport,
}: {
  candidate: CandidateResult;
  construct?: ConstructResult;
  /** The construct the engine withheld for this candidate, if it did. */
  withheld?: ConstructResult;
  enrichment?: AlleleEnrichment[];
  criteria: CriteriaSummary;
  offtargetSkipped: boolean;
  rank: number;
  /** How the design was ranked; without discrimination the allele views explain why. */
  mode?: RankMode;
  targetId: string;
  terminator?: string;
  startTab?: Tab;
  buildView?: BuildView;
  keepTab?: boolean;
  onOpenConstruct?: () => void;
  /** Saves the selected construct in the chosen format; the buttons are disabled without it. */
  onExport?: (format: ExportFormat) => void;
}) {
  const [tab, setTab] = useState<Tab>(startTab);
  // keepTab: stay on the current tab when another candidate is picked.
  useEffect(() => {
    if (!keepTab) setTab(startTab);
  }, [startTab, candidate.region, keepTab]);
  const overview = buildView === "overview";

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
          Region {candidate.region} &middot; Rank {rank}
        </span>
        <span className="spacer" />
        <span className="mini disc">
          {mode === "discriminating"
            ? `Discrimination rank #${rank}`
            : `Target position order #${rank}`}
        </span>
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
            {t.id === "build" && overview ? "Export" : t.label}
          </button>
        ))}
      </div>

      {tab === "align" && <AlignmentTrack model={model} targetId={targetId} mode={mode} />}
      {tab === "disc" && <DiscriminationBreakdown candidate={candidate} mode={mode} />}
      {tab === "score" && (
        <ScoreReport
          candidate={candidate}
          criteria={criteria}
          offtargetSkipped={offtargetSkipped}
          rank={rank}
          mode={mode}
        />
      )}
      {tab === "build" && overview && (
        <div className="vpane">
          <SirnaDuplex candidate={candidate} />
          <div className="divider" />
          <div className="subhead">
            shRNA hairpin <span className="note">- Optional; for a vector-expressed hairpin</span>
          </div>
          {construct ? (
            <>
              <dl className="kv exportkv">
                <dt>Loop</dt>
                <dd>
                  {cap(construct.loop_name)} <span className="mono">({construct.loop_seq})</span>
                </dd>
                <dt>Hairpin</dt>
                <dd className="num">{construct.hairpin_rna_5to3.length} nt</dd>
                <dt>Construct checks</dt>
                <dd className={qcIssues(construct) === 0 ? "ready" : "attn"}>
                  {qcIssues(construct) === 0
                    ? "All passed"
                    : `${qcIssues(construct)} ${qcIssues(construct) === 1 ? "check" : "checks"} to review`}
                </dd>
              </dl>
              <div className="subhead">
                Cloning oligos{" "}
                <span className="note">- Vector-agnostic; add your restriction sticky-end tails</span>
              </div>
              <OligoRow label={"Insert top 5′→3′"} seq={construct.insert_top_dna_5to3} />
              <OligoRow label={"Insert bottom 5′→3′"} seq={construct.insert_bottom_dna_5to3} />
            </>
          ) : withheld ? (
            <WithheldNotice construct={withheld} terminator={terminator} />
          ) : (
            <p className="hint">
              Not converted. The siRNA duplex above needs no conversion; convert it in the
              construct step if you want a hairpin to express from a vector.
            </p>
          )}
          <div className="actions">
            {construct && (
              <>
                <button
                  type="button"
                  className="hbtn ghost"
                  disabled={!onExport}
                  onClick={() => onExport?.("tsv")}
                >
                  Export construct (TSV)
                </button>
                <button
                  type="button"
                  className="hbtn ghost"
                  disabled={!onExport}
                  onClick={() => onExport?.("txt")}
                >
                  Export construct (TXT)
                </button>
              </>
            )}
            {onOpenConstruct && (
              <button type="button" className="hbtn" onClick={onOpenConstruct}>
                {construct || withheld ? "Open construct view" : "Convert to shRNA"} &rarr;
              </button>
            )}
          </div>
        </div>
      )}
      {tab === "build" && !overview && (
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
              <div className="buildblock">
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
              <div className="divider" />
              <div className="subhead">
                Cloning oligos{" "}
                <span className="note">- Vector-agnostic; add your restriction sticky-end tails</span>
              </div>
              <OligoRow label={"Insert top 5′→3′"} seq={construct.insert_top_dna_5to3} />
              <OligoRow
                label={"Insert bottom 5′→3′"}
                seq={construct.insert_bottom_dna_5to3}
              />
              <div className="actions">
                {onOpenConstruct && (
                  <button type="button" className="hbtn" onClick={onOpenConstruct}>
                    Open full construct view
                  </button>
                )}
                <button
                  type="button"
                  className="hbtn ghost"
                  disabled={!onExport}
                  onClick={() => onExport?.("tsv")}
                >
                  Export construct (TSV)
                </button>
                <button
                  type="button"
                  className="hbtn ghost"
                  disabled={!onExport}
                  onClick={() => onExport?.("txt")}
                >
                  Export construct (TXT)
                </button>
              </div>
            </>
          ) : (
            <>
              {withheld ? (
                <WithheldNotice construct={withheld} terminator={terminator} />
              ) : (
                <p className="hint" style={{ marginBottom: 14 }}>
                  No shRNA construct yet. Converting is optional: it builds a hairpin to express
                  from a vector, and the siRNA duplex below needs no conversion. Use the Convert
                  button to make one.
                </p>
              )}
              <SirnaDuplex candidate={candidate} />
            </>
          )}
        </div>
      )}
    </div>
  );
}
