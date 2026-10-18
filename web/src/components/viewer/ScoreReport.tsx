import type { CandidateResult, CriteriaSummary } from "../../api/types";
import { CLS, classOf, discReasoning, regionOf } from "../../lib/disc";
import { NO_ALLELES, type RankMode } from "../../lib/targeting";

const RULE_NAMES = ["Ui-Tei", "Reynolds", "Amarzguioui"];

function passedRules(code: string): string {
  const names = [...code].map((ch, i) => (ch !== "-" ? RULE_NAMES[i] : null)).filter(Boolean);
  return names.length ? names.join(", ") : "none";
}

function ordinal(n: number): string {
  const tens = n % 100;
  if (tens >= 11 && tens <= 13) return `${n}th`;
  return `${n}${["th", "st", "nd", "rd"][n % 10 > 3 ? 0 : n % 10]}`;
}

export function ScoreReport({
  candidate,
  criteria,
  offtargetSkipped,
  rank,
  mode = "discriminating",
}: {
  candidate: CandidateResult;
  criteria: CriteriaSummary;
  offtargetSkipped: boolean;
  rank: number;
  mode?: RankMode;
}) {
  const weakest = candidate.allele?.per_allele?.[0];
  const d = weakest?.discrimination ?? candidate.allele?.discrimination ?? null;
  const mismatchCount = candidate.allele?.spare_min_mismatch ?? weakest?.diff_count;
  const gcMaxPct = (criteria.gc_max ?? 0.52) * 100;
  const gcOk = candidate.gc_percent <= gcMaxPct;

  return (
    <div className="vpane">
      <div className="why-lead">
        {mode === "discriminating" ? (
          <>
            <b>Why this candidate ranks #{rank}.</b> It clears every enabled gate and, among the
            survivors, best distinguishes the allele it separates least well
            {weakest ? <> ({weakest.allele})</> : null}
            {d ? (
              <>
                {" "}
                - a {regionOf(d) === "post" ? "robust 3′-of-seed" : "seed-region"}{" "}
                {CLS[classOf(d)]?.split(" -")[0] ?? "mismatch"} change at guide position{" "}
                {d.guide_position}.
              </>
            ) : (
              "."
            )}
          </>
        ) : (
          <>
            <b>Why this candidate is #{rank}.</b> It clears every enabled gate. With{" "}
            {mode === "none" ? "no panel" : "nothing spared"} there is no allele to discriminate
            against, so candidates are ordered by position along the target: this is the{" "}
            {ordinal(rank)} window to pass, counting from the start of the sequence.
          </>
        )}
      </div>

      <div className="report">
        <div className="srow">
          <span className="vd pass" />
          <span className="sk">Classical rules</span>
          <span className="sv">
            Passes <b>{passedRules(candidate.rules_passed)}</b> - the union of Ui-Tei / Reynolds /
            Amarzguioui.
          </span>
        </div>

        <div className="srow">
          <span className={`vd ${criteria.seed_filter ? "pass" : "info"}`} />
          <span className="sk">Seed off-target</span>
          <span className="sv">
            {criteria.seed_filter
              ? `Seed-duplex Tm within ${criteria.seed_tm_threshold_c ?? 21.5} C.`
              : "Filter switched off for this design."}
          </span>
        </div>

        <div className="srow">
          <span className={`vd ${offtargetSkipped ? "info" : "pass"}`} />
          <span className="sk">Transcriptome off-target</span>
          <span className="sv">
            {offtargetSkipped ? (
              "Skipped - no reference index configured."
            ) : (
              <>
                Nearest off-target <b>&ge; {candidate.offtarget_min_mismatch} mismatches</b>.
                {candidate.offtarget_hits && (
                  <>
                    {" "}
                    Hits at 0/1/2/3&nbsp;mm:{" "}
                    <span className="mono">{candidate.offtarget_hits.join(" / ")}</span>.
                  </>
                )}
              </>
            )}
          </span>
        </div>

        <div className="srow">
          <span className={`vd ${!criteria.gc_filter || gcOk ? "pass" : "warn"}`} />
          <span className="sk">GC content</span>
          <span className="sv">
            <b>{candidate.gc_percent.toFixed(1)}%</b>
            {criteria.gc_filter
              ? ` - ${gcOk ? "within" : "above"} the ${gcMaxPct}% window.`
              : " - no GC window applied."}
          </span>
        </div>

        <div className="srow">
          <span className="vd pass" />
          <span className="sk">Homopolymer</span>
          <span className="sv">
            Longest run <b>{candidate.homopolymer_max_run ?? "-"}</b> &le; cap{" "}
            {criteria.max_homopolymer_run ?? 3}.
          </span>
        </div>

        {!candidate.allele && mode !== "discriminating" && (
          <div className="srow">
            <span className="vd info" />
            <span className="sk">Allele discrimination</span>
            <span className="sv">
              Not assessed. {NO_ALLELES[mode]}
              {mode === "none"
                ? " Load a panel and choose something to spare to compare alleles."
                : " Choose an allele-specific, target-specific or custom request to compare alleles."}
            </span>
          </div>
        )}

        {candidate.allele && (
          <div className="srow">
            <span className="vd info" />
            <span className="sk">Allele discrimination</span>
            <div className="sv">
              <p>
                The weakest spared allele{weakest ? <> ({weakest.allele})</> : null} still carries{" "}
                <b>
                  {mismatchCount ?? "-"} {mismatchCount === 1 ? "mismatch" : "mismatches"}
                </b>{" "}
                against this guide. The window is ranked by how well its best-placed difference
                from that allele discriminates.
              </p>
              {d && <p>{discReasoning(d).join(" ")}</p>}
            </div>
          </div>
        )}
      </div>
    </div>
  );
}
