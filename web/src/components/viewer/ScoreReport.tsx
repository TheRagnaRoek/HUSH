import type { CandidateResult, CriteriaSummary } from "../../api/types";
import { CLS, discToken } from "../../lib/disc";
import { ownedSet } from "../criteria/ownedSet";

const RULE_NAMES = ["Ui-Tei", "Reynolds", "Amarzguioui"];

function passedRules(code: string): string {
  const names = [...code].map((ch, i) => (ch !== "-" ? RULE_NAMES[i] : null)).filter(Boolean);
  return names.length ? names.join(", ") : "none";
}

export function ScoreReport({
  candidate,
  criteria,
  offtargetSkipped,
  rank,
}: {
  candidate: CandidateResult;
  criteria: CriteriaSummary;
  offtargetSkipped: boolean;
  rank: number;
}) {
  const weakest = candidate.allele?.per_allele?.[0];
  const d = weakest?.discrimination ?? candidate.allele?.discrimination ?? null;
  const gcMaxPct = (criteria.gc_max ?? 0.52) * 100;
  const gcOk = candidate.gc_percent <= gcMaxPct;
  // An active gate that owns nothing counted the target gene against this
  // candidate, so its mismatch verdict is not an off-target verdict.
  const own = ownedSet(criteria);
  const offtargetDot = offtargetSkipped ? "info" : own?.empty ? "warn" : "pass";

  return (
    <div className="vpane">
      <div className="why-lead">
        <b>Why this candidate ranks #{rank}.</b> It clears every enabled gate and, among the
        survivors, best distinguishes the allele it separates least well
        {weakest ? <> ({weakest.allele})</> : null}
        {d ? (
          <>
            {" "}
            - a {d.region === "post" ? "robust 3′-of-seed" : "seed-region"}{" "}
            {CLS[candidate.allele?.discrimination?.mismatch_class ?? "?"]?.split(" -")[0] ??
              "mismatch"}{" "}
            change at guide position {d.guide_position}.
          </>
        ) : (
          "."
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
              : "Filter off for this GC-rich family (documented HLA behaviour)."}
          </span>
        </div>

        <div className="srow">
          <span className={`vd ${offtargetDot}`} />
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
                {own &&
                  (own.empty ? (
                    <>
                      {" "}
                      <b>No transcript was excluded as this design&apos;s own</b>, so the
                      target gene itself is among those hits.
                    </>
                  ) : (
                    <>
                      {" "}
                      Excludes <b>{own.count}</b> own transcript
                      {own.count === 1 ? "" : "s"}.
                    </>
                  ))}
              </>
            )}
          </span>
        </div>

        <div className="srow">
          <span className={`vd ${gcOk ? "pass" : "warn"}`} />
          <span className="sk">GC content</span>
          <span className="sv">
            <b>{candidate.gc_percent}%</b> -{" "}
            {gcOk ? "within the window" : `above the ${gcMaxPct}% window, expected on HLA`}.
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

        {candidate.allele && (
          <div className="srow">
            <span className="vd info" />
            <span className="sk">Allele discrimination</span>
            <span className="sv">
              Weakest spared allele still carries <b>{candidate.allele.spare_min_mismatch}</b>{" "}
              mismatches
              {candidate.allele.discrimination && (
                <>
                  ; deciding change graded{" "}
                  <span className="mono">{discToken(candidate.allele.discrimination)}</span>
                </>
              )}
              .
            </span>
          </div>
        )}
      </div>
    </div>
  );
}
