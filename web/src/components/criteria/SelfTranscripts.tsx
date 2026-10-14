/**
 * The off-target gate's owned-transcript set, as the GUI reports it. Mirrors
 * what the CLI prints: one line on what was excluded and how it was found, and
 * a loud warning when an active filter owns nothing. See docs/design-rules.md
 * "What the design owns" for why an empty set is a broken run, not a clean one.
 */
import { useState } from "react";
import type { CriteriaSummary } from "../../api/types";
import { CopyButton } from "../common/CopyButton";
import { ownedSet, type OwnedSet } from "./ownedSet";

/** How the set was found: named by hand, resolved by homology, or both. */
function provenance(own: OwnedSet): string {
  const parts: string[] = [];
  if (own.declared > 0) parts.push(`${own.declared} named`);
  if (own.resolved > 0) {
    parts.push(
      own.coveragePercent != null
        ? `${own.resolved} by homology at ${own.coveragePercent}% k-mers`
        : `${own.resolved} by homology`,
    );
  }
  return parts.join("; ");
}

/**
 * The owned set for the last run. The accession list sits behind a disclosure
 * because homology resolution on a gene family returns well over a hundred.
 */
export function OwnedTranscripts({ criteria }: { criteria: CriteriaSummary }) {
  const [open, setOpen] = useState(false);
  const own = ownedSet(criteria);
  if (!own) return null;

  const how = provenance(own);
  return (
    <div className="owned">
      <div className="owned-row">
        <span className="owned-k">Own transcripts</span>
        {own.empty ? (
          <span className="owned-v empty">
            None excluded - every reference hit counts against a candidate.
          </span>
        ) : (
          <span className="owned-v">
            <b>{own.count}</b> excluded from the off-target gate
            {how ? ` (${how})` : ""}.
          </span>
        )}
        {!own.empty && (
          <button
            type="button"
            className="linkbtn"
            aria-expanded={open}
            onClick={() => setOpen((v) => !v)}
          >
            {open ? "Hide accessions" : `Show ${own.count} accessions`}
          </button>
        )}
      </div>
      {open && !own.empty && (
        <div className="owned-list">
          <div className="owned-acc mono">{own.accessions.join(", ")}</div>
          <CopyButton text={own.accessions.join(",")} label="Copy accessions" />
        </div>
      )}
    </div>
  );
}

/**
 * The GUI half of the CLI's NO SELF EXCLUDED warning. An active gate that owns
 * nothing measures on-target homology, so it reads as "no candidate is any
 * good" when the designs are fine.
 */
export function NoSelfExcludedCallout({
  criteria,
  onOpenCriteria,
}: {
  criteria: CriteriaSummary;
  onOpenCriteria?: () => void;
}) {
  const own = ownedSet(criteria);
  if (!own || !own.empty) return null;

  return (
    <div className="callout bad">
      <span className="ic">{"\u26A0"}</span>
      <span>
        <b>
          No self excluded - the off-target gate counted this design&apos;s own
          transcript against every candidate.
        </b>{" "}
        Nothing identified the transcripts the design owns, so the gate measured
        on-target homology rather than off-target risk. On a gene family that
        rejects nearly every candidate, which looks like a quality problem when
        there is nothing wrong with the designs.
        <br />
        Keep <b>resolve own</b> on in the off-target gate and lower its k-mer
        coverage if the family is diverse, then run again. The CLI can also name
        them outright with --offtarget-self-accession.
        {onOpenCriteria && (
          <>
            <br />
            <button type="button" className="mini-btn" onClick={onOpenCriteria}>
              Open the off-target gate
            </button>
          </>
        )}
      </span>
    </div>
  );
}
