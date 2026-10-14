/**
 * The off-target gate's owned-transcript set, reduced to what the UI needs.
 * Kept apart from SelfTranscripts.tsx so that file exports components only.
 * See docs/design-rules.md "What the design owns".
 */
import type { CriteriaSummary } from "../../api/types";

export interface OwnedSet {
  /** Normalised to upper case and sorted, as the CLI prints them. */
  accessions: string[];
  count: number;
  declared: number;
  resolved: number;
  coveragePercent: number | null;
  empty: boolean;
}

/** The owned set for a run, or null when the gate is off and it has no meaning. */
export function ownedSet(criteria: CriteriaSummary): OwnedSet | null {
  if (!criteria.offtarget_filter) return null;
  const accessions = [...(criteria.offtarget_self_accessions ?? [])]
    .map((a) => a.toUpperCase())
    .sort();
  const coverage = criteria.offtarget_self_coverage;
  return {
    accessions,
    count: accessions.length,
    declared: criteria.offtarget_self_declared ?? 0,
    resolved: criteria.offtarget_self_resolved ?? 0,
    coveragePercent: coverage != null ? Math.round(coverage * 100) : null,
    empty: accessions.length === 0,
  };
}
