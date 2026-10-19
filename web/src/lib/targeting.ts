/**
 * How a design was ranked, from its targeting summary. The API leaves
 * targeting empty when no panel was loaded, and marks a pan-allele request
 * (nothing spared) as not discriminating. Both are ordered by target position.
 */

import type { TargetingSummary } from "../api/types";

export type RankMode = "none" | "pan" | "discriminating";

export function rankMode(targeting: TargetingSummary | null | undefined): RankMode {
  if (!targeting) return "none";
  return targeting.discriminating ? "discriminating" : "pan";
}

/** Why there is nothing to compare against, for the modes that have no alleles. */
export const NO_ALLELES: Record<Exclude<RankMode, "discriminating">, string> = {
  none: "No panel is loaded, so there are no spared alleles to compare against.",
  pan: "Pan-allele targeting spares nothing, so there are no alleles to compare against.",
};

/** The ranking note shown above the results for the modes that have no discrimination. */
export const POSITION_ORDER: Record<Exclude<RankMode, "discriminating">, string> = {
  none: "No panel is loaded, so candidates are ordered by position along the target and carry no allele discrimination. Add a panel to rank them by how well they tell alleles apart.",
  pan: "Pan-allele targeting spares nothing, so candidates are ordered by position along the target and carry no allele discrimination.",
};
