/**
 * Decode a DiscriminationResult into the compact disc token and plain language.
 * Fields are explained in docs/design-rules.md (Discrimination ranking).
 */

import type { DiscriminationResult } from "../api/types";

export const REGION: Record<string, string> = {
  seed: "in the seed (2-8) - can be overwhelmed by high siRNA levels",
  post: "3' of the seed - discriminates robustly",
  anch: "position 1 anchor - unpaired",
};

export const HUANG: Record<string, string> = {
  "H+": "listed by Huang et al. for this exact base change",
  "H-": "not listed by Huang et al.",
  "H?": "ambiguous bases (no table row)",
};

export const TIER: Record<string, string> = {
  S: "a strongest position (P10 / P16)",
  C: "a merely consistent position",
  "-": "neither tier",
};

export const CLS: Record<string, string> = {
  pp: "purine:purine - discriminates best",
  pu: "purine:pyrimidine - milder",
  up: "pyrimidine:purine - milder",
  uu: "pyrimidine:pyrimidine - mildest",
  "?": "ambiguous",
};

export function huangOf(d: DiscriminationResult): "H+" | "H-" | "H?" {
  if (d.huang_recommended === true) return "H+";
  if (d.huang_recommended === false) return "H-";
  return "H?";
}

export function tierOf(d: DiscriminationResult): "S" | "C" | "-" {
  if (d.strongest_position) return "S";
  if (d.consistent_position) return "C";
  return "-";
}

export function classOf(d: DiscriminationResult): string {
  return d.mismatch_class ?? "?";
}

/** The one-line token, e.g. "P16 post H+ S pp". */
export function discToken(d: DiscriminationResult): string {
  return `P${d.guide_position} ${d.region} ${huangOf(d)} ${tierOf(d)} ${classOf(d)}`;
}

export interface DiscChip {
  label: string;
  value: string;
  emphasis?: boolean;
}

export function discChips(d: DiscriminationResult): DiscChip[] {
  return [
    { label: "pos", value: `P${d.guide_position}`, emphasis: true },
    { label: "region", value: d.region },
    { label: "huang", value: huangOf(d) },
    { label: "tier", value: tierOf(d) },
    { label: "class", value: classOf(d) },
  ];
}

/** Prose lines for a tooltip on a highlighted mismatch base. */
export function discProse(d: DiscriminationResult): string[] {
  return [
    `region: ${REGION[d.region] ?? d.region}`,
    `Huang: ${HUANG[huangOf(d)]}`,
    `tier: ${TIER[tierOf(d)]}`,
    `class: ${CLS[classOf(d)]}`,
  ];
}
