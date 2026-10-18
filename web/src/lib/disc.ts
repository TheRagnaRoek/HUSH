/**
 * Decode a DiscriminationResult into the compact disc token and plain language.
 * Fields are explained in docs/design-rules.md (Discrimination ranking).
 */

import type { DiscriminationResult } from "../api/types";
import type { SourceKey } from "../reference/sources";

// The API sends enum values ("post_seed", "purine:purine"); the CLI and the
// sample use short codes. Everything here works on the short codes.
const REGION_CODE: Record<string, string> = {
  anchor: "anch",
  anch: "anch",
  seed: "seed",
  post_seed: "post",
  post: "post",
};

const CLASS_CODE: Record<string, string> = {
  "purine:purine": "pp",
  "purine:pyrimidine": "pu",
  "pyrimidine:purine": "up",
  "pyrimidine:pyrimidine": "uu",
  pp: "pp",
  pu: "pu",
  up: "up",
  uu: "uu",
};

export const REGION: Record<string, string> = {
  seed: "in the seed (positions 2-8) - can be overwhelmed by high siRNA levels",
  post: "3' of the seed (position 9 or later) - discriminates robustly",
  anch: "at position 1, the anchor - unpaired",
};

export const HUANG: Record<string, string> = {
  "H+": "listed by Huang et al. for this exact base change",
  "H-": "not listed by Huang et al.",
  "H?": "ambiguous bases (no table row)",
};

export const TIER: Record<string, string> = {
  S: "one of the best-replicated positions (P10, P13, P16)",
  W: "inside the P9-P16 window where mismatches are least tolerated",
  "-": "outside the P9-P16 window",
};

export const CLS: Record<string, string> = {
  pp: "purine:purine - discriminates best",
  pu: "purine:pyrimidine - milder",
  up: "pyrimidine:purine - milder",
  uu: "pyrimidine:pyrimidine - mildest",
  "?": "ambiguous",
};

export function regionOf(d: DiscriminationResult): string {
  return REGION_CODE[d.region] ?? d.region;
}

export function huangOf(d: DiscriminationResult): "H+" | "H-" | "H?" {
  if (d.huang_recommended === true) return "H+";
  if (d.huang_recommended === false) return "H-";
  return "H?";
}

/** S = a best-replicated position, W = inside the P9-P16 window, - = outside. */
export function tierOf(d: DiscriminationResult): "S" | "W" | "-" {
  if (d.strongest_position) return "S";
  if (d.consistent_position) return "W";
  return "-";
}

export function classOf(d: DiscriminationResult): string {
  return (d.mismatch_class && CLASS_CODE[d.mismatch_class]) || "?";
}

/** The one-line token, e.g. "P16 post H+ S pp". */
export function discToken(d: DiscriminationResult): string {
  return `P${d.guide_position} ${regionOf(d)} ${huangOf(d)} ${tierOf(d)} ${classOf(d)}`;
}

export interface DiscChip {
  label: string;
  value: string;
  /** Hover text: what the qualifier is, then what this value means. */
  info: string;
  /** Sources named in the hover tip; every qualifier links to the ranking section. */
  cites: SourceKey[];
  emphasis?: boolean;
}

export function discChips(d: DiscriminationResult): DiscChip[] {
  const region = regionOf(d);
  return [
    {
      label: "Position",
      value: `P${d.guide_position}`,
      info: "Where the deciding mismatch sits on the guide, counted from its 5' end. The position sets the region and the tier.",
      cites: ["schwarz2006", "du2005", "huang2009"],
      emphasis: true,
    },
    {
      label: "Region",
      value: region,
      info: `Which part of the guide carries the mismatch. This one is ${REGION[region] ?? region}.`,
      cites: ["schwarz2006", "du2005"],
    },
    {
      label: "Huang",
      value: huangOf(d),
      info: `Whether Huang et al. list this exact base change at this position. Here: ${HUANG[huangOf(d)].replace(/\.$/, "")}.`,
      cites: ["huang2009"],
    },
    {
      label: "Tier",
      value: tierOf(d),
      info: `How the position sits in the published data: S for a best-replicated position, W for the rest of the P9-P16 window, a dash outside it. Here: ${TIER[tierOf(d)]}.`,
      cites: ["schwarz2006", "du2005"],
    },
    {
      label: "Class",
      value: classOf(d),
      info: `The chemistry of the mismatched base pair. Here: ${CLS[classOf(d)]}.`,
      cites: ["schwarz2006", "du2005"],
    },
  ];
}

/** Prose lines for a tooltip on a highlighted mismatch base. */
export function discProse(d: DiscriminationResult): string[] {
  return [
    `Region: ${REGION[regionOf(d)] ?? regionOf(d)}`,
    `Huang: ${HUANG[huangOf(d)]}`,
    `Tier: ${TIER[tierOf(d)]}`,
    `Class: ${CLS[classOf(d)]}`,
    "Sources: Huang 2009; Schwarz 2006; Du 2005",
  ];
}

const CLASS_NAME: Record<string, string> = {
  pp: "purine:purine",
  pu: "purine:pyrimidine",
  up: "pyrimidine:purine",
  uu: "pyrimidine:pyrimidine",
};

const CLASS_NOTE: Record<string, string> = {
  pp: "the class that discriminates best",
  pu: "a milder class",
  up: "a milder class",
  uu: "the mildest class",
};

/** Full sentences explaining the deciding mismatch, for the reasoning tab. */
export function discReasoning(d: DiscriminationResult): string[] {
  const p = d.guide_position;
  const region = regionOf(d);
  const cls = classOf(d);
  const out: string[] = [];

  if (region === "seed") {
    out.push(
      `The deciding difference sits at guide position ${p}, inside the seed (positions 2-8). A seed mismatch blocks binding but can be overcome by a high siRNA concentration, so it discriminates less reliably.`,
    );
  } else if (region === "post") {
    out.push(
      `The deciding difference sits at guide position ${p}, 3' of the seed. A mismatch there blocks cleavage of the target rather than binding, so it discriminates robustly.`,
    );
  } else {
    out.push(
      `The deciding difference sits at guide position ${p}, the anchor, which does not pair with the target.`,
    );
  }

  // The API writes the change as "spared>target" in DNA letters.
  const change = /([A-Za-z])\s*>\s*([A-Za-z])/.exec(d.mutation);
  if (CLASS_NAME[cls]) {
    out.push(
      (change
        ? `The spared allele has ${change[1].toUpperCase()} where the target has ${change[2].toUpperCase()}, so the guide meets it as a`
        : "The guide meets the spared allele as a") +
        ` ${CLASS_NAME[cls]} mismatch, ${CLASS_NOTE[cls]}.`,
    );
  } else {
    out.push("One of the bases is ambiguous, so the mismatch has no class.");
  }

  const huang = huangOf(d);
  out.push(
    huang === "H+"
      ? "Huang et al. list this exact base change at this position."
      : huang === "H-"
        ? "Huang et al. do not list this exact base change at this position."
        : "Huang et al. have no entry for ambiguous bases.",
  );

  const tier = tierOf(d);
  out.push(
    tier === "S"
      ? `Position ${p} is one of the best-replicated positions (P10, P13, P16), inside the P9-P16 window where mismatches are least tolerated.`
      : tier === "W"
        ? `Position ${p} lies inside the P9-P16 window where mismatches are least tolerated, but it is not one of the best-replicated positions (P10, P13, P16).`
        : `Position ${p} lies outside the P9-P16 window where mismatches are least tolerated.`,
  );

  return out;
}
