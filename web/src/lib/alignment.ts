/**
 * View model for the window alignment. The API returns only the deciding
 * mismatch and diff_count per spared allele, so an API-derived alignment is
 * marked complete: false. Sample data carries the full per-base window.
 *
 * TODO(api): return the aligned window or full mismatch list per panel member.
 */

import type { CandidateResult, DiscriminationResult } from "../api/types";

export interface AlnMismatch {
  /** 0-based column on the sense/target window. */
  col: number;
  /** The target base at that column. */
  from: string;
  /** The base the spared allele carries there (the one the row displays). */
  to: string;
  disc: DiscriminationResult;
}

export interface AlnAllele {
  name: string;
  diff: number;
  weak: boolean;
  mism: AlnMismatch[];
}

export interface AlignmentModel {
  target: string; // sense window core (= passenger core)
  guide: string; // guide core, 5'->3'
  alleles: AlnAllele[];
  complete: boolean; // false => only deciding mismatches are shown
}

/** Full per-base enrichment for an allele (sample data / post-API-addition). */
export interface AlleleEnrichment {
  name: string;
  diff: number;
  weak: boolean;
  mism: AlnMismatch[];
}

function core(seq: string, overhang: string): string {
  const n = overhang ? overhang.length : 0;
  return n > 0 ? seq.slice(0, seq.length - n) : seq;
}

/** Parse the API's "spared>target" mutation string into [target, spared]; tolerant of odd input. */
function parseMut(mutation: string): [string, string] {
  const m = /([A-Za-z])\s*[>/-]\s*([A-Za-z])/.exec(mutation);
  return m ? [m[2].toUpperCase(), m[1].toUpperCase()] : ["?", "?"];
}

export function buildAlignment(
  candidate: CandidateResult,
  enrichment?: AlleleEnrichment[],
): AlignmentModel {
  const target = core(candidate.passenger_5to3, candidate.passenger_overhang);
  const guide = core(candidate.guide_5to3, candidate.guide_overhang);
  const L = target.length;

  if (enrichment && enrichment.length > 0) {
    return { target, guide, alleles: enrichment, complete: true };
  }

  // API-derived: one deciding mismatch per spared allele, column inferred from
  // the guide position (guide is antisense; sense col = L - guide_position).
  const perAllele = candidate.allele?.per_allele ?? [];
  const alleles: AlnAllele[] = perAllele.map((pa, i) => {
    const mism: AlnMismatch[] = [];
    if (pa.discrimination) {
      const [from, to] = parseMut(pa.discrimination.mutation);
      // DNA letters from the API; the window is shown in RNA.
      const rna = (b: string) => (b === "T" ? "U" : b);
      const col = Math.min(Math.max(L - pa.discrimination.guide_position, 0), L - 1);
      mism.push({ col, from: rna(from), to: rna(to), disc: pa.discrimination });
    }
    return { name: pa.allele, diff: pa.diff_count, weak: i === 0, mism };
  });
  return { target, guide, alleles, complete: false };
}
