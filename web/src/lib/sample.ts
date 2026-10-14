/**
 * Representative sample data in the real API response shapes, so the app opens
 * in a realistic working state and runs without a backend. Guides and construct
 * oligos are real HUSH output for the HLA-A 474/476/478 windows, from an
 * A*01+A*02 run; the spared-allele scenario is embellished to a 4-allele panel
 * to exercise the per-allele view. Clearly sample - not the user's data.
 *
 * Predates the overhang change, so the constructs here carry none; regenerate
 * from a current run if the sample should show the shipped hairpin geometry.
 */

import type {
  ConstructResult,
  ConvertResponse,
  DesignResponse,
  DiscriminationResult,
  InspectResponse,
} from "../api/types";
import type { AlleleEnrichment } from "./alignment";

type Huang = "H+" | "H-" | "H?";
type Tier = "S" | "C" | "-";

function disc(
  P: number,
  region: string,
  huang: Huang,
  tier: Tier,
  cls: string,
  mutation: string,
  multi = false,
): DiscriminationResult {
  return {
    guide_position: P,
    mutation,
    mismatch_class: cls,
    huang_recommended: huang === "H+" ? true : huang === "H-" ? false : null,
    region,
    sources_agree: !(region === "seed" && huang === "H+"),
    strongest_position: tier === "S",
    consistent_position: tier === "C",
    multi_mismatch: multi,
  };
}

export const sampleInspect: InspectResponse = {
  target: {
    id: "HLA-A*01:01",
    description: "HLA-A*01:01 mRNA (CDS)",
    length: 1101,
    moltype: "RNA",
    is_nucleotide: true,
    gc_percent: 64.8,
  },
  target_ready: true,
  note: "GC-rich target: the seed filter may reject every candidate.",
  panel: {
    name: "HLA-A",
    count: 4,
    uniform: true,
    members: [
      { id: "HLA-A*01:01", description: "", length: 1101, moltype: "RNA", is_nucleotide: true },
      { id: "HLA-A*02:01", description: "", length: 1101, moltype: "RNA", is_nucleotide: true },
      { id: "HLA-A*03:01", description: "", length: 1098, moltype: "RNA", is_nucleotide: true },
      { id: "HLA-A*24:02", description: "", length: 1101, moltype: "RNA", is_nucleotide: true },
    ],
  },
};

export const sampleDesign: DesignResponse = {
  target: sampleInspect.target,
  length: 19,
  overhang: "natural",
  candidate_count: 1081,
  offtarget_skipped: false,
  criteria: {
    rules: ["Ui-Tei", "Reynolds", "Amarzguioui"],
    rule_mode: "any",
    seed_filter: false,
    offtarget_filter: true,
    offtarget_min_mismatch: 2,
    gc_filter: true,
    gc_min: 0.3,
    gc_max: 0.52,
    hairpin_filter: true,
    hairpin_tm_max_c: 20,
    homopolymer_filter: true,
    max_homopolymer_run: 3,
    asymmetry_filter: false,
  },
  funnel: {
    enumerated: 1081,
    "classical rules": 512,
    homopolymer: 468,
    "off-target": 455,
    "allele hit": 41,
    "allele spare": 3,
  },
  targeting: {
    mode: "allele-specific",
    hit: ["HLA-A*01:01"],
    spare: ["HLA-A*02:01", "HLA-A*03:01", "HLA-A*24:02"],
    discriminating: true,
  },
  selected: [
    {
      region: "474-492",
      start: 474,
      end: 492,
      guide_5to3: "AGCUGCCAUGUCCGCCGCGGU",
      guide_overhang: "GU",
      passenger_5to3: "CGCGGCGGACAUGGCAGCUCA",
      passenger_overhang: "CA",
      gc_percent: 73.7,
      rules_passed: "--A",
      offtarget_min_mismatch: 2,
      offtarget_hits: [0, 1, 2, 4],
      homopolymer_max_run: 2,
      allele: {
        hit_max_mismatch: 0,
        spare_min_mismatch: 2,
        spare_count: 3,
        discrimination: disc(6, "seed", "H-", "-", "pu", "G>A"),
        per_allele: [
          { allele: "HLA-A*02:01", diff_count: 2, discrimination: disc(6, "seed", "H-", "-", "pu", "G>A", true) },
          { allele: "HLA-A*03:01", diff_count: 3, discrimination: disc(15, "post", "H+", "C", "up", "C>U", true) },
          { allele: "HLA-A*24:02", diff_count: 4, discrimination: disc(10, "post", "H+", "S", "pp", "G>C", true) },
        ],
      },
    },
    {
      region: "476-494",
      start: 476,
      end: 494,
      guide_5to3: "UGAGCUGCCAUGUCCGCCGCG",
      guide_overhang: "CG",
      passenger_5to3: "CGGCGGACAUGGCAGCUCAGA",
      passenger_overhang: "GA",
      gc_percent: 68.4,
      rules_passed: "--A",
      offtarget_min_mismatch: 2,
      offtarget_hits: [0, 1, 2, 6],
      homopolymer_max_run: 2,
      allele: {
        hit_max_mismatch: 0,
        spare_min_mismatch: 2,
        spare_count: 3,
        discrimination: disc(16, "post", "H+", "S", "pp", "G>A"),
        per_allele: [
          { allele: "HLA-A*03:01", diff_count: 2, discrimination: disc(16, "post", "H+", "S", "pp", "G>A") },
          { allele: "HLA-A*02:01", diff_count: 3, discrimination: disc(16, "post", "H+", "S", "pp", "G>A", true) },
          { allele: "HLA-A*24:02", diff_count: 4, discrimination: disc(12, "post", "H+", "C", "pp", "G>A", true) },
        ],
      },
    },
    {
      region: "478-496",
      start: 478,
      end: 496,
      guide_5to3: "UCUGAGCUGCCAUGUCCGCCG",
      guide_overhang: "CG",
      passenger_5to3: "GCGGACAUGGCAGCUCAGACC",
      passenger_overhang: "CC",
      gc_percent: 63.2,
      rules_passed: "--A",
      offtarget_min_mismatch: 2,
      offtarget_hits: [0, 1, 2, 8],
      homopolymer_max_run: 2,
      allele: {
        hit_max_mismatch: 0,
        spare_min_mismatch: 2,
        spare_count: 3,
        discrimination: disc(8, "post", "H?", "C", "uu", "C>U"),
        per_allele: [
          { allele: "HLA-A*02:01", diff_count: 2, discrimination: disc(8, "post", "H?", "C", "uu", "C>U") },
          { allele: "HLA-A*03:01", diff_count: 3, discrimination: disc(5, "seed", "H+", "-", "up", "C>A", true) },
          { allele: "HLA-A*24:02", diff_count: 5, discrimination: disc(20, "post", "H+", "-", "pp", "C>G", true) },
        ],
      },
    },
  ],
};

/** Full per-base window enrichment (what a panel-windows API addition would
 * give). Keyed by region; demonstrates the complete alignment view. */
export const sampleEnrichment: Record<string, AlleleEnrichment[]> = {
  "474-492": [
    { name: "HLA-A*02:01", diff: 2, weak: true, mism: [{ col: 5, from: "G", to: "A", disc: disc(6, "seed", "H-", "-", "pu", "G>A") }] },
    {
      name: "HLA-A*03:01",
      diff: 3,
      weak: false,
      mism: [
        { col: 5, from: "G", to: "A", disc: disc(6, "seed", "H-", "-", "pu", "G>A") },
        { col: 14, from: "C", to: "U", disc: disc(15, "post", "H+", "C", "up", "C>U") },
      ],
    },
    {
      name: "HLA-A*24:02",
      diff: 4,
      weak: false,
      mism: [
        { col: 5, from: "G", to: "A", disc: disc(6, "seed", "H-", "-", "pu", "G>A") },
        { col: 9, from: "G", to: "C", disc: disc(10, "post", "H+", "S", "pp", "G>C") },
      ],
    },
  ],
  "476-494": [
    { name: "HLA-A*03:01", diff: 2, weak: true, mism: [{ col: 6, from: "G", to: "A", disc: disc(16, "post", "H+", "S", "pp", "G>A") }] },
    {
      name: "HLA-A*02:01",
      diff: 3,
      weak: false,
      mism: [{ col: 6, from: "G", to: "A", disc: disc(16, "post", "H+", "S", "pp", "G>A") }],
    },
    {
      name: "HLA-A*24:02",
      diff: 4,
      weak: false,
      mism: [
        { col: 6, from: "G", to: "A", disc: disc(16, "post", "H+", "S", "pp", "G>A") },
        { col: 2, from: "G", to: "A", disc: disc(12, "post", "H+", "C", "pp", "G>A") },
      ],
    },
  ],
  "478-496": [
    { name: "HLA-A*02:01", diff: 2, weak: true, mism: [{ col: 8, from: "C", to: "U", disc: disc(8, "post", "H?", "C", "uu", "C>U") }] },
    {
      name: "HLA-A*03:01",
      diff: 3,
      weak: false,
      mism: [
        { col: 8, from: "C", to: "U", disc: disc(8, "post", "H?", "C", "uu", "C>U") },
        { col: 16, from: "C", to: "A", disc: disc(5, "seed", "H+", "-", "up", "C>A") },
      ],
    },
    {
      name: "HLA-A*24:02",
      diff: 5,
      weak: false,
      mism: [
        { col: 8, from: "C", to: "U", disc: disc(8, "post", "H?", "C", "uu", "C>U") },
        { col: 1, from: "C", to: "G", disc: disc(20, "post", "H+", "-", "pp", "C>G") },
      ],
    },
  ],
};

/** Constructs (real HUSH convert output for these windows). */
const CONSTRUCTS: ConstructResult[] = [
  {
    region: "474-492",
    start: 474,
    end: 492,
    loop_name: "jensen",
    loop_seq: "UGUGCUU",
    orientation: "sa",
    hairpin_rna_5to3: "CGCGGCGGACAUGGCAGCUUGUGCUUAGCUGCCAUGUCCGCCGCGUUUUUU",
    insert_top_dna_5to3: "CGCGGCGGACATGGCAGCTTGTGCTTAGCTGCCATGTCCGCCGCGTTTTTT",
    insert_bottom_dna_5to3: "AAAAAACGCGGCGGACATGGCAGCTAAGCACAAGCTGCCATGTCCGCCGCG",
    folds_as_designed: true,
    fold_dg37: -41.82,
    fold_tm_c: 111.0,
    starts_with_purine: false,
    premature_terminator_free: true,
    max_homopolymer_run: 2,
    sense_overhang: "",
    guide_overhang: "",
  },
  {
    region: "476-494",
    start: 476,
    end: 494,
    loop_name: "jensen",
    loop_seq: "UGUGCUU",
    orientation: "sa",
    hairpin_rna_5to3: "CGGCGGACAUGGCAGCUCAUGUGCUUUGAGCUGCCAUGUCCGCCGUUUUUU",
    insert_top_dna_5to3: "CGGCGGACATGGCAGCTCATGTGCTTTGAGCTGCCATGTCCGCCGTTTTTT",
    insert_bottom_dna_5to3: "AAAAAACGGCGGACATGGCAGCTCAAAGCACATGAGCTGCCATGTCCGCCG",
    folds_as_designed: true,
    fold_dg37: -40.49,
    fold_tm_c: 109.17,
    starts_with_purine: false,
    premature_terminator_free: true,
    max_homopolymer_run: 3,
    sense_overhang: "",
    guide_overhang: "",
  },
  {
    region: "478-496",
    start: 478,
    end: 496,
    loop_name: "jensen",
    loop_seq: "UGUGCUU",
    orientation: "sa",
    hairpin_rna_5to3: "GCGGACAUGGCAGCUCAGAUGUGCUUUCUGAGCUGCCAUGUCCGCUUUUUU",
    insert_top_dna_5to3: "GCGGACATGGCAGCTCAGATGTGCTTTCTGAGCTGCCATGTCCGCTTTTTT",
    insert_bottom_dna_5to3: "AAAAAAGCGGACATGGCAGCTCAGAAAGCACATCTGAGCTGCCATGTCCGC",
    folds_as_designed: true,
    fold_dg37: -39.31,
    fold_tm_c: 107.06,
    starts_with_purine: true,
    premature_terminator_free: true,
    max_homopolymer_run: 3,
    sense_overhang: "",
    guide_overhang: "",
  },
];

export const sampleConvert: ConvertResponse = {
  loop_name: "jensen",
  loop_seq: "UGUGCUU",
  terminator_length: 6,
  construct_count: CONSTRUCTS.length,
  fold_ok_count: CONSTRUCTS.filter((c) => c.folds_as_designed).length,
  purine_start_count: CONSTRUCTS.filter((c) => c.starts_with_purine).length,
  converted_count: CONSTRUCTS.length,
  withheld_count: 0,
  constructs: CONSTRUCTS,
  withheld: [],
};

export function constructForRegion(region: string): ConstructResult | undefined {
  return CONSTRUCTS.find((c) => c.region === region);
}
