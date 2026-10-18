/**
 * Sample data in the API response shapes, so the app runs without a backend.
 * Guides and oligos are HUSH output for the HLA-A 474/476/478 windows (A*01 and
 * A*02 run); the 4-allele spared panel is embellished. The discrimination
 * fields are derived from the window sequences (see the builder below).
 *
 * The constructs come from shrna.convert on these windows, overhangs included.
 */


import type {
  CandidateResult,
  ConstructResult,
  ConvertResponse,
  DesignResponse,
  DiscriminationResult,
  InspectResponse,
} from "../api/types";
import type { AlleleEnrichment, AlnMismatch } from "./alignment";

// The discrimination fields are derived from the window sequences with the same
// rules as msa/discrimination.py, so the sample cannot contradict itself.
const RNA_COMPLEMENT: Record<string, string> = { A: "U", U: "A", G: "C", C: "G" };
const PURINES = new Set(["A", "G"]);
const STRONGEST = [10, 13, 16];

// Huang et al. 2009 Table 1, keyed "spared>target" (DNA letters) as the API writes it.
const HUANG_POSITIONS: Record<string, number[]> = {
  "A>C": [8, 9, 10, 11, 16],
  "A>G": [10, 11, 12, 16],
  "A>T": [9, 10, 14, 15, 16],
  "C>A": [11, 12, 13, 14, 16],
  "C>G": [9, 10, 12, 13, 16],
  "C>T": [9, 10, 11, 12, 14],
  "G>A": [10, 13, 15, 16],
  "G>C": [10, 11, 12, 13],
  "G>T": [9, 10, 11, 12, 14],
  "T>A": [9, 11, 12, 15],
  "T>C": [10, 11, 12, 13],
  "T>G": [10, 11, 12, 13, 16],
};

const toDna = (b: string) => (b === "U" ? "T" : b);

/** One difference at a guide position (1-based from the guide 5' end). */
function assess(sense: string, guidePos: number, spared: string, count: number): AlnMismatch {
  const col = sense.length - guidePos;
  const target = sense[col];
  const mutation = `${toDna(spared)}>${toDna(target)}`;
  const listed = HUANG_POSITIONS[mutation];
  if (!listed || spared === target) throw new Error(`Bad sample mismatch at P${guidePos}`);
  const guideBase = RNA_COMPLEMENT[target];
  const guidePurine = PURINES.has(guideBase);
  const sparedPurine = PURINES.has(spared);
  const cls = `${guidePurine ? "purine" : "pyrimidine"}:${sparedPurine ? "purine" : "pyrimidine"}`;
  const region = guidePos === 1 ? "anchor" : guidePos <= 8 ? "seed" : "post_seed";
  const huang = listed.includes(guidePos);
  return {
    col,
    from: target,
    to: spared,
    disc: {
      guide_position: guidePos,
      mutation,
      mismatch_class: cls,
      huang_recommended: huang,
      region,
      sources_agree: huang === (region === "post_seed"),
      strongest_position: STRONGEST.includes(guidePos),
      consistent_position: guidePos >= 9 && guidePos <= 16,
      multi_mismatch: count >= 2,
    },
  };
}

type Key = boolean[];

/** The five ranking booleans of one mismatch; higher is better. */
function keyOf(d: DiscriminationResult): Key {
  return [
    d.huang_recommended === true,
    d.consistent_position,
    d.strongest_position,
    d.region === "post_seed",
    d.mismatch_class === "purine:purine",
  ];
}

function compareKeys(a: Key, b: Key): number {
  for (let i = 0; i < a.length; i++) {
    if (a[i] !== b[i]) return a[i] ? 1 : -1;
  }
  return 0;
}

interface SampleWindow {
  region: string;
  start: number;
  end: number;
  guide: string;
  guideOverhang: string;
  passenger: string;
  passengerOverhang: string;
  gc: number;
  hits: number[];
  homopolymer: number;
  /** Per spared allele: [guide position, base the spared allele carries]. */
  spared: { allele: string; at: [number, string][] }[];
}

const SAMPLE_WINDOWS: SampleWindow[] = [
  {
    region: "474-492",
    start: 474,
    end: 492,
    guide: "AGCUGCCAUGUCCGCCGCGGU",
    guideOverhang: "GU",
    passenger: "CGCGGCGGACAUGGCAGCUCA",
    passengerOverhang: "CA",
    gc: 73.7,
    hits: [0, 1, 2, 4],
    homopolymer: 2,
    spared: [
      { allele: "HLA-A*02:01", at: [[6, "A"], [14, "G"]] },
      { allele: "HLA-A*03:01", at: [[6, "A"], [11, "C"], [17, "U"]] },
      { allele: "HLA-A*24:02", at: [[10, "A"], [3, "A"]] },
    ],
  },
  {
    region: "476-494",
    start: 476,
    end: 494,
    guide: "UGAGCUGCCAUGUCCGCCGCG",
    guideOverhang: "CG",
    passenger: "CGGCGGACAUGGCAGCUCAGA",
    passengerOverhang: "GA",
    gc: 68.4,
    hits: [0, 1, 2, 6],
    homopolymer: 2,
    spared: [
      { allele: "HLA-A*02:01", at: [[13, "G"], [9, "U"], [5, "A"]] },
      { allele: "HLA-A*03:01", at: [[16, "A"], [8, "A"]] },
      { allele: "HLA-A*24:02", at: [[10, "C"], [7, "A"]] },
    ],
  },
  {
    region: "478-496",
    start: 478,
    end: 496,
    guide: "UCUGAGCUGCCAUGUCCGCCG",
    guideOverhang: "CG",
    passenger: "GCGGACAUGGCAGCUCAGACC",
    passengerOverhang: "CC",
    gc: 63.2,
    hits: [0, 1, 2, 8],
    homopolymer: 2,
    spared: [
      { allele: "HLA-A*02:01", at: [[8, "G"], [12, "C"]] },
      { allele: "HLA-A*03:01", at: [[13, "C"], [5, "C"], [3, "G"]] },
      { allele: "HLA-A*24:02", at: [[10, "A"], [16, "A"]] },
    ],
  },
];

interface BuiltWindow {
  candidate: CandidateResult;
  enrichment: AlleleEnrichment[];
  key: Key;
}

function buildWindow(w: SampleWindow): BuiltWindow {
  const sense = w.passenger.slice(0, w.passenger.length - w.passengerOverhang.length);
  // Stable sorts, as in the engine: ties keep the order they were given in.
  const alleles = w.spared
    .map((s) => {
      const found = s.at.map(([p, base]) => assess(sense, p, base, s.at.length));
      const ranked = [...found].sort((a, b) => compareKeys(keyOf(b.disc), keyOf(a.disc)));
      const disc: DiscriminationResult = {
        ...ranked[0].disc,
        all_guide_positions: ranked.map((m) => m.disc.guide_position),
      };
      return { allele: s.allele, found, disc };
    })
    .sort((a, b) => compareKeys(keyOf(a.disc), keyOf(b.disc)));

  const weakest = alleles[0];
  const minDiff = Math.min(...alleles.map((a) => a.found.length));
  const candidate: CandidateResult = {
    region: w.region,
    start: w.start,
    end: w.end,
    guide_5to3: w.guide,
    guide_overhang: w.guideOverhang,
    passenger_5to3: w.passenger,
    passenger_overhang: w.passengerOverhang,
    gc_percent: w.gc,
    rules_passed: "--A",
    offtarget_min_mismatch: 2,
    offtarget_hits: w.hits,
    homopolymer_max_run: w.homopolymer,
    allele: {
      hit_max_mismatch: 0,
      spare_min_mismatch: minDiff,
      spare_count: alleles.length,
      discrimination: { ...weakest.disc, multi_mismatch: minDiff >= 2 },
      per_allele: alleles.map((a) => ({
        allele: a.allele,
        diff_count: a.found.length,
        discrimination: a.disc,
      })),
    },
  };
  const enrichment: AlleleEnrichment[] = alleles.map((a, i) => ({
    name: a.allele,
    diff: a.found.length,
    weak: i === 0,
    mism: [...a.found].sort((x, y) => x.col - y.col),
  }));
  return { candidate, enrichment, key: keyOf(weakest.disc) };
}

// Ranked as the engine does: best key first, ties by target position.
const BUILT = SAMPLE_WINDOWS.map(buildWindow).sort(
  (a, b) => compareKeys(b.key, a.key) || a.candidate.start - b.candidate.start,
);

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
  selected: BUILT.map((b) => b.candidate),
};

/** Full per-base window enrichment (what a panel-windows API addition would
 * give). Keyed by region; demonstrates the complete alignment view. */
export const sampleEnrichment: Record<string, AlleleEnrichment[]> = Object.fromEntries(
  BUILT.map((b) => [b.candidate.region, b.enrichment]),
);

/** Constructs (real HUSH convert output for these windows). */
const CONSTRUCTS: ConstructResult[] = [
  {
    region: "474-492",
    start: 474,
    end: 492,
    loop_name: "jensen",
    loop_seq: "UGUGCUU",
    orientation: "sa",
    hairpin_rna_5to3: "CGCGGCGGACAUGGCAGCUCAUGUGCUUAGCUGCCAUGUCCGCCGCGGUUUUUUU",
    insert_top_dna_5to3: "CGCGGCGGACATGGCAGCTCATGTGCTTAGCTGCCATGTCCGCCGCGGTTTTTTT",
    insert_bottom_dna_5to3: "AAAAAAACCGCGGCGGACATGGCAGCTAAGCACATGAGCTGCCATGTCCGCCGCG",
    folds_as_designed: true,
    fold_dg37: -41.32,
    fold_tm_c: 109.9,
    starts_with_purine: false,
    premature_terminator_free: true,
    max_homopolymer_run: 2,
    sense_overhang: "CA",
    guide_overhang: "GU",
  },
  {
    region: "476-494",
    start: 476,
    end: 494,
    loop_name: "jensen",
    loop_seq: "UGUGCUU",
    orientation: "sa",
    hairpin_rna_5to3: "CGGCGGACAUGGCAGCUCAGAUGUGCUUUGAGCUGCCAUGUCCGCCGCGUUUUUU",
    insert_top_dna_5to3: "CGGCGGACATGGCAGCTCAGATGTGCTTTGAGCTGCCATGTCCGCCGCGTTTTTT",
    insert_bottom_dna_5to3: "AAAAAACGCGGCGGACATGGCAGCTCAAAGCACATCTGAGCTGCCATGTCCGCCG",
    folds_as_designed: true,
    fold_dg37: -40.33,
    fold_tm_c: 108.83,
    starts_with_purine: false,
    premature_terminator_free: true,
    max_homopolymer_run: 3,
    sense_overhang: "GA",
    guide_overhang: "CG",
  },
  {
    region: "478-496",
    start: 478,
    end: 496,
    loop_name: "jensen",
    loop_seq: "UGUGCUU",
    orientation: "sa",
    hairpin_rna_5to3: "GCGGACAUGGCAGCUCAGACCUGUGCUUUCUGAGCUGCCAUGUCCGCCGUUUUUU",
    insert_top_dna_5to3: "GCGGACATGGCAGCTCAGACCTGTGCTTTCTGAGCTGCCATGTCCGCCGTTTTTT",
    insert_bottom_dna_5to3: "AAAAAACGGCGGACATGGCAGCTCAGAAAGCACAGGTCTGAGCTGCCATGTCCGC",
    folds_as_designed: true,
    fold_dg37: -38.42,
    fold_tm_c: 105.12,
    starts_with_purine: true,
    premature_terminator_free: true,
    max_homopolymer_run: 3,
    sense_overhang: "CC",
    guide_overhang: "CG",
  },
];

export const sampleConvert: ConvertResponse = {
  loop_name: "jensen",
  loop_seq: "UGUGCUU",
  stem_overhang: "keep",
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
