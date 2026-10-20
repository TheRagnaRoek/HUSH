/**
 * TSV export in the column layout of the CLI (app/cli.py), so a file saved from
 * the browser reads back with 'hush convert'.
 */

import type { CandidateResult, ConstructResult, DesignResponse } from "../api/types";

const DESIGN_BASE = [
  "region",
  "start",
  "end",
  "guide_5to3",
  "guide_overhang",
  "passenger_5to3",
  "passenger_overhang",
  "gc_percent",
  "rules_passed",
  "seed_max_tm_c",
  "hairpin_tm_c",
  "homopolymer_max_run",
  "offtarget_checked",
  "offtarget_min_mismatch",
  "offtarget_hits_0",
  "offtarget_hits_1",
  "offtarget_hits_2",
  "offtarget_hits_3",
];

const DESIGN_ALLELE = [
  "allele_hit_max_mismatch",
  "allele_spare_min_mismatch",
  "allele_spare_count",
  "disc_guide_position",
  "disc_mutation",
  "disc_mismatch_class",
  "disc_huang_recommended",
  "disc_region",
  "disc_sources_agree",
  "disc_strongest_position",
  "disc_in_window",
  "disc_multi_mismatch",
  "disc_all_positions",
];

const CONVERT_COLUMNS = [
  "region",
  "start",
  "end",
  "loop_name",
  "loop_seq",
  "orientation",
  "hairpin_rna_5to3",
  "insert_top_dna_5to3",
  "insert_bottom_dna_5to3",
  "folds_as_designed",
  "fold_dg37",
  "fold_tm_c",
  "starts_with_purine",
  "premature_terminator_free",
  "max_homopolymer_run",
  "sense_overhang",
  "guide_overhang",
];

// Python-style booleans, as the CLI writes them.
const pyBool = (b: boolean) => (b ? "True" : "False");
const fixed = (n: number | null | undefined, digits: number) =>
  n == null ? "" : n.toFixed(digits);
const text = (v: number | string | null | undefined) => (v == null ? "" : String(v));

function designRow(c: CandidateResult, offtargetChecked: boolean, withAllele: boolean): string[] {
  const hits = offtargetChecked ? (c.offtarget_hits ?? []) : [];
  const row = [
    c.region,
    text(c.start),
    text(c.end),
    c.guide_5to3,
    c.guide_overhang,
    c.passenger_5to3,
    c.passenger_overhang,
    c.gc_percent.toFixed(1),
    c.rules_passed,
    fixed(c.seed_max_tm_c, 2),
    fixed(c.hairpin_tm_c, 2),
    text(c.homopolymer_max_run),
    pyBool(offtargetChecked),
    offtargetChecked ? text(c.offtarget_min_mismatch) : "",
    ...[0, 1, 2, 3].map((i) => text(hits[i])),
  ];
  if (!withAllele) return row;

  const a = c.allele;
  if (!a) return [...row, ...DESIGN_ALLELE.map(() => "")];
  const d = a.discrimination;
  return [
    ...row,
    text(a.hit_max_mismatch),
    text(a.spare_min_mismatch),
    text(a.spare_count),
    text(d?.guide_position),
    text(d?.mutation),
    text(d?.mismatch_class),
    d?.huang_recommended == null ? "" : pyBool(d.huang_recommended),
    text(d?.region),
    d ? pyBool(d.sources_agree) : "",
    d ? pyBool(d.strongest_position) : "",
    d ? pyBool(d.consistent_position) : "",
    d ? pyBool(d.multi_mismatch) : "",
    (d?.all_guide_positions ?? []).map((p) => `P${p}`).join(","),
  ];
}

/** The selected candidates as a 'hush design' TSV. */
export function designTsv(design: DesignResponse): string {
  const withAllele = design.selected.some((c) => c.allele != null);
  const header = withAllele ? [...DESIGN_BASE, ...DESIGN_ALLELE] : DESIGN_BASE;
  const checked = !design.offtarget_skipped;
  const rows = design.selected.map((c) => designRow(c, checked, withAllele));
  return [header, ...rows].map((r) => r.join("\t")).join("\n") + "\n";
}

/** Constructs as a 'hush convert' TSV. */
export function constructTsv(constructs: ConstructResult[]): string {
  const rows = constructs.map((c) => [
    c.region,
    text(c.start),
    text(c.end),
    c.loop_name,
    c.loop_seq,
    c.orientation,
    c.hairpin_rna_5to3,
    c.insert_top_dna_5to3,
    c.insert_bottom_dna_5to3,
    pyBool(c.folds_as_designed),
    fixed(c.fold_dg37, 2),
    fixed(c.fold_tm_c, 2),
    pyBool(c.starts_with_purine),
    pyBool(c.premature_terminator_free),
    text(c.max_homopolymer_run),
    c.sense_overhang,
    c.guide_overhang,
  ]);
  return [CONVERT_COLUMNS, ...rows].map((r) => r.join("\t")).join("\n") + "\n";
}

/** Both formats hold the same tab-separated table; .txt just opens anywhere. */
export type ExportFormat = "tsv" | "txt";

const MIME: Record<ExportFormat, string> = {
  tsv: "text/tab-separated-values",
  txt: "text/plain",
};

/** Sample data is exported too, but the file name says so. */
const sampleTag = (offline: boolean) => (offline ? "_sample" : "");

/** A target id made safe for a file name ("HLA-A*01:01" becomes "HLA-A_01_01"). */
function fileSafe(name: string): string {
  const safe = name
    .trim()
    .replace(/[^A-Za-z0-9._-]+/g, "_")
    .replace(/^_+|_+$/g, "")
    .slice(0, 60);
  return safe || "target";
}

export function exportDesign(
  design: DesignResponse,
  offline: boolean,
  format: ExportFormat = "tsv",
): void {
  downloadText(
    `hush_design_${fileSafe(design.target.id)}${sampleTag(offline)}.${format}`,
    designTsv(design),
    format,
  );
}

export function exportConstruct(
  construct: ConstructResult,
  targetId: string,
  offline: boolean,
  format: ExportFormat = "tsv",
): void {
  downloadText(
    `hush_construct_${fileSafe(targetId)}_${construct.region}${sampleTag(offline)}.${format}`,
    constructTsv([construct]),
    format,
  );
}

/** Every kept construct of a conversion, as the CLI's 'convert -o' file. */
export function exportConstructs(
  constructs: ConstructResult[],
  targetId: string,
  offline: boolean,
  format: ExportFormat = "tsv",
): void {
  downloadText(
    `hush_constructs_${fileSafe(targetId)}${sampleTag(offline)}.${format}`,
    constructTsv(constructs),
    format,
  );
}

/** Save text as a file through the browser's download. */
export function downloadText(
  filename: string,
  content: string,
  format: ExportFormat = "tsv",
): void {
  const url = URL.createObjectURL(new Blob([content], { type: MIME[format] }));
  const link = document.createElement("a");
  link.href = url;
  link.download = filename;
  document.body.appendChild(link);
  link.click();
  link.remove();
  URL.revokeObjectURL(url);
}
