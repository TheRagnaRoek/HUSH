/**
 * TypeScript mirror of api/schemas.py. Keep in sync with it; 'npm run gen:api'
 * generates schema.d.ts from /openapi.json to diff against.
 */

// ---- requests ---- //

export interface TargetingSchema {
  pan_allele?: boolean;
  allele_specific?: string | null;
  target_specific?: boolean;
  hit?: string[] | null;
  spare?: string[] | null;
}

export interface DesignParams {
  length?: number;
  overhang?: string;

  use_ui_tei?: boolean;
  use_reynolds?: boolean;
  use_amarzguioui?: boolean;
  rule_mode?: "any" | "all";
  reynolds_threshold?: number;
  amarzguioui_threshold?: number;

  use_seed_filter?: boolean;
  seed_tm_threshold?: number;

  use_offtarget_filter?: boolean;
  offtarget_db?: string | null;
  offtarget_min_mismatch?: number;
  offtarget_max_mismatch?: number;
  offtarget_self_accessions?: string[];
  resolve_offtarget_self?: boolean;
  offtarget_self_coverage?: number; // fraction in (0, 1]
  bowtie_path?: string;

  use_gc_filter?: boolean;
  gc_min?: number; // fraction in [0, 1]
  gc_max?: number; // fraction in [0, 1]

  use_homopolymer_filter?: boolean;
  max_homopolymer_run?: number;

  use_asymmetry_filter?: boolean;
  asymmetry_window?: number;
  asymmetry_temp_c?: number;

  use_hairpin_filter?: boolean;
  hairpin_tm_max?: number;

  thermo_params?: string;

  targeting?: TargetingSchema | null;
  min_spare_mismatches?: number;
  sort?: "auto" | "position" | "discrimination";
}

export type StemOverhang = "keep" | "trim";

export interface ConvertOptions {
  loop?: string;
  stem_overhang?: StemOverhang;
  include_terminator?: boolean;
  terminator_length?: number;
  /** Withhold constructs with a premature Pol III terminator (the default). */
  drop_premature_terminators?: boolean;
}

export interface InspectRequest {
  target_fasta: string;
  panel_fasta?: string | null;
}

export interface DesignRequest {
  target_fasta: string;
  panel_fasta?: string | null;
  params?: DesignParams;
}

export interface CreateSessionRequest {
  target_fasta: string;
  panel_fasta?: string | null;
}

// ---- responses ---- //

export interface SequenceSummary {
  id: string;
  description: string;
  length: number;
  moltype: string;
  is_nucleotide: boolean;
  gc_percent?: number | null;
}

export interface PanelSummary {
  name: string;
  count: number;
  uniform: boolean;
  members: SequenceSummary[];
}

export interface InspectResponse {
  target: SequenceSummary;
  target_ready: boolean;
  note?: string | null;
  panel?: PanelSummary | null;
}

export interface CriteriaSummary {
  rules: string[];
  rule_mode: string;
  seed_filter: boolean;
  seed_tm_threshold_c?: number | null;
  offtarget_filter: boolean;
  offtarget_min_mismatch?: number | null;
  offtarget_self_accessions?: string[];
  offtarget_self_declared?: number;
  offtarget_self_resolved?: number;
  offtarget_self_coverage?: number | null;
  gc_filter: boolean;
  gc_min?: number | null;
  gc_max?: number | null;
  hairpin_filter: boolean;
  hairpin_tm_max_c?: number | null;
  homopolymer_filter: boolean;
  max_homopolymer_run?: number | null;
  asymmetry_filter: boolean;
}

export interface TargetingSummary {
  mode: string;
  hit: string[];
  spare: string[];
  discriminating: boolean;
}

export interface DiscriminationResult {
  guide_position: number;
  mutation: string;
  mismatch_class?: string | null;
  huang_recommended?: boolean | null;
  region: string;
  sources_agree: boolean;
  strongest_position: boolean;
  /** Inside the P9-P16 window; the field name is kept from the old wire format. */
  consistent_position: boolean;
  multi_mismatch: boolean;
  /** Every difference to that allele, in ranking order. */
  all_guide_positions?: number[];
}

export interface PerAlleleDiscrimination {
  allele: string;
  diff_count: number;
  discrimination?: DiscriminationResult | null;
}

export interface AlleleResult {
  hit_max_mismatch: number;
  spare_min_mismatch?: number | null;
  spare_count: number;
  discrimination?: DiscriminationResult | null;
  per_allele: PerAlleleDiscrimination[];
}

export interface CandidateResult {
  region: string;
  start: number;
  end: number;
  guide_5to3: string;
  guide_overhang: string;
  passenger_5to3: string;
  passenger_overhang: string;
  gc_percent: number;
  rules_passed: string;
  seed_max_tm_c?: number | null;
  offtarget_min_mismatch?: number | null;
  offtarget_hits?: number[] | null;
  hairpin_tm_c?: number | null;
  homopolymer_max_run?: number | null;
  allele?: AlleleResult | null;
}

export interface DesignResponse {
  target: SequenceSummary;
  length: number;
  overhang: string;
  candidate_count: number;
  offtarget_skipped: boolean;
  criteria: CriteriaSummary;
  funnel: Record<string, number>;
  targeting?: TargetingSummary | null;
  selected: CandidateResult[];
}

export interface ConstructResult {
  region: string;
  start: number;
  end: number;
  loop_name: string;
  loop_seq: string;
  orientation: string;
  hairpin_rna_5to3: string;
  insert_top_dna_5to3: string;
  insert_bottom_dna_5to3: string;
  folds_as_designed: boolean;
  fold_dg37?: number | null;
  fold_tm_c?: number | null;
  starts_with_purine: boolean;
  premature_terminator_free: boolean;
  max_homopolymer_run: number;
  sense_overhang: string;
  guide_overhang: string;
}

export interface ConvertResponse {
  loop_name: string;
  loop_seq: string;
  stem_overhang?: StemOverhang;
  terminator_length?: number | null;
  construct_count: number;
  converted_count: number;
  withheld_count: number;
  fold_ok_count: number;
  purine_start_count: number;
  constructs: ConstructResult[];
  withheld: ConstructResult[];
}

export interface SessionResponse {
  id: string;
  target: SequenceSummary;
  panel?: PanelSummary | null;
  has_last_design: boolean;
}

export interface HealthResponse {
  status: string;
  service: string;
  version: string;
  /** Absent on an engine older than this field. */
  bowtie_available?: boolean;
  offtarget_reference?: boolean;
}
