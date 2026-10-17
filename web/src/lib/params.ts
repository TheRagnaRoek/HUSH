import type { DesignParams } from "../api/types";

/**
 * Default parameters, matching DesignParamsSchema in api/schemas.py (defaults in
 * docs/cli.md). GC bounds are fractions, as the API expects.
 */
export const DEFAULT_PARAMS: DesignParams = {
  length: 19,
  overhang: "natural",

  use_ui_tei: true,
  use_reynolds: true,
  use_amarzguioui: true,
  rule_mode: "any",

  use_seed_filter: true,
  seed_tm_threshold: 21.5,

  use_offtarget_filter: true,
  offtarget_min_mismatch: 2,

  use_homopolymer_filter: true,
  max_homopolymer_run: 3,

  use_gc_filter: false,
  gc_min: 0.3,
  gc_max: 0.52,

  use_hairpin_filter: false,
  hairpin_tm_max: 20,

  use_asymmetry_filter: false,
  asymmetry_window: 4,

  min_spare_mismatches: 1,
  sort: "auto",
};

/** What the "relaxed HLA preset" button applies: drop the seed filter (too strict
 * on GC-rich families) and add the standalone GC + hairpin gates. */
export const RELAXED_HLA_PATCH: Partial<DesignParams> = {
  use_seed_filter: false,
  use_gc_filter: true,
  use_hairpin_filter: true,
};
