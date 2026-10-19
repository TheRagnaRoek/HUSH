/**
 * Gate catalogue for the criteria step. Each gate maps to fields on DesignParams
 * (api/types.ts); GateRow renders a toggle + controls bound to those fields.
 * Help text is distilled from docs/cli.md.
 */
import type { DesignParams } from "../../api/types";
import type { RefId } from "../../reference/content";
import type { SourceKey } from "../../reference/sources";

export type ControlKind =
  | "rules"
  | "seedtm"
  | "homo"
  | "offtarget"
  | "gc"
  | "hairpin"
  | "asym"
  | "duplex";

/** DesignParams boolean fields a simple gate toggles. */
export type BoolKey =
  | "use_seed_filter"
  | "use_homopolymer_filter"
  | "use_offtarget_filter"
  | "use_gc_filter"
  | "use_hairpin_filter"
  | "use_asymmetry_filter";

export interface GateDef {
  id: string;
  name: string;
  sub: string;
  what: string;
  opt: string;
  control: ControlKind;
  /** Sources named in the info panel, and the reference section it links to. */
  cites: SourceKey[];
  ref: RefId;
  /** "simple": enabledKey is the toggled boolean. "rules": toggles the three
   * rule-set flags together. "none": no toggle (settings always apply). */
  toggle: "simple" | "rules" | "none";
  enabledKey?: BoolKey;
}

export const GATES: GateDef[] = [
  {
    id: "rules",
    cites: ["uitei2004", "reynolds2004", "amarzguioui2004"],
    ref: "rules",
    name: "Classical rules",
    sub: "Ui-Tei / Reynolds / Amarzguioui guide-efficacy rules",
    what: "The three 2004 rule sets for guide-strand efficacy, combined by union so each is weighted equally.",
    opt: "Optimal: keep all three on 'any'. Tighten to 'all' only when you have far too many candidates.",
    control: "rules",
    toggle: "rules",
  },
  {
    id: "seed",
    cites: ["naito2009", "uitei2008"],
    ref: "seed",
    name: "Seed-duplex off-target Tm",
    sub: "Seed-mediated off-target proxy (Naito et al. 2009), both strands",
    what: "The seed-duplex melting-temperature proxy of Naito et al. 2009 for miRNA-like off-targeting.",
    opt: "Turn OFF on a GC-rich family like HLA, where it rejects everything. The 21.5 C cut-off is only valid with the fixed parameters of Naito et al. 2009.",
    control: "seedtm",
    toggle: "simple",
    enabledKey: "use_seed_filter",
  },
  {
    id: "homopolymer",
    cites: [],
    ref: "homopolymer",
    name: "Homopolymer run",
    sub: "Rejects runs like GGGG that hurt stability / cause read-through",
    what: "Rejects a core with a single-nucleotide run longer than the cap. Poly-runs hurt stability and cause Pol III read-through.",
    opt: "Optimal: keep on, cap 3.",
    control: "homo",
    toggle: "simple",
    enabledKey: "use_homopolymer_filter",
  },
  {
    id: "offtarget",
    cites: ["naito2009"],
    ref: "offtarget",
    name: "Transcriptome off-target",
    sub: "Bowtie search of the reference transcriptome built by 'hush reference build'",
    what: "A bowtie search of the reference transcriptome; rejects a guide whose nearest off-target has fewer than N mismatches.",
    opt: "Optimal: on, min 2 (matches Naito et al. 2009 on the PDCD1 benchmark).",
    control: "offtarget",
    toggle: "simple",
    enabledKey: "use_offtarget_filter",
  },
  {
    id: "gc",
    cites: ["reynolds2004"],
    ref: "gc",
    name: "GC-content window",
    sub: "Standalone GC window, separate from Reynolds' internal criterion",
    what: "A standalone GC window, separate from Reynolds' internal GC criterion.",
    opt: "Add on a GC-rich target to keep duplex loading in a sane range.",
    control: "gc",
    toggle: "simple",
    enabledKey: "use_gc_filter",
  },
  {
    id: "hairpin",
    cites: ["reynolds2004", "xia1998"],
    ref: "hairpin",
    name: "Internal hairpin",
    sub: "Rejects a sense strand that folds on itself",
    what: "Rejects a sense strand that folds into a stable internal hairpin competing with duplex formation.",
    opt: "Worth adding on structured / GC-rich targets.",
    control: "hairpin",
    toggle: "simple",
    enabledKey: "use_hairpin_filter",
  },
  {
    id: "asymmetry",
    cites: ["khvorova2003", "schwarz2003"],
    ref: "asymmetry",
    name: "Loading asymmetry",
    sub: "Favours guide-strand loading (sign only)",
    what: "Requires duplex-end thermodynamic asymmetry favouring guide-strand loading (sign only, no magnitude cut-off).",
    opt: "Optional; leave off unless you want the extra stringency.",
    control: "asym",
    toggle: "simple",
    enabledKey: "use_asymmetry_filter",
  },
  {
    id: "duplex",
    cites: ["elbashir2001"],
    ref: "duplex",
    name: "Duplex length & overhang",
    sub: "19-nt paired core, rule sets defined on it",
    what: "The paired core length and the 3' overhang style. Overhangs are tracked apart from the core and never reach a filter.",
    opt: "Optimal: 19-nt core, natural overhangs (the true offset Tuschl duplex).",
    control: "duplex",
    toggle: "none",
  },
];

/** True when any classical rule set is enabled. */
export function rulesEnabled(p: DesignParams): boolean {
  return Boolean(p.use_ui_tei || p.use_reynolds || p.use_amarzguioui);
}
