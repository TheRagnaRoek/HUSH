/**
 * Gate catalogue for the criteria step. Help text is distilled from
 * docs/design_options.md. `control` selects a cosmetic control renderer in
 * GateRow; TODO: bind these controls to DesignParams so a run uses them.
 */
export type ControlKind =
  | "rules"
  | "seedtm"
  | "homo"
  | "offtarget"
  | "gc"
  | "hairpin"
  | "asym"
  | "duplex";

export interface GateDef {
  id: string;
  name: string;
  def: string;
  defaultOn: boolean;
  sub: string;
  what: string;
  opt: string;
  control: ControlKind;
}

export const GATES: GateDef[] = [
  {
    id: "rules",
    name: "Classical rules",
    def: "on",
    defaultOn: true,
    sub: "Ui-Tei / Reynolds / Amarzguioui guide-efficacy rules",
    what: "The three 2004 rule sets for guide-strand efficacy, combined by union so each is weighted equally.",
    opt: "Optimal: keep all three on 'any'. Tighten to 'all' only when you have far too many candidates.",
    control: "rules",
  },
  {
    id: "seed",
    name: "Seed-duplex off-target Tm",
    def: "on",
    defaultOn: false,
    sub: "siDirect seed-mediated off-target proxy, both strands",
    what: "The siDirect seed-duplex melting-temperature proxy for miRNA-like off-targeting.",
    opt: "Turned OFF here: on a GC-rich family like HLA it rejects everything. The 21.5 C cut-off is only valid with the fixed siDirect parameters.",
    control: "seedtm",
  },
  {
    id: "homopolymer",
    name: "Homopolymer run",
    def: "on",
    defaultOn: true,
    sub: "rejects runs like GGGG that hurt stability / cause read-through",
    what: "Rejects a core with a single-nucleotide run longer than the cap. Poly-runs hurt stability and cause Pol III read-through.",
    opt: "Optimal: keep on, cap 3.",
    control: "homo",
  },
  {
    id: "offtarget",
    name: "Transcriptome off-target",
    def: "on",
    defaultOn: true,
    sub: "bowtie search of the reference transcriptome, served server-side",
    what: "A bowtie search of the reference transcriptome; rejects a guide whose nearest off-target has fewer than N mismatches.",
    opt: "Optimal: on, min 2 (matches siDirect on the PDCD1 benchmark).",
    control: "offtarget",
  },
  {
    id: "gc",
    name: "GC-content window",
    def: "off",
    defaultOn: true,
    sub: "standalone GC window, separate from Reynolds' internal criterion",
    what: "A standalone GC window, separate from Reynolds' internal GC criterion.",
    opt: "Added here for the GC-rich target to keep duplex loading in a sane range.",
    control: "gc",
  },
  {
    id: "hairpin",
    name: "Internal hairpin",
    def: "off",
    defaultOn: true,
    sub: "rejects a sense strand that folds on itself",
    what: "Rejects a sense strand that folds into a stable internal hairpin competing with duplex formation.",
    opt: "Worth adding on structured / GC-rich targets.",
    control: "hairpin",
  },
  {
    id: "asymmetry",
    name: "Loading asymmetry",
    def: "off",
    defaultOn: false,
    sub: "favours guide-strand loading (sign only)",
    what: "Requires duplex-end thermodynamic asymmetry favouring guide-strand loading (sign only, no magnitude cut-off).",
    opt: "Optional; leave off unless you want the extra stringency.",
    control: "asym",
  },
  {
    id: "duplex",
    name: "Duplex length & overhang",
    def: "19 / natural",
    defaultOn: true,
    sub: "19-nt paired core, rule sets defined on it",
    what: "The paired core length and the 3' overhang style. Overhangs are tracked apart from the core and never reach a filter.",
    opt: "Optimal: 19-nt core, natural overhangs (the true offset Tuschl duplex).",
    control: "duplex",
  },
];
