/**
 * The in-app reference: a brief overview of each design rule and why it exists,
 * condensed from docs/design-rules.md and docs/shrna-conversion.md (the full
 * reference). Section ids are the targets of the "Read more" links in the
 * (i) panels and hover tips.
 */

import type { SourceKey } from "./sources";

export type RefId =
  | "overview"
  | "duplex"
  | "rules"
  | "seed"
  | "offtarget"
  | "homopolymer"
  | "gc"
  | "hairpin"
  | "asymmetry"
  | "thermo"
  | "allele"
  | "discrimination"
  | "shrna"
  | "loops"
  | "qc";

export type RefGroup = "Selection" | "Allele targeting" | "Hairpin construct";

export interface RefSection {
  id: RefId;
  title: string;
  group: RefGroup;
  summary: string;
  points: string[];
  table?: { head: string[]; rows: string[][] };
  legend?: [string, string][];
  defaults?: string;
  cites: SourceKey[];
}

export const GROUPS: RefGroup[] = ["Selection", "Allele targeting", "Hairpin construct"];

export const SECTIONS: RefSection[] = [
  {
    id: "overview",
    title: "How selection works",
    group: "Selection",
    summary:
      "Selection is a filter, not a score. Every criterion is an independent, optional gate, and a candidate is selected only if it clears every enabled one.",
    points: [
      "Nothing is combined into a single number. A composite score hides which property disqualified a candidate, and the three classical rule sets were each published with their own thresholds, which a weighted blend would invalidate.",
      "The gates test the paired core of the duplex. Overhangs are stored separately and never scored.",
      "With a panel, the allele gate is combined (AND) with every other gate, and the survivors are ranked by how well they discriminate between alleles.",
    ],
    cites: [],
  },
  {
    id: "duplex",
    title: "The duplex",
    group: "Selection",
    summary:
      "The canonical duplex has a 19-nt paired core and symmetric 2-nt 3′ overhangs, so the synthesised strands are 21 nt.",
    points: [
      "The length setting is the length of the core, not of the strand. All three rule sets are defined on the core, so scoring follows the source papers whatever the overhang policy.",
      "Natural (default): the passenger overhang is the 2 nt 3′ of the core in the target, and the guide overhang is the reverse complement of the 2 nt 5′ of the core, which gives the offset duplex of Elbashir et al. Where a flank runs past the end of the transcript, dTdT is used for that overhang.",
      "dTdT is 2′-deoxythymidine, which is tolerated at the 3′ end and is the common synthetic default. UU is ribo-uridine.",
      "The two overhangs come from opposite flanks, so the guide overhang is not the reverse complement of the passenger overhang.",
    ],
    cites: ["elbashir2001"],
  },
  {
    id: "rules",
    title: "Classical rules",
    group: "Selection",
    summary:
      "Three published efficacy rule sets, reimplemented as published. By default a candidate passes if any enabled rule set accepts it; the stricter mode requires all of them.",
    points: [
      "Ui-Tei: four pass/fail rules. A or U at the guide 5′ end, G or C at the sense 5′ end, at least 4 A/U among the first 7 nt of the guide, and no G/C run longer than 9 nt. The first three rules express one principle: RISC keeps the strand whose 5′ end lies at the less stable end of the duplex.",
      "Reynolds: a score out of 10 on the sense strand, from GC content of 30-52%, no internal repeat, and base preferences at positions 3, 10, 13 and 15-19, with penalties for G at 13 and G or C at 19. The range is -2 to +10 and the default cut-off is 6.",
      "Amarzguioui: an additive score on the sense strand. Motif points (G or C at 1, A at 6 and A or U at 19 count +1; U at 1 and G at 19 count -1) plus DT3, the A/U count in the last 3 nt minus the first 3 nt. A positive DT3 favours loading of the guide. The default cut-off is 3, inclusive.",
      "Each rule set keeps the thresholds it was published with. The constants are fixed for a 19-mer and are not rescaled for longer cores.",
    ],
    defaults: "On by default, combined with 'any'.",
    cites: ["uitei2004", "reynolds2004", "amarzguioui2004", "khvorova2003", "schwarz2003"],
  },
  {
    id: "seed",
    title: "Seed-duplex off-target Tm",
    group: "Selection",
    summary:
      "Rejects a candidate unless the seed duplex of both strands has a melting temperature below the threshold. The seed is positions 2-8 of a strand.",
    points: [
      "Seed-mediated off-target silencing correlates with the thermodynamic stability of the seed duplex: a weakly pairing seed causes few off-target effects even where matching sites exist.",
      "Both strands are filtered, because the passenger is loaded often enough to cause its own seed-dependent off-target silencing. Position 1 lies in the Argonaute MID pocket and does not pair with the target.",
      "The filter is strict. On a very GC-rich target such as HLA it can reject every candidate. That is expected, and the filter can be switched off.",
      "The 21.5 C cut-off is empirical and only valid with the fixed parameters of Naito et al. (Freier stacking values, 100 uM strands, 100 mM Na+). Only the threshold can be changed.",
      "The seed duplex does not form in solution: Argonaute changes the binding properties of its guide. HUSH keeps the filter because it is the standard published one, and treats it as a calibrated empirical proxy, not as a mechanism.",
    ],
    defaults: "On by default, threshold 21.5 C.",
    cites: ["naito2009", "uitei2008", "salomon2015", "wee2012", "freier1986"],
  },
  {
    id: "offtarget",
    title: "Transcriptome off-target",
    group: "Selection",
    summary:
      "Checks whether a guide has a near match elsewhere in the human transcriptome. A candidate is rejected when its nearest off-target has fewer than the minimum number of mismatches.",
    points: [
      "Both strands are searched, because either can be loaded into RISC. Each candidate contributes two 19-nt cores; overhangs are not searched. Only the mismatch count is used, not the positions, because the seed filter already covers the seed region.",
      "Only the forward strand of each transcript is searched, because transcripts are single-stranded sense mRNA.",
      "The design's own transcripts are excluded from the hits, so the gate measures off-target risk and not homology to the target. They are the union of declared accessions, the accession in the target's FASTA header, and transcripts that match the target and panel by sequence.",
      "The search uses bowtie1 in its mismatch-bounded mode over a prebuilt index, because an end-to-end search bounded by a mismatch count is the question this gate asks. Build the reference once with 'hush reference build'.",
      "Without bowtie or the reference the check is skipped. Guides from such a run have not been screened, and an empty hit count must not be read as 'no off-targets found'.",
    ],
    defaults: "On by default, minimum 2 mismatches, which matches Naito et al. 2009 on the PDCD1 benchmark.",
    cites: ["naito2009"],
  },
  {
    id: "homopolymer",
    title: "Homopolymer runs",
    group: "Selection",
    summary:
      "A core with a single-nucleotide run longer than the maximum is rejected.",
    points: [
      "Homopolymer runs reduce siRNA stability and cause polymerase read-through.",
      "The gate does not depend on the base and checks the sense core only.",
      "U runs in an expressed hairpin are checked separately, after the hairpin is built (see construct QC).",
    ],
    defaults: "On by default, maximum run 3.",
    cites: [],
  },
  {
    id: "gc",
    title: "GC window",
    group: "Selection",
    summary: "An optional, tunable limit on the GC content of the core.",
    points: [
      "It stands alongside the GC criterion inside the Reynolds score, which keeps its published form and cannot be changed.",
      "Worth adding on a GC-rich target, to keep duplex loading in a sane range.",
    ],
    defaults: "Off by default, window 30-52%.",
    cites: ["reynolds2004"],
  },
  {
    id: "hairpin",
    title: "Internal hairpin",
    group: "Selection",
    summary:
      "Rejects a sense strand that folds into a stable self-structure, which would compete with duplex formation.",
    points: [
      "It is the tunable counterpart to the internal-repeat criterion of Reynolds.",
      "The structure is the most stable simple hairpin: a perfectly paired stem of at least 3 bp closing a loop of at least 3 nt. The cut-off is compared with the melting temperature of the fold with the lowest free energy, and a strand that folds into nothing passes.",
      "It is an approximation, not a reproduction of the folding method of the Reynolds paper, which the paper does not specify.",
    ],
    defaults: "Off by default, melting temperature at most 20 C.",
    cites: ["reynolds2004", "xia1998", "groebe1988"],
  },
  {
    id: "asymmetry",
    title: "Loading asymmetry",
    group: "Selection",
    summary:
      "Requires the less stable end of the duplex to favour loading of the guide.",
    points: [
      "RISC keeps the strand whose 5′ end lies at the less stable end of the duplex.",
      "Only the sign of the end-stability difference is used. There is no magnitude threshold, because converting kcal/mol into a cut-off needs calibration data that HUSH does not have.",
    ],
    defaults: "Off by default.",
    cites: ["khvorova2003", "schwarz2003"],
  },
  {
    id: "thermo",
    title: "Thermodynamic model",
    group: "Selection",
    summary:
      "The asymmetry and hairpin gates use a nearest-neighbour model of RNA stacking energies.",
    points: [
      "Xia et al. is the default parameter set, and Freier et al. can be chosen for the asymmetry and hairpin read-outs. A published threshold is only valid with the parameter set it was calibrated on, so the set is always named explicitly; the seed filter always uses Freier.",
      "Free energy is reported at 37 C. Mismatches, G-U pairs and dangling ends are not modelled, and there is no salt correction.",
      "Hairpin loops use the loop destabilisation values of Groebe & Uhlenbeck.",
    ],
    cites: ["xia1998", "freier1986", "groebe1988"],
  },
  {
    id: "allele",
    title: "Allele targeting",
    group: "Allele targeting",
    summary:
      "With a panel, each window is checked for whether it hits the alleles to silence and discriminates against the alleles to spare. The result is combined (AND) with every other gate.",
    points: [
      "Four requests: pan-allele (conserved sites across the family), allele-specific (hit one member and spare the rest), target-specific (spare the whole panel, for a target loaded apart from its reference panel) and custom hit/spare sets.",
      "Coverage is required. A window in which a targeted or spared allele is not covered cannot satisfy the request: missing data is never treated as evidence.",
      "A window is scored by its weakest spared allele, because an siRNA is only as specific as the allele it distinguishes least well.",
      "The gate counts differences, at least a minimum per spared allele. It does not weigh them by type or position; placement only enters through the ranking.",
      "Panel members are aligned pairwise to the target and projected onto target coordinates. The method suits near-identical members, such as the HLA-A alleles that are 95-98% identical, and loses accuracy as members diverge.",
    ],
    cites: [],
  },
  {
    id: "discrimination",
    title: "Discrimination ranking",
    group: "Allele targeting",
    summary:
      "When a request spares an allele, the windows that pass the allele gate are ordered by how reliably their difference from the spared allele is expected to discriminate. The ranking never removes a window.",
    table: {
      head: ["#", "Criterion", "True when", "Source"],
      rows: [
        ["1", "Listed for this base change", "The table of Huang et al. lists the substitution at that position", "Huang 2009"],
        ["2", "Inside the low-tolerance window", "Guide position P9-P16", "Du 2005, Schwarz, Huang"],
        ["3", "At a replicated position", "P10, P13 or P16", "Schwarz 2006"],
        ["4", "3′ of the seed", "Position 9 or later", "Schwarz"],
        ["5", "Purine:purine mismatch", "Guide and spared base are both purines", "Schwarz, Du"],
      ],
    },
    points: [
      "The key is a lexicographic tuple of five booleans, not a weighted score. The criteria are checked in order; the first one that separates two windows decides, and nothing is weighted, so nothing needs calibration.",
      "A window is ranked by its best-placed difference from its weakest spared allele. Ties are broken by target position.",
      "A seed mismatch (P2-P8) blocks binding and can be overcome by a high siRNA concentration. A mismatch 3′ of the seed blocks cleavage and cannot.",
      "The number of differences is not ranked, because no source says whether several poorly placed differences outweigh one well-placed one. The count is reported next to the rank so you can weigh it.",
      "A high rank does not guarantee that the spared allele escapes silencing: Du et al. found that 73% of singly mismatched targets were still significantly silenced. A difference is necessary for discrimination but not sufficient.",
      "The order of the five criteria is only partly supported by evidence. Listing first is well founded; the window before single positions, and position before class, are judgements, and the class criterion is the one most open to reordering.",
    ],
    legend: [
      ["Position", "Guide position, 1-based from the 5′ end of the guide (P1, P2 ...)."],
      ["Region", "Seed (P2-P8), post (P9 or later, 3′ of the seed) or anchor (P1, which does not pair)."],
      ["Huang", "H+ if the table of Huang et al. lists this base change at this position, H- if not, H? for ambiguous bases."],
      ["Tier", "S for P10, P13 or P16; W for the rest of P9-P16; a dash outside P9-P16."],
      ["Class", "The guide base against the spared base: pp (purine:purine) discriminates best, pu and up are milder, uu (pyrimidine:pyrimidine) is the mildest."],
    ],
    cites: ["huang2009", "schwarz2006", "du2005"],
  },
  {
    id: "shrna",
    title: "shRNA conversion",
    group: "Hairpin construct",
    summary:
      "Converts a candidate into a cloneable stem-loop hairpin, as top and bottom DNA oligos in sense-loop-antisense order with a Pol III poly-T terminator. The oligos carry no vector-specific restriction tails; add the sticky ends for your vector in the lab.",
    points: [
      "The design follows the canonical pSUPER system: a Pol III promoter, an inverted repeat, a loop, the complementary repeat and a poly-T terminator.",
      "Both strands are encoded in full, each as paired core plus 3′ overhang, so the hairpin expresses exactly the guide the design stage validated. The sense overhang lies between stem and loop and is unpaired, so it counts as part of the loop. The guide overhang is a single-stranded 3′ end and joins the terminator.",
      "Sense overhang, keep or trim: keep (default) encodes the full 21-nt sense strand; trim drops the sense overhang so the loop joins the core directly. The guide overhang is kept either way. Published designs disagree: Dicer cuts more precisely when the cleavage site is 2 nt from a loop or bulge, so the choice can change the guide that is produced. HUSH exposes the choice and reports the same checks for both.",
      "dTdT overhangs are encoded as ribo-U, because an expressed hairpin cannot contain 2′-deoxythymidine. The ordered DNA oligo contains T in both cases.",
      "Only Dicer-processed shRNA is supported.",
    ],
    cites: ["brummelkamp2002", "gu2012"],
  },
  {
    id: "loops",
    title: "Loops",
    group: "Hairpin construct",
    summary: "Three loops are available, each from a published source.",
    table: {
      head: ["Loop", "Sequence", "Source"],
      rows: [
        ["Jensen (default)", "UGUGCUU", "Jensen 2012"],
        ["pSUPER", "UUCAAGAGA", "Brummelkamp 2002"],
        ["TRC / pLKO.1", "CUCGAG", "Moffat 2006"],
      ],
    },
    points: [
      "Jensen is the default because it was selected from randomised retroviral libraries for robust knockdown with the least dependence on the stem sequence. That suits a tool that does not know the stem in advance.",
      "The TRC loop contains an XhoI site (CTCGAG as DNA).",
      "Loop length affects knockdown: 8-nt loops perform best overall, and a 4-nt loop gives less than 50% knockdown.",
    ],
    cites: ["jensen2012", "brummelkamp2002", "moffat2006", "wang2013"],
  },
  {
    id: "qc",
    title: "Construct checks",
    group: "Hairpin construct",
    summary:
      "Each construct is folded and checked; it is not inserted into a template without looking at its structure.",
    points: [
      "Folds as designed: the body is folded and the most stable structure is compared with the intended one, a stem equal to the core and a loop equal to the sense overhang plus the loop sequence. No means a competing structure is more stable. Known interaction: the Jensen loop ends in UU, so an AA sense overhang pairs with it and extends the stem to 21 bp. HUSH reports this and does not correct it. Stems of 19 bp or shorter behave differently and can be very potent.",
      "Pol III terminator: a run of T (U in the transcript) after the guide overhang, 6 by default and at least 4. T4 is the minimal Pol III termination signal, and full efficiency needs six or more.",
      "No premature terminator: Pol III terminates within a run of four or more T, so a UUUU anywhere in the hairpin body truncates the transcript. Constructs with such a run are withheld by default, and the count is always reported. The planned poly-U tail is not counted.",
      "Pol III start (5′-A/G): a human U6 or 7SK promoter initiates efficiently, with a precise start site, on A or G at +1, so a C or U start is flagged. It is a flag, not a gate, because the right response depends on the promoter you use.",
      "Max run: the longest run of a single nucleotide in the hairpin body, not counting the poly-U tail.",
    ],
    cites: ["bogenhagen1981", "gao2018", "gao2017", "ge2010"],
  },
];
