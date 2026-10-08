/**
 * The stem-loop loop registry, mirroring shrna/loops.py. Jensen is the default:
 * library-selected for robust knockdown with the least stem-sequence dependence,
 * the safest choice for an automated converter.
 */
export interface LoopInfo {
  name: string;
  label: string;
  seq: string;
  citation: string;
}

export const LOOPS: LoopInfo[] = [
  {
    name: "jensen",
    label: "Jensen (UGUGCUU)",
    seq: "UGUGCUU",
    citation: "Jensen et al. 2012 - least stem-dependent; the safe default.",
  },
  {
    name: "psuper",
    label: "pSUPER (UUCAAGAGA)",
    seq: "UUCAAGAGA",
    citation: "Brummelkamp et al. 2002, Science 296:550.",
  },
  {
    name: "trc",
    label: "TRC / pLKO.1 (CUCGAG)",
    seq: "CUCGAG",
    citation: "Moffat et al. 2006, Cell 124:1283 (pLKO.1 / TRC).",
  },
];

export const DEFAULT_LOOP = "jensen";
export const DEFAULT_TERMINATOR_LENGTH = 6;
export const MIN_TERMINATOR_LENGTH = 4;

export function loopInfo(name: string): LoopInfo | undefined {
  return LOOPS.find((l) => l.name === name);
}

/** The Pol III poly-U terminator shown in the schematic for a given length. */
export function terminatorString(length: number | null | undefined): string {
  return length && length > 0 ? "U".repeat(length) : "";
}
