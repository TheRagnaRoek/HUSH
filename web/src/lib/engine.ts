/**
 * Whether the off-target check can run on this engine. It needs bowtie on the
 * PATH and the transcriptome reference; /health reports both (older engines
 * report neither, which is treated as unknown rather than as a problem).
 */

import type { HealthResponse } from "../api/types";

export type OfftargetStatus = "ready" | "no-bowtie" | "no-reference" | "unknown";

export function offtargetStatus(engine: HealthResponse | null): OfftargetStatus {
  if (!engine || engine.bowtie_available === undefined) return "unknown";
  if (!engine.bowtie_available) return "no-bowtie";
  if (engine.offtarget_reference === false) return "no-reference";
  return "ready";
}

export const OFFTARGET_SHORT: Record<"no-bowtie" | "no-reference", string> = {
  "no-bowtie": "Bowtie not found - off-target screening unavailable.",
  "no-reference": "No off-target reference built - off-target check is skipped.",
};

export const OFFTARGET_LONG: Record<"no-bowtie" | "no-reference", string> = {
  "no-bowtie":
    "Bowtie is not installed or not on the PATH, so the off-target check cannot run and is skipped. Install bowtie1 (the HUSH container includes it), then build the reference once with 'hush reference build'.",
  "no-reference":
    "No transcriptome reference has been built, so the off-target check is skipped. Build it once with 'hush reference build'.",
};

/** The reason alone, for a sentence that already says the check was skipped. */
export const OFFTARGET_REASON: Record<"no-bowtie" | "no-reference", string> = {
  "no-bowtie": "Bowtie is not installed or not on the PATH. Install bowtie1, then build the reference with 'hush reference build'.",
  "no-reference": "No transcriptome reference has been built. Build it once with 'hush reference build'.",
};
