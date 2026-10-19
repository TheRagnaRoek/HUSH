import { createContext } from "react";
import type {
  CandidateResult,
  ConstructResult,
  ConvertOptions,
  ConvertResponse,
  DesignParams,
  DesignResponse,
  HealthResponse,
  InspectResponse,
} from "../api/types";
import type { AlleleEnrichment } from "../lib/alignment";

export type Step = "input" | "targeting" | "criteria" | "candidates" | "construct";

export const STEPS: { id: Step; label: string }[] = [
  { id: "input", label: "Target & panel" },
  { id: "targeting", label: "Targeting" },
  { id: "criteria", label: "Criteria" },
  { id: "candidates", label: "Candidates" },
  { id: "construct", label: "shRNA construct" },
];

export interface SessionValue {
  step: Step;
  setStep: (s: Step) => void;
  expert: boolean;
  setExpert: (v: boolean) => void;

  targetFasta: string;
  setTargetFasta: (v: string) => void;
  panelFasta: string;
  setPanelFasta: (v: string) => void;

  /** The design criteria sent to the API. The criteria step binds to this. */
  params: DesignParams;
  updateParams: (patch: Partial<DesignParams>) => void;

  /** The shRNA scaffold options. The construct step binds to this. */
  convertOptions: ConvertOptions;
  updateConvertOptions: (patch: Partial<ConvertOptions>) => void;

  inspect: InspectResponse | null;
  design: DesignResponse | null;
  convert: ConvertResponse | null;

  selectedRegion: string | null;
  setSelectedRegion: (r: string) => void;
  selectedCandidate: CandidateResult | undefined;
  selectedConstruct: ConstructResult | undefined;
  /** The selected candidate's construct when the engine withheld it. */
  selectedWithheld: ConstructResult | undefined;
  /** Full per-base enrichment for the selected region, when available (sample /
   * post-API-addition); undefined means the alignment shows deciding mismatches
   * only. */
  enrichment: AlleleEnrichment[] | undefined;

  /** The engine's health report; null when no engine is reachable. */
  engine: HealthResponse | null;
  offline: boolean;
  busy: boolean;
  error: string | null;

  runDesign: (override?: Partial<DesignParams>) => Promise<void>;
  runConvert: () => Promise<void>;
}

export const SessionContext = createContext<SessionValue | null>(null);
