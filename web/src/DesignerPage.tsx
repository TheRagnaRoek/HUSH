import type { ReactElement } from "react";
import { ExpertPage } from "./components/expert/ExpertPage";
import { CandidatesStep } from "./components/steps/CandidatesStep";
import { ConstructStep } from "./components/steps/ConstructStep";
import { CriteriaStep } from "./components/steps/CriteriaStep";
import { InputStep } from "./components/steps/InputStep";
import { TargetingStep } from "./components/steps/TargetingStep";
import type { Step } from "./session/context";
import { useSession } from "./session/useSession";

const COMPONENTS: Record<Step, (props: { compact?: boolean }) => ReactElement | null> = {
  input: InputStep,
  targeting: TargetingStep,
  criteria: CriteriaStep,
  candidates: CandidatesStep,
  construct: ConstructStep,
};

export function DesignerPage() {
  const { step, expert } = useSession();

  // Expert mode is one overview page; guided shows one step at a time.
  if (expert) return <ExpertPage />;

  const Current = COMPONENTS[step];
  return <Current />;
}
