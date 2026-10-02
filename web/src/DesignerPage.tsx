import type { ReactElement } from "react";
import { CandidatesStep } from "./components/steps/CandidatesStep";
import { ConstructStep } from "./components/steps/ConstructStep";
import { CriteriaStep } from "./components/steps/CriteriaStep";
import { InputStep } from "./components/steps/InputStep";
import { TargetingStep } from "./components/steps/TargetingStep";
import type { Step } from "./session/context";
import { useSession } from "./session/useSession";

const COMPONENTS: Record<Step, () => ReactElement | null> = {
  input: InputStep,
  targeting: TargetingStep,
  criteria: CriteriaStep,
  candidates: CandidatesStep,
  construct: ConstructStep,
};

export function DesignerPage() {
  const { step, expert } = useSession();

  // Expert mode stacks every step on one scrollable page; guided shows one.
  if (expert) {
    return (
      <>
        <InputStep />
        <TargetingStep />
        <CriteriaStep />
        <CandidatesStep />
        <ConstructStep />
      </>
    );
  }

  const Current = COMPONENTS[step];
  return <Current />;
}
