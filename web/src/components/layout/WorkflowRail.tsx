import { STEPS } from "../../session/context";
import { useSession } from "../../session/useSession";

export function WorkflowRail() {
  const { step, setStep } = useSession();
  const currentIndex = STEPS.findIndex((s) => s.id === step);
  return (
    <nav className="flow" aria-label="Design workflow">
      {STEPS.map((s, i) => {
        const current = s.id === step;
        const done = i < currentIndex;
        return (
          <button
            key={s.id}
            type="button"
            className={`step${done ? " done" : ""}`}
            aria-current={current || undefined}
            onClick={() => setStep(s.id)}
          >
            <span className="n">{i + 1}</span>
            <span className="lbl">{s.label}</span>
          </button>
        );
      })}
    </nav>
  );
}
