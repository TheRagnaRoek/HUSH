import type { ReactNode } from "react";
import type { Step } from "../../session/context";
import { useSession } from "../../session/useSession";

/** Heidi, the design wizard. Her line changes per step. */
const LINES: Record<Step, ReactNode> = {
  input: (
    <>
      <b>Step 1.</b> Drop in the transcript to silence. Add a panel to unlock allele-aware
      targeting - I&apos;ll flag if it&apos;s GC-rich.
    </>
  ),
  targeting: (
    <>
      <b>Step 2.</b> Tell me what to hit and what to spare. Sparing something turns on
      discrimination ranking automatically.
    </>
  ),
  criteria: (
    <>
      <b>Step 3.</b> These are best-practice defaults. Open any <b>?</b> to learn a gate - or take
      my one-click HLA preset.
    </>
  ),
  candidates: (
    <>
      <b>Hi, I&apos;m Heidi.</b> Pick a candidate on the left and I&apos;ll walk you through why
      it&apos;s a strong pick.
    </>
  ),
  construct: (
    <>
      <b>Step 5.</b> Here&apos;s the cloneable construct. Green QC chips are good; amber ones are
      worth a look before you order oligos.
    </>
  ),
};

export function GuideBar() {
  const { step } = useSession();
  return (
    <div className="guide">
      <img className="heidi" src="/heidi.png" alt="Heidi, the HUSH design wizard" />
      <span className="gt">{LINES[step]}</span>
    </div>
  );
}
