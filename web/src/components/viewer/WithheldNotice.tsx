import type { ConstructResult } from "../../api/types";

/** Why the engine withheld a construct, with the offending U run marked in its hairpin. */
export function WithheldNotice({
  construct,
  terminator,
}: {
  construct: ConstructResult;
  terminator: string;
}) {
  const hairpin = construct.hairpin_rna_5to3;
  // The planned poly-U tail is meant to terminate; only the body is checked.
  const tail = terminator && hairpin.endsWith(terminator) ? terminator : "";
  const body = hairpin.slice(0, hairpin.length - tail.length);
  const parts = body.split(/(U{4,})/);

  return (
    <div className="callout warn">
      <span className="ic">{"⚠"}</span>
      <span>
        <b>Withheld.</b> This hairpin has a run of four or more U in its body (marked below). Pol
        III stops at such a run, so the transcript would be cut short before the hairpin is
        complete. Pick another candidate, try another loop or the natural overhang, or untick
        &quot;Withhold premature terminators&quot; in the scaffold options and convert again to
        inspect it anyway.
        <span className="mono hp">
          {parts.map((p, i) => (i % 2 === 1 ? <mark key={i}>{p}</mark> : p))}
          {tail && <span className="hp-tail">{tail}</span>}
        </span>
      </span>
    </div>
  );
}
