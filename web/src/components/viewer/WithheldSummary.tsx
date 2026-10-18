import { useSession } from "../../session/useSession";

/** After a conversion, says how many constructs the engine withheld and why. */
export function WithheldSummary() {
  const { convert } = useSession();
  if (!convert || convert.withheld_count === 0) return null;
  return (
    <div className="banner">
      {convert.withheld_count} of {convert.converted_count} constructs were withheld: their hairpin
      body contains a run of four or more U, where Pol III would stop. They are marked
      &quot;Withheld&quot; in the candidate list. Another loop or the natural overhang usually
      avoids it.
    </div>
  );
}
