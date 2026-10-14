import type { ConstructResult } from "../../api/types";

function chip(ok: boolean, label: string, value: string) {
  return (
    <span className={`qcchip ${ok ? "ok" : "warn"}`} key={label}>
      {label} <span className="v">{value}</span>
    </span>
  );
}

export function QcChips({ construct }: { construct: ConstructResult }) {
  return (
    <div className="qc">
      {chip(construct.folds_as_designed, "folds as designed", construct.folds_as_designed ? "yes" : "no")}
      {construct.fold_dg37 != null && (
        <span className="qcchip ok">
          {"ΔG₃₇"} <span className="v num">{construct.fold_dg37.toFixed(2)}</span>
        </span>
      )}
      {construct.fold_tm_c != null && (
        <span className="qcchip ok">
          stem Tm <span className="v num">{construct.fold_tm_c.toFixed(1)}&deg;</span>
        </span>
      )}
      {chip(construct.starts_with_purine, "Pol III 5′-A/G", construct.starts_with_purine ? "yes" : "no")}
      {chip(
        construct.premature_terminator_free,
        "no premature term",
        construct.premature_terminator_free ? "yes" : "no",
      )}
      {chip(construct.max_homopolymer_run <= 3, "max run", String(construct.max_homopolymer_run))}
    </div>
  );
}
