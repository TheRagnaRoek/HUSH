import type { ConstructResult } from "../../api/types";
import type { SourceKey } from "../../reference/sources";
import { RefNote } from "../common/RefLink";

function Chip({
  ok,
  label,
  value,
  info,
  cites,
  num,
}: {
  ok: boolean;
  label: string;
  value: string;
  info: string;
  cites: SourceKey[];
  num?: boolean;
}) {
  return (
    <span className={`qcchip has-tip ${ok ? "ok" : "warn"}`} tabIndex={0}>
      {label} <span className={num ? "v num" : "v"}>{value}</span>
      <span className="dtip" role="tooltip">
        <span className="dtt">{label}</span>
        {info}
        <RefNote cites={cites} to="qc" />
      </span>
    </span>
  );
}

const yesNo = (b: boolean) => (b ? "Yes" : "No");

export function QcChips({ construct }: { construct: ConstructResult }) {
  return (
    <div className="qc">
      <Chip
        ok={construct.folds_as_designed}
        label="Folds as designed"
        cites={["ge2010"]}
        value={yesNo(construct.folds_as_designed)}
        info="Whether the most stable predicted fold of the hairpin body is the intended one: the stem equals the guide core and the loop is the scaffold loop. No means a competing structure is more stable."
      />
      {construct.fold_dg37 != null && (
        <Chip
          ok
          num
          label={"ΔG₃₇"}
          cites={["xia1998", "groebe1988"]}
          value={construct.fold_dg37.toFixed(2)}
          info="Predicted free energy of the intended fold at 37 °C. The more negative it is, the more stable the hairpin."
        />
      )}
      {construct.fold_tm_c != null && (
        <Chip
          ok
          num
          label="Stem Tm"
          cites={["xia1998", "groebe1988"]}
          value={`${construct.fold_tm_c.toFixed(1)}°`}
          info="Predicted melting temperature of the hairpin fold. The higher it is, the more stable the stem."
        />
      )}
      <Chip
        ok={construct.starts_with_purine}
        label="Pol III 5′-A/G"
        cites={["gao2017"]}
        value={yesNo(construct.starts_with_purine)}
        info="Whether the transcript starts with A or G. A human U6 or 7SK promoter starts most efficiently on a purine at position +1. This is a flag only; the construct is kept either way."
      />
      <Chip
        ok={construct.premature_terminator_free}
        label="No premature terminator"
        cites={["bogenhagen1981", "gao2018"]}
        value={yesNo(construct.premature_terminator_free)}
        info="Whether the hairpin body is free of a run of four or more U (T in the DNA). Pol III stops at such a run and would cut the hairpin short. The planned poly-U tail is not counted."
      />
      <Chip
        ok={construct.max_homopolymer_run <= 3}
        num
        label="Max run"
        cites={[]}
        value={String(construct.max_homopolymer_run)}
        info="The longest run of a single nucleotide in the hairpin body, not counting the poly-U tail. A run above 3 is flagged."
      />
    </div>
  );
}
