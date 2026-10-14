import { useState } from "react";
import type { DesignParams } from "../../api/types";
import { rulesEnabled, type GateDef } from "./gates";

type Update = (patch: Partial<DesignParams>) => void;

/** Number input that writes back only valid numbers (clearing keeps the prior
 * value rather than pushing NaN into params). */
function NumberInput({
  value,
  onCommit,
  step,
  label,
  disabled,
}: {
  value: number;
  onCommit: (n: number) => void;
  step?: number;
  label: string;
  disabled?: boolean;
}) {
  return (
    <input
      type="number"
      step={step}
      aria-label={label}
      value={value}
      disabled={disabled}
      onChange={(e) => {
        const n = e.target.valueAsNumber;
        if (!Number.isNaN(n)) onCommit(n);
      }}
    />
  );
}

function GateControl({
  gate,
  params,
  update,
  expert,
}: {
  gate: GateDef;
  params: DesignParams;
  update: Update;
  expert: boolean;
}) {
  switch (gate.control) {
    case "rules":
      return (
        <span className="ctl">
          <label className="ruleset">
            <input
              type="checkbox"
              checked={params.use_ui_tei ?? true}
              onChange={(e) => update({ use_ui_tei: e.target.checked })}
            />{" "}
            UT
          </label>
          <label className="ruleset">
            <input
              type="checkbox"
              checked={params.use_reynolds ?? true}
              onChange={(e) => update({ use_reynolds: e.target.checked })}
            />{" "}
            Rey
          </label>
          <label className="ruleset">
            <input
              type="checkbox"
              checked={params.use_amarzguioui ?? true}
              onChange={(e) => update({ use_amarzguioui: e.target.checked })}
            />{" "}
            Am
          </label>
          <select
            aria-label="combine mode"
            value={params.rule_mode ?? "any"}
            onChange={(e) => update({ rule_mode: e.target.value as "any" | "all" })}
          >
            <option value="any">any</option>
            <option value="all">all</option>
          </select>
        </span>
      );
    case "seedtm":
      return (
        <span className="ctl">
          <NumberInput
            value={params.seed_tm_threshold ?? 21.5}
            onCommit={(n) => update({ seed_tm_threshold: n })}
            step={0.5}
            label="seed Tm threshold"
          />{" "}
          &deg;C
        </span>
      );
    case "homo":
      return (
        <span className="ctl">
          max{" "}
          <NumberInput
            value={params.max_homopolymer_run ?? 3}
            onCommit={(n) => update({ max_homopolymer_run: Math.round(n) })}
            label="max homopolymer run"
          />
        </span>
      );
    case "offtarget": {
      // Own-transcript resolution is expert-only: the defaults are the measured
      // ones, and a wizard user changing them silently breaks the gate.
      const resolve = params.resolve_offtarget_self ?? true;
      return (
        <span className="ctl">
          min mm{" "}
          <NumberInput
            value={params.offtarget_min_mismatch ?? 2}
            onCommit={(n) => update({ offtarget_min_mismatch: Math.round(n) })}
            label="min mismatch"
          />
          {expert && (
            <>
              <label className="ruleset">
                <input
                  type="checkbox"
                  checked={resolve}
                  onChange={(e) => update({ resolve_offtarget_self: e.target.checked })}
                />{" "}
                resolve own
              </label>
              <NumberInput
                value={Math.round((params.offtarget_self_coverage ?? 0.1) * 100)}
                onCommit={(n) => update({ offtarget_self_coverage: n / 100 })}
                label="own transcript k-mer coverage percent"
                disabled={!resolve}
              />
              % k-mers
            </>
          )}
        </span>
      );
    }
    case "gc":
      return (
        <span className="ctl">
          <NumberInput
            value={Math.round((params.gc_min ?? 0.3) * 100)}
            onCommit={(n) => update({ gc_min: n / 100 })}
            label="GC min percent"
          />
          &ndash;
          <NumberInput
            value={Math.round((params.gc_max ?? 0.52) * 100)}
            onCommit={(n) => update({ gc_max: n / 100 })}
            label="GC max percent"
          />
          %
        </span>
      );
    case "hairpin":
      return (
        <span className="ctl">
          Tm &le;{" "}
          <NumberInput
            value={params.hairpin_tm_max ?? 20}
            onCommit={(n) => update({ hairpin_tm_max: n })}
            label="hairpin Tm max"
          />
        </span>
      );
    case "asym":
      return (
        <span className="ctl">
          win{" "}
          <NumberInput
            value={params.asymmetry_window ?? 4}
            onCommit={(n) => update({ asymmetry_window: Math.round(n) })}
            label="asymmetry window"
          />
        </span>
      );
    case "duplex":
      return (
        <span className="ctl">
          <NumberInput
            value={params.length ?? 19}
            onCommit={(n) => update({ length: Math.round(n) })}
            label="core length"
          />
          <select
            aria-label="overhang"
            value={params.overhang ?? "natural"}
            onChange={(e) => update({ overhang: e.target.value })}
          >
            <option value="natural">natural</option>
            <option value="dtdt">dTdT</option>
            <option value="uu">UU</option>
          </select>
        </span>
      );
  }
}

export function GateRow({
  gate,
  params,
  update,
  expert = false,
}: {
  gate: GateDef;
  params: DesignParams;
  update: Update;
  expert?: boolean;
}) {
  const [open, setOpen] = useState(false);

  const enabled =
    gate.toggle === "rules"
      ? rulesEnabled(params)
      : gate.toggle === "simple" && gate.enabledKey
        ? Boolean(params[gate.enabledKey])
        : true;

  function toggle() {
    const next = !enabled;
    if (gate.toggle === "rules") {
      update({ use_ui_tei: next, use_reynolds: next, use_amarzguioui: next });
    } else if (gate.toggle === "simple" && gate.enabledKey) {
      update({ [gate.enabledKey]: next } as Partial<DesignParams>);
    }
  }

  return (
    <div className="gate">
      {gate.toggle === "none" ? (
        <span style={{ width: 40 }} aria-hidden="true" />
      ) : (
        <button
          type="button"
          className="switch"
          role="switch"
          aria-checked={enabled}
          aria-label={`Toggle ${gate.name}`}
          onClick={toggle}
        />
      )}
      <div className="name">
        {gate.name}
        <span className="default">{gate.def}</span>
      </div>
      <div className="right">
        <GateControl gate={gate} params={params} update={update} expert={expert} />
        <button
          type="button"
          className="why"
          aria-label={`Explain ${gate.name}`}
          aria-expanded={open}
          onClick={() => setOpen((v) => !v)}
        >
          ?
        </button>
      </div>
      <div className="sub">{gate.sub}</div>
      {open && (
        <div className="info">
          <b>What it does.</b> {gate.what}
          <div className="opt">{gate.opt}</div>
        </div>
      )}
    </div>
  );
}
