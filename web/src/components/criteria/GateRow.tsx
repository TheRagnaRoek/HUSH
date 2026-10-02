import { useState } from "react";
import type { ControlKind, GateDef } from "./gates";

/** Cosmetic controls for the scaffold. TODO: bind to DesignParams. */
function Control({ kind }: { kind: ControlKind }) {
  switch (kind) {
    case "rules":
      return (
        <span className="ctl">
          <select defaultValue="all" aria-label="rule sets">
            <option>all</option>
            <option>Ui-Tei</option>
          </select>
          <select defaultValue="any" aria-label="combine mode">
            <option>any</option>
            <option>all</option>
          </select>
        </span>
      );
    case "seedtm":
      return (
        <span className="ctl">
          <input type="number" defaultValue={21.5} step={0.5} aria-label="seed Tm" /> &deg;C
        </span>
      );
    case "homo":
      return (
        <span className="ctl">
          max <input type="number" defaultValue={3} aria-label="max run" />
        </span>
      );
    case "offtarget":
      return (
        <span className="ctl">
          min mm <input type="number" defaultValue={2} aria-label="min mismatch" />
        </span>
      );
    case "gc":
      return (
        <span className="ctl">
          <input type="number" defaultValue={30} aria-label="gc min" />&ndash;
          <input type="number" defaultValue={52} aria-label="gc max" />%
        </span>
      );
    case "hairpin":
      return (
        <span className="ctl">
          Tm &le; <input type="number" defaultValue={20} aria-label="hairpin Tm" />
        </span>
      );
    case "asym":
      return (
        <span className="ctl">
          win <input type="number" defaultValue={4} aria-label="asymmetry window" />
        </span>
      );
    case "duplex":
      return (
        <span className="ctl">
          <input type="number" defaultValue={19} aria-label="core length" />
          <select defaultValue="natural" aria-label="overhang">
            <option>natural</option>
            <option>dTdT</option>
          </select>
        </span>
      );
  }
}

export function GateRow({
  gate,
  on,
  onToggle,
}: {
  gate: GateDef;
  on: boolean;
  onToggle: (v: boolean) => void;
}) {
  const [open, setOpen] = useState(false);
  return (
    <div className="gate">
      <button
        type="button"
        className="switch"
        role="switch"
        aria-checked={on}
        aria-label={`Toggle ${gate.name}`}
        onClick={() => onToggle(!on)}
      />
      <div className="name">
        {gate.name}
        <span className="default">{gate.def}</span>
      </div>
      <div className="right">
        <Control kind={gate.control} />
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
