import { useMemo } from "react";
import { FastaInput } from "../common/FastaInput";
import { validateFasta } from "../../lib/fasta";
import { useSession } from "../../session/useSession";

export function InputStep() {
  const { inspect, targetFasta, setTargetFasta, panelFasta, setPanelFasta, setStep } = useSession();

  const tSum = useMemo(
    () => (targetFasta.trim() ? validateFasta(targetFasta, "target") : null),
    [targetFasta],
  );
  const pSum = useMemo(
    () => (panelFasta.trim() ? validateFasta(panelFasta, "panel") : null),
    [panelFasta],
  );

  // Prefer a live preview of the user's own input; fall back to the sample.
  const usingOwn = tSum != null;
  const targetId = usingOwn ? (tSum.records[0]?.id ?? "-") : (inspect?.target.id ?? "-");
  const moltype = usingOwn ? tSum.moltype : inspect?.target.is_nucleotide ? "nucleotide" : "-";
  const length = usingOwn ? tSum.length : inspect?.target.length;
  const gc = usingOwn ? tSum.gcPercent : inspect?.target.gc_percent;
  const panelLabel = usingOwn
    ? pSum
      ? `loaded · ${pSum.recordCount}`
      : "none"
    : inspect?.panel
      ? `${inspect.panel.name} · ${inspect.panel.count}`
      : "none";
  const ready = usingOwn ? tSum.ok : (inspect?.target_ready ?? false);
  const gcRich = (gc ?? 0) >= 60;
  const blockContinue = usingOwn && !tSum.ok;

  return (
    <div className="workspace">
      <div className="whead">
        <h2>Load a target, and optionally a panel</h2>
        <p>
          Give HUSH the transcript you want to silence. Add a panel of related sequences - an HLA
          locus, a gene family - to turn on allele-aware targeting.
        </p>
      </div>

      <div className="intake">
        <div className="sources">
          <FastaInput
            kind="target"
            icon={"\u{1F9EC}"}
            title="Target sequence"
            hint="Click to browse, drag a FASTA here, or paste. mRNA/CDS; a protein is flagged and must be resolved to its transcript."
            value={targetFasta}
            onChange={setTargetFasta}
          />
          <FastaInput
            kind="panel"
            icon={"☷"}
            altAccent
            title={
              <>
                Related-sequence panel{" "}
                <span style={{ color: "var(--muted)", fontWeight: 400 }}>(optional)</span>
              </>
            }
            hint="FASTA of the allele family. Its presence engages the MSA / allele stage."
            value={panelFasta}
            onChange={setPanelFasta}
          />
        </div>

        <div className="inspect">
          <h4>Inspect {usingOwn ? "" : "(sample)"}</h4>
          <dl className="kv">
            <dt>Target</dt>
            <dd className="mono">{targetId}</dd>
            <dt>Molecule</dt>
            <dd>{moltype}</dd>
            <dt>Length</dt>
            <dd className="num">{length != null ? `${length} nt` : "-"}</dd>
            <dt>GC content</dt>
            <dd className="num">{gc != null ? `${gc}%` : "-"}</dd>
            <dt>Panel</dt>
            <dd className="mono">{panelLabel}</dd>
            <dt>Readiness</dt>
            <dd className={ready ? "ready" : ""}>{ready ? "ready" : usingOwn ? "fix errors" : "-"}</dd>
          </dl>
        </div>
      </div>

      {gcRich && moltype === "nucleotide" && (
        <>
          <div className="divider" />
          <div className="callout cyan">
            <span className="ic">{"ⓘ"}</span>
            <span>
              <b>GC-rich target.</b> At {gc}% GC this behaves like a classic HLA case - the seed
              filter may reject every candidate. HUSH will suggest relaxing it in the criteria step.
            </span>
          </div>
        </>
      )}

      <div className="actions">
        <button
          type="button"
          className="hbtn"
          disabled={blockContinue}
          onClick={() => setStep("targeting")}
        >
          Continue &rarr;
        </button>
        {blockContinue && <span className="hint">Fix the target sequence errors to continue.</span>}
      </div>
    </div>
  );
}
