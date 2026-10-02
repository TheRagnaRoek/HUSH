import { useState } from "react";
import { useSession } from "../../session/useSession";

export function InputStep() {
  const { inspect, targetFasta, setTargetFasta, panelFasta, setPanelFasta, setStep } = useSession();
  const [openT, setOpenT] = useState(false);
  const [openP, setOpenP] = useState(false);
  const t = inspect?.target;
  const gcRich = (t?.gc_percent ?? 0) >= 60;

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
          <button type="button" className="drop" onClick={() => setOpenT((v) => !v)}>
            <span className="ic">{"\u{1F9EC}"}</span>
            <span>
              <b>Target sequence</b>
              <small>
                mRNA or CDS, as FASTA or pasted. A protein (PDB/mmCIF) is flagged and must be
                resolved to its transcript first.
              </small>
            </span>
            <span className="pick">browse / paste</span>
          </button>
          {openT && (
            <textarea
              className="os"
              style={{ width: "100%", minHeight: 90, marginBottom: 10 }}
              placeholder=">HLA-A*01:01&#10;ATG..."
              value={targetFasta}
              onChange={(e) => setTargetFasta(e.target.value)}
            />
          )}

          <button type="button" className="drop alt" onClick={() => setOpenP((v) => !v)}>
            <span className="ic">{"☷"}</span>
            <span>
              <b>
                Related-sequence panel{" "}
                <span style={{ color: "var(--muted)", fontWeight: 400 }}>(optional)</span>
              </b>
              <small>FASTA of the allele family. Its presence engages the MSA / allele stage.</small>
            </span>
            <span className="pick">browse / paste</span>
          </button>
          {openP && (
            <textarea
              className="os"
              style={{ width: "100%", minHeight: 90 }}
              placeholder=">HLA-A*02:01&#10;ATG..."
              value={panelFasta}
              onChange={(e) => setPanelFasta(e.target.value)}
            />
          )}
        </div>

        <div className="inspect">
          <h4>Inspect</h4>
          <dl className="kv">
            <dt>Target</dt>
            <dd className="mono">{t?.id ?? "-"}</dd>
            <dt>Molecule</dt>
            <dd>{t ? `${t.moltype}${t.is_nucleotide ? " (CDS)" : ""}` : "-"}</dd>
            <dt>Length</dt>
            <dd className="num">{t ? `${t.length} nt` : "-"}</dd>
            <dt>GC content</dt>
            <dd className="num">{t?.gc_percent != null ? `${t.gc_percent}%` : "-"}</dd>
            <dt>Panel</dt>
            <dd className="mono">{inspect?.panel ? `${inspect.panel.name} · ${inspect.panel.count}` : "none"}</dd>
            <dt>Alignment</dt>
            <dd>{inspect?.panel?.uniform ? "uniform" : "-"}</dd>
            <dt>Readiness</dt>
            <dd className={inspect?.target_ready ? "ready" : ""}>
              {inspect?.target_ready ? "ready" : "-"}
            </dd>
          </dl>
        </div>
      </div>

      {gcRich && (
        <>
          <div className="divider" />
          <div className="callout cyan">
            <span className="ic">{"ⓘ"}</span>
            <span>
              <b>GC-rich target.</b> At {t?.gc_percent}% GC this behaves like a classic HLA case -
              the seed filter may reject every candidate. HUSH will suggest relaxing it in the
              criteria step.
            </span>
          </div>
        </>
      )}

      <div className="actions">
        <button type="button" className="hbtn" onClick={() => setStep("targeting")}>
          Continue &rarr;
        </button>
      </div>
    </div>
  );
}
