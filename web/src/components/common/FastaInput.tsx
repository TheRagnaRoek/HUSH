import { useMemo, useRef, useState, type DragEvent, type ReactNode } from "react";
import { validateFasta } from "../../lib/fasta";

const MAX_BYTES = 10 * 1024 * 1024; // 10 MB
const ACCEPT = ".fasta,.fa,.fna,.ffn,.frn,.fas,.mpfa,.txt";

export function FastaInput({
  kind,
  icon,
  title,
  hint,
  value,
  onChange,
  altAccent,
}: {
  kind: "target" | "panel";
  icon: ReactNode;
  title: ReactNode;
  hint: string;
  value: string;
  onChange: (text: string) => void;
  altAccent?: boolean;
}) {
  const fileRef = useRef<HTMLInputElement>(null);
  const [drag, setDrag] = useState(false);
  const [pasteOpen, setPasteOpen] = useState(false);
  const [fileName, setFileName] = useState<string | null>(null);
  const [readError, setReadError] = useState<string | null>(null);

  const summary = useMemo(
    () => (value.trim() ? validateFasta(value, kind) : null),
    [value, kind],
  );

  function readFile(file: File) {
    setReadError(null);
    if (file.size > MAX_BYTES) {
      setReadError("File is larger than 10 MB.");
      return;
    }
    const reader = new FileReader();
    reader.onload = () => {
      onChange(String(reader.result ?? ""));
      setFileName(file.name);
    };
    reader.onerror = () => setReadError("Could not read the file.");
    reader.readAsText(file);
  }

  function onDrop(e: DragEvent) {
    e.preventDefault();
    setDrag(false);
    const file = e.dataTransfer.files?.[0];
    if (file) readFile(file);
  }

  return (
    <div className="fastafield">
      <div
        className={`drop${altAccent ? " alt" : ""}${drag ? " drag" : ""}`}
        onClick={() => fileRef.current?.click()}
        onDragOver={(e) => {
          e.preventDefault();
          setDrag(true);
        }}
        onDragLeave={() => setDrag(false)}
        onDrop={onDrop}
      >
        <span className="ic">{icon}</span>
        <div className="dropbody">
          <b>{title}</b>
          <small>{fileName ? `Loaded: ${fileName}` : hint}</small>
        </div>
        <div className="dropactions">
          <span className="pick">browse</span>
          <button
            type="button"
            className="pastelink"
            onClick={(e) => {
              e.stopPropagation();
              setPasteOpen((v) => !v);
            }}
          >
            {pasteOpen ? "hide" : "paste"}
          </button>
        </div>
      </div>

      <input
        ref={fileRef}
        type="file"
        accept={ACCEPT}
        hidden
        onChange={(e) => {
          const file = e.target.files?.[0];
          if (file) readFile(file);
          e.target.value = ""; // allow re-picking the same file
        }}
      />

      {pasteOpen && (
        <textarea
          className="fasta-textarea mono"
          placeholder={kind === "target" ? ">HLA-A*01:01\nATG..." : ">HLA-A*02:01\nATG...\n>HLA-A*03:01\nATG..."}
          value={value}
          onChange={(e) => {
            onChange(e.target.value);
            setFileName(null);
          }}
        />
      )}

      {readError ? (
        <div className="fstatus err">{readError}</div>
      ) : summary ? (
        summary.ok ? (
          <div className="fstatus ok">
            {"✓"} {summary.recordCount} sequence{summary.recordCount === 1 ? "" : "s"}
            {" · "}
            {summary.length} {summary.moltype === "protein" ? "aa" : "nt"} {"· "}
            {summary.moltype}
            {summary.gcPercent != null ? ` · ${summary.gcPercent}% GC` : ""}
            {summary.warnings.map((w) => (
              <span className="fwarn" key={w}>
                {" ⚠ "}
                {w}
              </span>
            ))}
          </div>
        ) : (
          <div className="fstatus err">
            {"✗"} {summary.errors[0]}
          </div>
        )
      ) : null}
    </div>
  );
}
