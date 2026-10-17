/**
 * Client-side FASTA parsing and pre-validation for the input step. The backend
 * (core.inputs) decides what is accepted.
 */

export interface FastaRecord {
  id: string;
  description: string;
  seq: string; // letters only, uppercased, whitespace stripped
}

export type MolType = "nucleotide" | "protein" | "unknown";

export interface FastaSummary {
  ok: boolean;
  records: FastaRecord[];
  recordCount: number;
  /** Length of the first record (the one a single-target load uses). */
  length: number;
  moltype: MolType;
  gcPercent?: number;
  errors: string[];
  warnings: string[];
}

const SEQ_CHARS = /^[A-Za-z*.-]+$/;
const NUC = new Set("ACGTUN");

function classify(seq: string): MolType {
  if (!seq) return "unknown";
  let nuc = 0;
  for (const ch of seq) if (NUC.has(ch)) nuc++;
  return nuc / seq.length >= 0.85 ? "nucleotide" : "protein";
}

function gcPercent(seq: string): number | undefined {
  let gc = 0;
  let acgt = 0;
  for (const ch of seq) {
    if (ch === "G" || ch === "C") {
      gc++;
      acgt++;
    } else if (ch === "A" || ch === "T" || ch === "U") {
      acgt++;
    }
  }
  if (acgt === 0) return undefined;
  return Math.round((gc / acgt) * 1000) / 10;
}

function parseRecords(text: string): { records: FastaRecord[]; errors: string[] } {
  const trimmed = text.trim();
  if (!trimmed) return { records: [], errors: ["No sequence provided."] };

  if (!trimmed.startsWith(">")) {
    if (trimmed.includes(">")) {
      return {
        records: [],
        errors: ['Sequence data appears before the first ">" header. A FASTA record must start with ">".'],
      };
    }
    // Bare sequence, no header.
    const seq = trimmed.replace(/\s+/g, "").toUpperCase();
    return { records: [{ id: "(unnamed)", description: "", seq }], errors: [] };
  }

  const records: FastaRecord[] = [];
  let cur: FastaRecord | null = null;
  for (const raw of trimmed.split(/\r?\n/)) {
    if (raw.startsWith(">")) {
      if (cur) records.push(cur);
      const header = raw.slice(1).trim();
      const sp = header.search(/\s/);
      cur = {
        id: sp === -1 ? header : header.slice(0, sp),
        description: sp === -1 ? "" : header.slice(sp + 1).trim(),
        seq: "",
      };
    } else if (cur) {
      cur.seq += raw.trim().toUpperCase();
    }
  }
  if (cur) records.push(cur);
  return { records, errors: [] };
}

export function validateFasta(text: string, kind: "target" | "panel"): FastaSummary {
  const { records, errors } = parseRecords(text);
  const warnings: string[] = [];

  if (errors.length > 0) {
    return { ok: false, records: [], recordCount: 0, length: 0, moltype: "unknown", errors, warnings };
  }

  records.forEach((r, i) => {
    const label = r.id && r.id !== "(unnamed)" ? `"${r.id}"` : `record ${i + 1}`;
    if (!r.seq) {
      errors.push(`${label} has a header but no sequence.`);
      return;
    }
    if (!SEQ_CHARS.test(r.seq)) {
      const bad = Array.from(new Set(r.seq.replace(/[A-Za-z*.-]/g, "")));
      errors.push(`${label} contains invalid characters: ${bad.slice(0, 6).join(" ")}`);
    }
  });

  if (kind === "target" && records.length > 1) {
    warnings.push(`${records.length} sequences found; the first is used as the target.`);
  }
  if (kind === "panel" && records.length === 1) {
    warnings.push("A panel usually holds 2+ related sequences.");
  }

  const first = records[0];
  const moltype = classify(first?.seq ?? "");
  if (moltype === "protein") {
    // The backend rejects this only at design time; warn early.
    errors.push(
      kind === "target"
        ? "This looks like a protein sequence. Design needs the coding transcript, so resolve it to its mRNA/CDS first."
        : "Panel members must be nucleotide sequences; this looks like protein.",
    );
  }
  return {
    ok: errors.length === 0,
    records,
    recordCount: records.length,
    length: first?.seq.length ?? 0,
    moltype,
    gcPercent: moltype === "nucleotide" ? gcPercent(first.seq) : undefined,
    errors,
    warnings,
  };
}
