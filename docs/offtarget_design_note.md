# Design note: the transcriptome homology off-target filter

Status: **implemented.** `core/offtarget.py` provides HUSH's homology off-target
search - the counterpart to siDirect's transcriptome BLAST step, distinct from the
seed-Tm proxy in `core/seed.py`. This note records the decisions behind it so they
do not have to be re-derived.

Date of decision: 2026-08-15.

---

## 1. Why bowtie1, and why not BLAST

The task is exact: for each guide, enumerate every near-complementary site in a
reference transcriptome at 0-3 mismatches over the full length, on the forward
strand, and report the nearest one plus per-level counts - siDirect's "minimum
number of mismatches against off-targets" column.

- **siDirect uses GGGenome**, a fast k-mismatch full-text (FM-index) search. It is
  not open, which is why HUSH does not depend on it. The open, offline equivalent
  of a k-mismatch FM-index search is **bowtie1 in `-v` mode** (`bowtie -v N -a`):
  end-to-end alignment with at most N mismatches, all hits reported. This maps 1:1
  onto the requirement.
- **BLAST was rejected.** `blastn-short` is a *local heuristic*: it seeds on exact
  words and prunes by E-value, so a 19-mer with scattered mismatches is not
  guaranteed to be found, and it reports local HSPs rather than a clean full-guide
  mismatch count. It cannot reproduce siDirect's mismatch column reliably. The
  README's earlier phrase "siDirect BLASTs each guide" was loose wording.
- **bowtie1, not bowtie2.** bowtie2 replaced the fixed-mismatch `-v` mode with a
  local/gapped, quality-aware model, which is the wrong shape here. bowtie1 keeps
  `-v` (capped at 3 mismatches, which is enough). A convenient side effect: bowtie1
  is under the Artistic License, lighter than bowtie2's GPLv3.

## 2. Dependency posture

bowtie is a **binary**, and HUSH is otherwise pip-only. This is contained:

- The binary is reached via stdlib `subprocess`, so `import core.offtarget` never
  needs it - there is no import-time dependency.
- HUSH invokes bowtie as a **separate process** and does not link, bundle, or
  redistribute it. bowtie's license therefore does not propagate; HIDE can license
  HUSH however it wishes. bioconda is only a distribution channel and imposes
  nothing. (Standard footing for bioinformatics orchestrators.)
- The pure-Python `BruteForceSearcher` means the core library, the tests, and the
  small-panel / allele-aware path all work with **zero** binary dependency. bowtie
  is only needed for transcriptome-scale search, which - given the reference is
  hundreds of MB - is a server-side or desktop-bundled job anyway, never a
  fully-in-browser one.

## 3. Reference: RefSeq human mRNA, pinned

RefSeq (open) over Ensembl/GENCODE because its NM_/NR_ accessions match the
benchmark fixtures already in the repo, making on-target exclusion by accession and
the siDirect column comparison exact. `tools/build_offtarget_db.py` downloads the
`H_sapiens/mRNA_Prot` RNA set, builds the bowtie index, and writes a manifest
(RefSeq release, source URLs, sha256, bowtie version, build date). Point its
`--base-url` at an archived release directory for a permanent pin.

## 4. Strand semantics

A guide silences transcripts complementary to it, i.e. those carrying its target
site - the **sense** core. The passenger, if mis-loaded, silences those carrying
the **antisense** core. So each candidate contributes two reads (sense for the
guide side, antisense for the passenger side), searched forward-strand-only
(`--norc`) because transcripts are single-stranded sense mRNA. Both strands are
scored because either can be RISC-loaded (the same reason `core.seed` scores both).

## 5. Gate semantics and calibration

Reject if the nearest non-self off-target has fewer than `offtarget_min_mismatch`
mismatches. Default **2**, calibrated on the PDCD1 benchmark: siDirect keeps
candidates whose nearest off-target is >= 2 mismatches away and drops only the
near-perfect (0/1) ones. A threshold of 3 would drop siDirect's own 2-mismatch
picks, so 2 is the value that reproduces its decision. Search depth defaults to 3
(the bowtie `-v` cap) so the reported nearest-mismatch and counts line up with
siDirect's columns. The gate is flat, not position-weighted - seed weighting is
already `core.seed`'s job.

## 6. On-target exclusion

Hits to the candidate's own source transcript are dropped, matching siDirect
(the whole source accession, version-insensitive, not just the on-target position).
The accession is auto-detected from the target FASTA header (a RefSeq accession
regex; e.g. `NM_005018.3` in the PDCD1 fixture) and overridable with
`--offtarget-self-accession`.

## 7. Allele-aware caveat (deferred)

For a polymorphic family like HLA, intra-family cross-reactivity can be the design
goal (pan-allele knockdown), not a hazard. A single global gate that treats every
non-source hit as off-target cannot express that. v1 excludes only the source
accession; the "intended vs unintended cross-reactivity" distinction belongs with
the `msa/` allele module (the same point raised in
[`seed_offtarget_design_note.md`](seed_offtarget_design_note.md), section 4). The
`BruteForceSearcher` over a small curated panel is the natural tool for the
targeted version of the question.

## 8. Out of scope

Genomic (non-transcript) off-targets; seed-mediated scoring beyond siDirect
(`core.seed`); the expression-weighted expected-off-target-load model proposed in
the seed design note (a separate, later research task with its own validation gate).

## References

- Naito, Y., Yoshimura, J., Morishita, S., Ui-Tei, K. (2009). siDirect 2.0.
  *BMC Bioinformatics* 10, 392. doi:10.1186/1471-2105-10-392
- Langmead, B., Trapnell, C., Pop, M., Salzberg, S.L. (2009). Ultrafast and
  memory-efficient alignment of short DNA sequences to the human genome.
  *Genome Biology* 10, R25. doi:10.1186/gb-2009-10-3-r25
- O'Leary, N.A. et al. (2016). Reference sequence (RefSeq) database at NCBI.
  *Nucleic Acids Research* 44(D1), D733-D745. doi:10.1093/nar/gkv1189
