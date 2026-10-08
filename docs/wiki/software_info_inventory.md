# Wiki software page - information inventory

Everything the repo currently documents, sorted into your page structure. This is
raw material, not prose: each bullet is a fact plus where it lives, so you can pull
what you need and ignore the rest.

Sources mined: `README.md`, `docker/README.md`, `docs/HUSH_concept_dossier.md`,
`docs/design_notes.md`, `docs/design_options.md`, `docs/offtarget_design_note.md`,
`docs/seed_offtarget_design_note.md`, `docs/shrna_design_note.md`,
`docs/offtarget_module_issue.md`, `pyproject.toml`, `docker/Dockerfile`,
`docker/environment.yml`, and every module docstring in `core/ msa/ shrna/ service/
app/ api/`.

Flags used below:
- **[GAP]** - something the page will want that the repo does not currently state.
- **[STALE]** - repo text that is out of date against the code.
- **[LATER]** - belongs to a chapter after Modules; parked here so it is not lost.

---

## Backstory

### The name and the project frame
- HUSH = **H**uman genome **U**niversal **S**ilencing **H**airpin-designer (`README.md:3`).
- Developed as part of **Project HIDE** (HLA Interference through Disguise and Erasure) (`README.md:16-17`).
- Self-description: open-source, from-sequence RNAi design tool; "modern, allele-aware siRNA/shRNA design pipeline for polymorphic gene families" (`README.md:3-11`).
- MIT licensed, copyright 2026 TheRagnaRoek (`LICENSE`).

### The field has been frozen since 2004
- Rational siRNA design "effectively frozen since 2004": the Ui-Tei, Reynolds and Amarzguioui rule sets - empirical thermodynamic-asymmetry heuristics - are still the backbone of nearly every modern tool (`dossier:8`).
- All three were defined and validated on the **19-nt duplex**, with position-specific criteria indexed against it (`design_notes.md:594-600`).

### siDirect 2.0: the one real advance, then abandoned
- siDirect 2.0 is the one implementation that meaningfully advanced the field beyond the 2004 rules: it added seed-duplex Tm-based off-target filtering (`dossier:8`).
- It has otherwise gone unmaintained: no shRNA output, no chemical-modification guidance, no modern interface, no path from siRNA candidate to a cloneable construct (`dossier:8`, `dossier:22`).
- siDirect relies on **GGGenome** for its transcriptome search, which is not open - the reason HUSH cannot simply reuse it (`offtarget_design_note.md:19-23`).
- siDirect's own caveat: it cannot design for **5.3% of RefSeq mRNAs** under default parameters, with multigene families named as the typical example. HLA is such a family (`seed_offtarget_design_note.md:122-125`).

### What was built instead, and why it does not replace siDirect
Two categories that do not overlap with siDirect's role (`dossier:10`):
- **Efficacy confirmation tools** (DeepSilencer, RNAiSpline and similar ML models): they *score* a given siRNA rather than generate candidates - they require a candidate to already exist (`dossier:10,27`).
- **Precomputed genome-scale resources** (BIOPREDsi): one scored siRNA per human gene, fixed and non-customizable, no allele or off-target-set flexibility (`dossier:10,28`).
- Neither replaces a flexible, from-sequence design pipeline (`dossier:10`).

### The prior-art table (ready to lift as a figure/table)
From `dossier:20-29` - tool, category, strength, gap:

| Tool | Category | Strength | Gap relative to the target |
|---|---|---|---|
| siDirect 2.0 | siRNA design | Seed-duplex Tm off-target filtering; still the most mechanistically sound classical tool | Unmaintained UI, no shRNA output, no chemistry guidance, no allele-awareness |
| DSIR | siRNA design | Purpose-built for Dicer-substrate siRNAs | Narrow scope, not general-purpose |
| SSD / siRNA Wizard | siRNA to shRNA conversion | Templates mRNA input directly into mono-/bi-/trivalent shRNA cassettes | Naive templating on outdated siRNA cores; no seed-Tm or scaffold-processing optimization |
| pssRNAit | Integrated RNAi pipeline (plant) | Chains rule-based + ML + genome-wide off-target scoring | Plant-specific biology and databases |
| SIREN (2025) | Integrated dsRNA pipeline | Thermodynamically-informed cumulative off-target modelling, modern and open-source | Built for long dsRNA, not discrete siRNA/shRNA |
| DeepSilencer / RNAiSpline | Efficacy confirmation (ML) | State-of-the-art scoring of a given candidate | Not generative |
| BIOPREDsi | Precomputed genome-scale | One scored siRNA per human gene | Fixed, non-customizable |
| HLA-Knockout (2026) | Allele-specific CRISPR design | Solves the HLA-polymorphism targeting problem - for CRISPR | No RNAi equivalent exists |

Extra detail on SSD worth having: the paper publishes **neither loop sequences nor folding rules**; it stamps its siRNA picks into mono/bi/trivalent cassettes (`shrna_design_note.md:72-75`).

### The HLA-specific hook
- HLA is among the most polymorphic loci in the human genome; that allelic diversity is what motivates the allele-aware stage (`README.md:13-15`).
- HLA's extreme polymorphism means a generic tool cannot simply be pointed at "the HLA gene": conserved-region targeting (broad/pan-allele) and allele-specific targeting both need sequence-set-aware logic that **no current RNAi tool implements** (`dossier:12`).
- The same problem was solved **for CRISPR** in 2026: HLA-Knockout (Rutgers) pulls allele sequences from IPD-IMGT/HLA and uses mismatch/PAM-disruption filtering for allele-specific knockout without disturbing non-target alleles. No RNAi equivalent exists (`dossier:14`).
- Prior HLA-knockdown shRNA work (lentiviral, pan-Class I and allele-specific) proves the biology works - but as bespoke, manually designed reagents, not a generalizable design tool (`dossier:14`).
- Application framing: hypoimmunogenic/allogeneic cell engineering and immune-evasion research (`dossier:12`).

### Dated decision record (useful for a development-timeline narrative)
- **2026-07-24** - the seed-Tm filter decision: keep siDirect-faithful, defer the replacement (`seed_offtarget_design_note.md:9`).
- **2026-08-15** - transcriptome off-target filter implemented; bowtie1 chosen over BLAST (`offtarget_design_note.md:8`, `offtarget_module_issue.md:3`).
- **2026-09-06** - deployment posture fixed: hosted, browser-based, server-side bowtie and index (`offtarget_design_note.md:50`).
- **2026-09-12** - shRNA scope fixed to the stem-loop scaffold; miR-30/miR-E explicitly descoped (`shrna_design_note.md:9,15-22`).

---

## Goal

### The thesis
- No existing tool combines a modernized accessible siRNA design core, native shRNA scaffold conversion, and allele-polymorphism-aware targeting in one open-source pipeline; and no RNAi tool addresses HLA's polymorphism problem at all, despite a directly analogous CRISPR tool existing and validated wet-lab precedent for HLA shRNA knockdown (`dossier:43`).
- Positioning claim: "a real, currently-unfilled gap, scoped to what a graduate team can credibly deliver and validate within an iGEM timeline" (`dossier:43`).
- Build a modernized, open-source **successor to siDirect's core logic** that closes three gaps simultaneously rather than one (`dossier:35`).

### The three pillars (`dossier:37-39`)
1. **Modernized design core** - keep the classical rule filter as a fast pre-screen, paired with a lightweight, *interpretable* ML scoring layer (logistic regression / gradient boosting trained on public benchmark data such as Huesken) rather than a deep architecture. Justification: feasible within the timeline, and per recent comparative work competitive with deep learning on the same benchmarks.
   - **[STALE]** The ML layer is **not implemented**. `scikit-learn>=1.4` is a declared dependency in `pyproject.toml:19` but is imported nowhere in `core/ msa/ shrna/ service/ app/ api/ tools/ tests/`. The page should present pillar 1 as "classical rules faithfully reimplemented and validated, ML layer designed but not yet built", or the dossier claim needs rewording.
2. **Native shRNA conversion** - seed-duplex-aware siRNA candidates carried through to *scored* shRNA scaffold design, rather than naive template-stamping. Dossier originally named stem-loop **and** miR-30/miR-E pri-miRNA context; miR-E was later descoped (see below). Status: **implemented** for stem-loop (`README.md:27`).
3. **Allele-awareness** - an MSA-driven module that aligns a user-selected allele panel and scores candidates for either pan-allele conserved-region targeting or single-allele-specific targeting. The tool is **sequence-agnostic**; HLA (via IPD-IMGT/HLA integration) is the proof-of-concept panel. Explicitly framed as "the RNAi-appropriate analog of what HLA-Knockout does for CRISPR" (`dossier:39`). Status: **implemented** (`README.md:26`).

### The design philosophy (the most quotable part of the repo)
- **Filter, not ranking.** No single combined score, because the signals (three rule sets, a GC window, a seed-duplex Tm, a loading-asymmetry sign) live on **incommensurable scales** and no defensible weighting exists without efficacy data HUSH does not have. Inventing weights would be a black box that needs tuning and validation; a transparent, option-rich filter does not (`design_notes.md:495-505`).
- **Every verdict is legible.** Each candidate keeps its full per-criterion breakdown; the run prints the criteria and a funnel of how many candidates clear each gate, so an empty result shows exactly which gate removed them (`design_notes.md:504-505`, `design_options.md:198-199`).
- **A published threshold travels with the parameter set it was fit on.** Hence two nearest-neighbor tables exposed by explicit name and no single "current" default (`design_notes.md:61-69`, `seed_offtarget_design_note.md:55-57`).
- **Fidelity over modernisation-for-its-own-sake.** The rules keep their published form; the tunable counterparts are separate standalone filters (`design_notes.md:544-556`).
- **Reproduce the incumbent exactly before claiming to beat it.** Three reasons given: it is the baseline to beat, it labels existing experimental data, and a published comparison will be asked to show the incumbent was implemented correctly (`seed_offtarget_design_note.md:25-34`).
- **No invented calibration.** Discrimination ranking is a lexicographic tuple of published facts, not a weighted score (`msa/discrimination.py` docstring).
- **Missing data is never evidence.** Coverage is required by both allele gates; folding "uncovered" into "mismatch" would turn missing data into evidence - "a silent error, not a design choice" (`design_notes.md:666-674`).
- **Best-practice defaults, not the strictest or the loosest** (`design_options.md:17-25`, `README.md:123`).

### Explicit non-goals / deferrals
- miR-30 / miR-E pri-miRNA scaffold: out of scope for the current timeline - demands research effort and wet-lab infrastructure (single-copy expression, backbone processing validation) the team does not have. Fellmann (2013) kept as context only (`shrna_design_note.md:15-22`).
- Genomic (non-transcript) off-targets (`offtarget_design_note.md:111`).
- The expression-weighted expected-off-target-load model (a separate, later research task with its own validation gate) (`offtarget_design_note.md:112-113`, `seed_offtarget_design_note.md:129-147`).
- A general sequence-motif constraint on the `design` step (the Pol III 5'-A/G start preference that motivated it is now a QC flag in the shRNA stage instead) (`design_options.md:201-210`).
- No asymmetry magnitude threshold - the filter uses the sign only (`design_options.md:151-153`).
- IUPAC ambiguity codes are not resolved in the MSA stage (`design_notes.md:752-753`).
- AgoshRNA regime: a boundary condition the module must decide deliberately rather than land on by accident (`shrna_design_note.md:60-68`).

---

## Use for the Community

### Open and reusable
- Open source, MIT licensed, pure-pip core (`LICENSE`, `README.md:39`).
- **Sequence-agnostic**: any panel of related sequences works - HLA is not baked in, it is the proof-of-concept panel (`README.md:10-15`).
- Three frontends over **one** pipeline, all driving it through the same `service.HushSession`, so the pipeline is wired together in exactly one place (`service/__init__.py` docstring):
  - one-shot CLI (`hush inspect` / `design` / `convert`) - the way to script HUSH;
  - interactive REPL (`hush shell`) - load once, tune criteria, re-run; alignment computed once and reused;
  - HTTP API (FastAPI) - stateless or session mode, with `/docs` and `/openapi.json` so a **typed TypeScript client can be generated** (`README.md:211-245`).
- Extension points designed in, each behind a named protocol or registry:
  - `AlignmentBackend` - a progressive aligner can replace the star backend without touching `variation`, `windows`, or the CLI (`design_notes.md:721-724`);
  - `OffTargetSearcher` - two interchangeable backends (pure-Python brute force, bowtie1) (`core/offtarget.py` docstring);
  - the loop registry - literature-backed loops, each with its citation (`shrna/loops.py`);
  - `Orientation` enum - modelled as an enum so antisense-loop-sense can be added later without changing the construct's shape (`shrna/convert.py`).

### Zero-friction for the end user
- Deployment posture (decided 2026-09-06): HUSH ships as a **hosted, browser-based tool with a server backend**, so bowtie and the RefSeq index are "our responsibility, never the end user's" (`offtarget_design_note.md:50-59`).
- Users only hit the web UI; they never install bowtie, download RefSeq, or run `bowtie-build` (`offtarget_design_note.md:55-57`, `docker/README.md:102-108`).
- The pure-Python core, the tests, and the small-panel / allele-aware path all run with **zero** binary dependency (`offtarget_design_note.md:45-47`).

### What a user actually gets out
- Four targeting intents, covering the full space of allele requests (`README.md:158-166`, `design_options.md:81-96`):
  - `--pan-allele` - conserved sites across the whole family (default when a panel is given);
  - `--allele-specific ID` - single one panel member out, spare the rest;
  - `--target-specific` - spare the whole panel, hit none of it (for a target loaded apart from its reference panel);
  - `--hit` (+ `--spare`) - explicit subsets.
- Per-candidate legibility: the `minDiff` and `disc` columns describe the **weakest** spared allele; `--explain` (or `assess <region>` in the shell) prints the full per-allele breakdown (`README.md:183-188`).
- The `disc` string encodes five published criteria in one cell, e.g. `P16 post H+ S pp` - guide position, region (post/seed/anch), Huang listing (H+/H-/H?), position tier (S/C/-), mismatch class (pp/pu/up/uu/?) (`design_options.md:108-124`).
- **Vector-agnostic** shRNA output: the bare hairpin insert as top/bottom DNA oligos plus a Pol III poly-T terminator - the wet-lab adds sticky-end tails for whichever vector they clone into (`README.md:192-199`, `shrna/convert.py`).
- Each construct is **folded and QC-flagged**, not template-stamped: `folds_as_designed`, `starts_with_purine` (U6 Pol III +1 preference), `premature_terminator_free` (no 4-U run in the hairpin body, including the loop/stem junctions the upstream homopolymer gate cannot see), `max_homopolymer_run` (`shrna/convert.py`).
- Outputs: terminal tables, TSV (`-o`), JSON (API). `save <dir>` in the shell writes `design.tsv` and `shrna.tsv` (`README.md:79`).

### Reproducibility as a community good
- Deterministic, offline, reproducible: a local reference DB with a pinned build (`offtarget_module_issue.md:80-81`).
- `tools/build_offtarget_db.py` writes a manifest recording RefSeq release, source URLs, sha256 checksums, bowtie version and build date; `--base-url` can point at an archived release directory for a permanent pin (`offtarget_design_note.md:63-68`).
- `tools/compare_sidirect.py` diffs a HUSH design TSV against a siDirect export, so the comparison against the gold standard is reproducible by anyone (`README.md:246-255`).
- Benchmark fixtures shipped in the repo: `data/siDirect_PDCD1-CDS.ods`, `data/pdcd1_hush_relaxed.tsv`, `tests/data/amarzguioui_table3.tsv` (`offtarget_module_issue.md:144-145`, `design_notes.md:481`).
- **[GAP]** No CITATION file, no Zenodo DOI, no public repo URL stated anywhere in the docs. The wiki will want at least a repo link, and a DOI is cheap and makes the tool citable.

### Honest limitations to state (better from you than from a judge)
- The seed filter is deliberately strict and can reject **every** candidate on a very GC-rich family such as HLA; the documented remedy is `--no-seed` (`README.md:141-143`, `design_options.md:23-25`).
- The star aligner is adequate while members are near-identical (an HLA locus runs ~98-99% identity) and degrades as they diverge; divergent families want a progressive backend (`design_notes.md:717-721`).
- API sessions are an in-memory warm cache, lost on restart, not a database (`README.md:244`).
- Gaps are placed without regard to reading frame - no codon awareness (`msa/align_star.py` docstring).
- Insertions are recorded but not scored (would make the guide bulge, a strong discriminating feature; needs mismatch thermodynamics `core.thermo` does not carry) (`design_notes.md:705-707`).

---

## Architecture

### The five-stage pipeline (`README.md:22-28`)

| Stage | Module | Status |
|---|---|---|
| 1. Input Layer | `core.inputs` + `core.sequence` | implemented |
| 2. siRNA Design | `core.sirna` / `rules` / `thermo` / `seed` / `offtarget` / `asymmetry` / `hairpin` / `selection` | implemented, validated against siDirect |
| 3. MSA / Allele | `msa` | implemented, wired into `hush design --panel` |
| 4. shRNA Conversion | `shrna` | implemented (stem-loop scaffold), wired into `hush convert` |
| 5. Output Layer | `core` / `app` | planned |

### Layering and the dependency rule
- Project layout, one directory per concern (`README.md:257-271`): `core/` domain model + siRNA engine; `msa/` allele stage; `shrna/` conversion; `service/` orchestration (`HushSession`); `app/` CLI + shell; `api/` HTTP; `data/`, `tools/`, `docs/`, `tests/`.
- **The dependency direction is `msa -> core`, and `core.selection` must never import `msa`.** The CLI (`app.cli`) is the composition root that may import both, and that is where the allele verdict is ANDed with the siRNA selection verdict (`design_notes.md:726-748`).
- Consequence, stated as the reason for the choice: the allele stage and the siRNA engine each stay testable in isolation, and the funnel reports the allele gate as just one more independent gate rather than a special case (`design_notes.md:745-748`).
- `msa/msa_errors.py` keeps the stage's own failure vocabulary but every type derives from `core.errors.HushError`, so "catch anything HUSH raised" keeps working and the stage stays catchable as a unit (`msa/msa_errors.py` docstring).
- The **service layer is the seam**: the pure pipeline lives in `core/msa/shrna`, the CLI and API are thin renderers, and `HushSession` sits between them holding loaded state and returning structured results - "never printed text, never JSON" (`service/__init__.py`).
- Results are **lossless carriers**: a result never knows about terminals, colours, or wire formats (`service/results.py`).
- `HushSession` raises domain exceptions (`HushError`, or `ValueError` for an incoherent request) and prints nothing; turning a result into a table or JSON, and a `ValueError` into an exit code or a 4xx, is the caller's job (`service/session.py`).
- The API mirrors the CLI pair-for-pair: `api.serialize` is "the counterpart to the CLI's printers" - same columns, as JSON (`api/serialize.py`).

### The architectural keystone: the coordinate contract
- Every panel member is projected onto **ungapped target coordinates**: a string exactly as long as the target, index *i* carrying whatever residue the member has at target position *i*, with columns the alignment inserted into the target removed (kept as insertion payloads, not dropped) (`design_notes.md:652-658`).
- Why it holds the stage together: window position *p* on the target is the same *p* that `core.sirna` enumerates a candidate at, so the two lists **join on position without either side knowing about the other**. No unit ever converts between alignment columns and target positions, because after projection there are no alignment columns left to convert (`design_notes.md:656-659`).
- `PanelAlignment` *enforces* the invariant (every row exactly the target's width) rather than trusting it, so a backend that violates it fails loudly at the hand-off instead of silently shifting which positions look polymorphic (`design_notes.md:661-664`).
- The window length is bound to the siRNA core length (`--length`) so the two lists stay one-for-one (`design_notes.md:739-741`).

### Cross-cutting structural decisions worth naming
- **Two thermodynamic entry points, deliberately separate** (`duplex_thermo()` vs `stacking_terms()`): using the wrong one silently shifts every number by the initiation term, so the choice is made explicit at the call site rather than hidden behind a flag (`design_notes.md:47-59`).
- **Ten stored NN values, six derived by symmetry** (`_expand_by_symmetry`), so no value is typed twice and the two spellings of a step can never disagree (`design_notes.md:71-81`).
- **The overhang is modelled as a separate feature**, tracked apart from the core, and **no filter ever sees it** (`design_notes.md:596-600`).
- **Two off-target angles kept orthogonal**: `core.seed` asks a sequence-intrinsic question and never looks at a transcript; `core.offtarget` asks a database question. They compose - a candidate can be seed-clean yet have a near-perfect full-length off-target, or vice versa (`design_notes.md:196-205`). The homology gate is deliberately **flat, not position-weighted**, because seed weighting is already `core.seed`'s job (`design_notes.md:260-261`).
- **Four intents, one structure**: pan-allele / allele-specific / target-specific / subset are not four types but one structure (a *hit* set and a disjoint *spare* set) with different values, so a single code path serves all of them. `spare` is stored explicitly rather than inferred as "whatever else is in the panel", because a panel is data that can change (`design_notes.md:680-691`). `TargetingMode` is *derived* from the two sets rather than stored, so a reported mode can never disagree with its contents (`msa/targeting.py`).
- **Gate and ordering are separate layers**: the gate ranks nothing; `msa.discrimination` orders survivors but never admits or rejects one (`design_notes.md:693-707`).
- **One bowtie process per run, not per candidate**: `precompute()` aligns every strand of every candidate in a single bowtie run and returns a lookup-backed `PrecomputedSearcher`, so the per-candidate selection loop pays no process cost (`design_notes.md:232-238`).
- **The brute-force searcher is the oracle** the bowtie backend is checked against; both engines are ungapped, Hamming-distance and forward-orientation, so their results agree exactly (`design_notes.md:229-238`).

### Compatibility

#### Language and dependencies (`pyproject.toml`)
- Python **>= 3.11**; `target-version = "py311"`, ruff line length 100.
- Core runtime dependencies, all pip: `biopython>=1.83`, `numpy>=1.26`, `pandas>=2.0`, `scikit-learn>=1.4`, `typer>=0.12`.
  - **[STALE]** `scikit-learn` is declared but imported nowhere - it is a placeholder for the unbuilt ML layer.
- Optional extras: `[api]` = `fastapi>=0.110`, `uvicorn>=0.29`; `[dev]` = `pytest>=8.0`, `pytest-cov>=5.0`, `httpx>=0.27`.
- Declared packages: `core`, `msa`, `shrna`, `service`, `app`, `api`. Console script: `hush = app.cli:app`. Version `0.0.1`.
- `core.sequence` is deliberately **dependency-light** - pure standard library, no Biopython - because every stage speaks in its types (`core/sequence.py`).
- Biopython is confined to low-level file parsing in `core.inputs`, which adapts its records into validated domain objects (`core/inputs.py`).

#### The one non-pip piece, and how it is contained
- **bowtie1**, used only for transcriptome-scale off-target search (`README.md:40-47`).
- Reached via stdlib `subprocess`, so `import core.offtarget` never needs the binary - **no import-time dependency** (`offtarget_design_note.md:38-39`).
- HUSH invokes bowtie as a **separate process** and does not link, bundle or redistribute it, so **bowtie's license does not propagate** - HIDE can license HUSH however it wishes. bioconda is only a distribution channel and imposes nothing. Called "standard footing for bioinformatics orchestrators" (`offtarget_design_note.md:40-43`).
- Licensing detail: bowtie1 is under the **Artistic License**, lighter than bowtie2's GPLv3 (`offtarget_design_note.md:32`, `design_notes.md:220-224`).
- Why bowtie1 and not bowtie2: bowtie2 replaced the fixed-mismatch `-v` mode with a local/gapped, quality-aware model, which is the wrong shape; bowtie1 keeps `-v` (capped at 3 mismatches, which is enough) (`offtarget_design_note.md:29-32`).
- Why not BLAST: `blastn-short` is a *local heuristic* that seeds on exact words and prunes by E-value, so a 19-mer with scattered mismatches is not guaranteed to be found, and it reports local HSPs rather than a clean full-guide mismatch count - it cannot reproduce siDirect's mismatch column reliably (`offtarget_design_note.md:24-28`).

#### Containerisation (`docker/`)
- Base image `mambaorg/micromamba:1.5.8`; conda env pins `python=3.11` and `bowtie=1.3.1` from bioconda/conda-forge; HUSH pip-installed on top (`docker/Dockerfile`, `docker/environment.yml`).
- **One image, two jobs**: build the RefSeq index once, and serve designs from the same image (`docker/README.md:8-13`).
- Because it is a Linux container it removes the friction of installing bowtie1 via bioconda on **native Windows** - everything runs the same on any host with Docker (`docker/README.md:15-17`, `offtarget_design_note.md:57-58`).
- Index build: `docker compose run --rm build-index`; output lands in the `hush_ref` Docker volume alongside `refseq_human_mrna.manifest.json`. Network-heavy, runs once (`docker/README.md:19-31`).
- Ready-made CLI wrappers that need nothing installed on the host: `docker/hush.ps1` (PowerShell) and `docker/hush.sh` (macOS/Linux). Both mount the current directory at `/work` and the reference volume at `/data/reference` (`docker/README.md:33-55`).
- Overridable via `HUSH_IMAGE` and `HUSH_REF_VOLUME` env vars (volume defaults to `hush_hush_ref`) (`docker/README.md:58-62`).
- Shell-function recipes for both PowerShell and bash/zsh so the wrapper "looks like a plain `hush` command" (`docker/README.md:64-77`).
- The wrappers pass an interactive terminal through when run from a console so the **REPL stays open**, and drop the interactive flags automatically when input/output is redirected (a pipe, or CI) (`docker/README.md:96-99`).
- Dockerfile detail worth noting as a portability fix: the reference dir is pre-created owned by the non-root run user, because otherwise a fresh named volume mounts root-owned and the index build fails with EACCES (`docker/Dockerfile`).
- Non-container launchers: `./hush` (POSIX) and `hush.bat` (Windows) call the venv's interpreter directly, so `hush` works from the project folder **without activating the venv** (`README.md:49-55`).

#### Input and output compatibility
- Input formats: **FASTA** (primary - a target mRNA/CDS, or a multi-record allele panel), **PDB**, **mmCIF**; `--panel-format` can override autodetection (`core/inputs.py`, `design_options.md:84`).
- A **protein** target (e.g. parsed from a PDB structure) is parsed and labelled faithfully as protein but flagged: siRNA/shRNA design needs the coding transcript, enforced downstream via `DesignInput.require_nucleotide_target` (`README.md:106-108`, `core/inputs.py`).
- The API takes target and panel as **FASTA text in the JSON body**, not files; a target may be a full FASTA record or a bare sequence, a panel must be FASTA since allele targeting addresses members by id (`api/inputs.py`).
- Output: terminal tables, TSV, JSON. API enum-valued fields are plain strings on the wire so the client never has to know the Python enums (`api/schemas.py`).
- API host/port/reload read from `HUSH_API_HOST`, `HUSH_API_PORT`, `HUSH_API_RELOAD` so a container can point it anywhere without code changes (`api/__main__.py`).

#### Reference data
- **RefSeq human mRNA**, pinned per build, chosen over Ensembl/GENCODE because its NM_/NR_ accessions match the benchmark fixtures already in the repo, making on-target exclusion by accession and the siDirect column comparison exact (`offtarget_design_note.md:61-68`).
- `tools/build_offtarget_db.py` downloads the `H_sapiens/mRNA_Prot` RNA set, runs `bowtie-build`, and writes the manifest (`offtarget_design_note.md:65-68`).
- Reference data lives under `data/reference/`, gitignored - the transcriptome is never committed (`offtarget_module_issue.md:100-102`).
- Hosting: publish the volume's `.ebwt` shards plus manifest to a mounted disk, object store, or a release asset the runtime downloads on startup (`docker/README.md:102-108`).

### Modules

Line counts are current and give a sense of weight. Every module's own docstring carries its contract and scope; `docs/design_notes.md` carries the derivations.

#### `core/` - domain model and the siRNA design engine
- **`core.sequence`** (323 lines) - the domain model every stage speaks: `MolType` (DNA/RNA/protein/unknown), `BioSequence` (immutable, validated, with transcription / reverse-complement / GC content), `SequencePanel` (a named set of related sequences), `DesignInput` (target + optional panel). Pure stdlib, heavily validated.
- **`core.inputs`** (191 lines) - stage 1, the Input Layer. FASTA/PDB parsing into the domain model; Biopython does low-level parsing, this module adapts and validates.
- **`core.sirna`** (243 lines) - stage 2 substrate. Turns a target mRNA into the full set of sliding-window candidates; **nothing is scored here**. Canonical duplex: 19-nt fully paired core with a 2-nt 3' overhang on each strand (21-nt synthesized strands, the Elbashir/Tuschl design). `OverhangPolicy`: `natural` (default, mRNA-derived flanks - the true *offset* Tuschl duplex, the two strands' overhangs come from opposite flanks so they are not a blunt reverse-complement pair), `dtdt` (lowercase `tt` to flag the DNA residues), `uu`. A candidate whose natural flank runs past a transcript end falls back to `dtdt`.
  - Worth stating on the wiki: **19 / 21 / 23 get conflated in the literature** - 19 is the paired core, 21 is the strand a vendor synthesizes, 23 is a different measurement again (an mRNA target-site footprint, or a Dicer-substrate design) (`design_notes.md:586-592`).
- **`core.thermo`** (379 lines) - the nearest-neighbor primitive every other unit builds on: sequence to `(dH, dS)`, hence dG at any temperature and Tm. Two entry points (`duplex_thermo()` / `stacking_terms()`), two parameter sets by explicit name (Xia 1998 default, Freier 1986). Watson-Crick only, 1 M NaCl, no salt correction. Scope deferral: mismatch, G-U wobble and dangling-end parameters are not included - allele-discrimination scoring will need them.
- **`core.rules`** (461 lines) - the three 2004 rule sets.
  - **Ui-Tei**: strict 4-rule pass/fail (A/U at guide 5' terminus; G/C at sense 5' terminus; >= 4 A/U in the guide's 5'-terminal 7 nt; no G/C stretch longer than 9 nt); "highly effective" requires all four.
  - **Reynolds**: scored **out of 10** following the original algorithm, not the common "8 equal binary criteria" simplification. +1 each for GC 30-52%, absence of internal repeat, A at 19, A at 3, U at 10; **+1 per A/U at positions 15-19** (additive 0 to +5, not a threshold); -1 for G at 13 and for G/C at 19. Range -2 to +10, >= 6 predicts functional. Position 19 carries three criteria at once, giving a clean A(+2) / U(+1) / G,C(-1) gradient. For cores longer than 19 nt the 3'-anchored parts follow the actual 3' terminus while internal positional parts stay at their absolute positions from the 5' end.
  - **Amarzguioui & Prydz**: additive motif score - +1 for G/C at 1, A at 6, A/U at 19; -1 for U at 1, G at 19 - **plus the raw A/U differential itself** (the paper's "DT3", over the three terminal base pairs at each end, sense strand reference, 3' end minus 5' end, range -3 to +3), not a lookup. Threshold >= 3.
- **`core.seed`** (223 lines) - the siDirect 2.0 seed-duplex Tm off-target filter, reproduced exactly: Freier 1986 parameters, initiation constant -10.8, 100 uM strands, 100 mM Na+ via `16.6*log10[Na]`, gas constant **1.987** (the paper's rounded value), no terminal A-U penalty, no symmetry correction, threshold Tm < 21.5 C on **both** strands. Formula: `Tm = 1000*dH / (A + dS + R*ln(Ct/4)) - 273.15 + 16.6*log10([Na+])`. Deliberately not built on `core.thermo` - all three divergences are named, because the 21.5 C threshold is only meaningful in combination with exactly these choices.
  - Mechanism to explain on the wiki: an siRNA also behaves like a microRNA - the seven nucleotides at guide positions 2-8 are enough to knock down any transcript carrying the complementary 3' UTR site. There are only 16384 possible 7-mers against tens of thousands of 3' UTRs, so seed matches **cannot be avoided by sequence choice alone**; Ui-Tei 2008's way out is that off-target strength tracks seed-duplex stability, so a weakly-pairing seed is largely off-target free even where its matches exist (`design_notes.md:116-126`).
  - Both strands are filtered because the passenger is loaded into RISC often enough to drive its own seed-dependent silencing (`design_notes.md:128-130`).
- **`core.offtarget`** (410 lines) - the transcriptome homology filter, HUSH's open offline equivalent of siDirect's GGGenome step. Two searchers behind one `OffTargetSearcher` protocol: `BruteForceSearcher` (pure-Python Hamming scan, the oracle, suited to small curated panels and the tests) and `BowtieSearcher` (bowtie1 `-v`, scales to the transcriptome, one process per run via `precompute()`). Strand semantics: each candidate contributes two reads - the guide side searches the **sense** core, the passenger side the **antisense** core - forward strand only (`--norc`) because transcripts are single-stranded sense mRNA. The candidate's own source transcript is excluded by accession (version-insensitive), auto-detected from the FASTA header via a RefSeq accession regex and overridable.
- **`core.asymmetry`** (183 lines) - duplex-end thermodynamic asymmetry, measured directly in kcal/mol instead of through the 2004 A/U proxies: `differential = dG(guide 5' end window) - dG(sense 5' end window)`. **Positive** means the guide's 5' end is the looser end, the orientation that loads the intended strand; negative warns the passenger is the thermodynamically preferred guide, which wastes the duplex and puts the passenger's seed in charge of off-target behaviour. Window length is a parameter (default 4 bp) because no canonical value exists - Khvorova used 5-nt windows, Schwarz emphasised the few terminal base pairs, later tools differ again. **No threshold is defined**: the module returns a quantity and its sign.
- **`core.hairpin`** (250 lines) - intramolecular hairpin Tm, supplying Reynolds' internal-structure criterion (Tm <= 20 C). A hairpin is **not** a free duplex: it melts intramolecularly (Groebe & Uhlenbeck showed the transition is concentration-independent, identical profiles from 1 to 10 uM), so `Tm = dH / dS` with no concentration term; and it carries a **loop penalty** depending mostly on loop size. Model: Watson-Crick stem stacking (Xia 1998) with no initiation or terminal-A-U term because the loop penalty plays initiation's role; loop dG37 from Groebe & Uhlenbeck Table 3, composition-weighted and size-interpolated; loop destabilisation treated as **purely entropic** (`dH_loop = 0`) because Table 3 gives dG37 only - "the one unavoidable modelling assumption, and it is the standard one". WC stems only, so internal structure is if anything under-estimated - the conservative direction for a filter that rewards its absence. Loop/stem minimums (3 nt, 3 bp) are steric conventions, not from the source, and are tunable.
- **`core.selection`** (479 lines) - the option-rich filter step. Two groups: the three rule sets (toggled independently, combined by `ALL` or `ANY`) and the independent filters ANDed on top (GC window, seed Tm, transcriptome homology, standalone hairpin, homopolymer run, asymmetry).
  - **The default paradigm and its justification** (`design_notes.md:522-540`): all three rule sets on, combined by **`ANY`** - each classical algorithm treated as an equal vote, because `ALL` would demand three independent 2004 heuristics agree, "which erases almost everything". Seed on because it is the siDirect industry standard. Homopolymer on because poly-runs are disruptive to both siRNA stability and polymerase read-through. GC off because it is counter-productively strict on the GC-rich polymorphic families HUSH targets; asymmetry off because it is a niche, conservative signal not found in comparable tools. Transcriptome off-target on **at the CLI** but off in `SelectionCriteria` itself, because the filter has a mandatory searcher argument with no sensible default.
  - **Options are not rule edits**: the standalone GC filter is separate from Reynolds' internal GC criterion (pinned at the published 30-52%), and the standalone hairpin gate is separate from Reynolds' internal-repeat sub-criterion (pinned to Xia 1998 at Tm <= 20 C). The rule keeps its published form; the standalone filter is the tunable one.
  - **Where the parameter set may not go**: `--params` is scoped to the asymmetry and hairpin readouts only, deliberately not wired into the seed filter, because exposing it there would let a caller silently invalidate the 21.5 C threshold.
- **`core.errors`** (31 lines) - shared exception types (`HushError` root).

#### `msa/` - the allele stage (stage 3)
- **`msa.alignment`** (245 lines) - the result type every backend produces and every consumer reads (`PanelAlignment` / `AlignedRow`). Nothing is aligned here. Three states kept distinct per position: a residue, `GAP` (no residue - "not the same as different"), and an insertion (extra residues with no target coordinate, recorded in `AlignedRow.insertions` rather than dropped). The projection is lossy by construction.
- **`msa.align_star`** (269 lines) - the star backend: every member aligned pairwise against the target, the target as the fixed centre. **Gusfield's center-star construction**, which carries a factor-2 sum-of-pairs bound when the distance obeys the triangle inequality - but here the centre is **fixed by the pipeline** (candidates are enumerated from the target), so the bound's premise does not apply and no optimality is claimed. "What the construction buys us is a consistent frame, not a quality guarantee." Gap costs default to making an indel buy back at least four mismatches before it pays for itself, which keeps the aligner from scattering frame-shifting gaps to chase a few mismatches - the failure mode that would move which positions look polymorphic. Where several alignments score equally the first is taken.
- **`msa.variation`** (228 lines) - per-position match / mismatch / uncovered. States stored as a **string** per member rather than a list of objects, so a whole panel is cheap (an IPD-IMGT/HLA locus can hold thousands of members over a few thousand bases) and a window is a slice.
- **`msa.windows`** (484 lines) - the two gates, in the order a rejection is reported: targeted alleles must be **covered and identical** across the window (the guide is built from the target, so a difference there is a mismatch against the very allele it should silence); every spared allele must be **covered and carry at least `min_spare_mismatches` differences**. Plus region merging. Windows indexed by 0-based start and length - the same coordinates `core.sirna` gives its candidates.
- **`msa.targeting`** (256 lines) - the request model `AlleleTargeting`: two disjoint id sets (`hit`, `spare`) with the four named cases provided as **constructors, not separate types**; admits combinations the named cases do not cover (e.g. sparing only part of the panel and leaving the remainder unconstrained).
- **`msa.discrimination`** (404 lines) - how *well* a single mismatch discriminates. Two levels of evidence used in order: **Huang et al. 2009** profiled every mismatch at every position and published a lookup keyed by the exact base change (Table 1) - the most specific statement available, measured rather than inferred; **Schwarz et al. 2006** established the general principles (purine:purine mismatches discriminate best; mismatches 3' of the seed discriminate robustly while seed mismatches are overwhelmed by high siRNA concentration). Above both sits the **number** of mismatches: both papers address the single-difference case, and **Ohnishi et al. 2008** showed a second mismatch - even inside the seed - markedly improves discrimination. Ranking is a lexicographic tuple of published facts, never a weighted score. A window is ranked by its **weakest** spared allele (the same "weakest link" rule `min_spare_mismatches` follows).
  - Honest caveats recorded in the docstring, all usable as wiki text: the tuple's order is only partly evidenced (Huang before Schwarz is well founded, specific beats general; region before position rests on a difference in kind; position before class and the count threshold are judgements). Huang and Schwarz disagree in places - A>C lists guide position 8, inside the seed - and disagreements are **reported** via `MismatchAssessment.sources_agree` rather than resolved. Schwarz lists P5 among consistently discriminating positions although it is in the seed; the published set is reproduced intact and region checked first, placing P5 in the seed tier. Within purine:purine there is further structure (G:G and A:G more selective than A:A and G:A at P10; A:G alone held up at 20 nM) which is not modelled. **Guide position 1 is held in Argonaute's MID pocket and does not pair with the target, so a "mismatch" there is not one** - reported as `GuideRegion.ANCHOR`. IUPAC codes yield `None` rather than a guess.
- **`msa.msa_errors`** (40 lines) - the stage's exception vocabulary, all deriving from `HushError`.

#### `shrna/` - conversion (stage 4)
- **`shrna.loops`** (82 lines) - the loop registry, each entry carrying its sequence (RNA 5'->3') and citation: `jensen` = **UGUGCUU** (default), `psuper` = **UUCAAGAGA** (Brummelkamp 2002), `trc` = **CUCGAG** (pLKO.1 / RNAi Consortium). Jensen is the default because it was library-selected specifically for robust knockdown with the **least stem-sequence dependence** - "the safest property for an automated tool that does not know the stem in advance".
- **`shrna.convert`** (203 lines) - the conversion and scoring. The 19-bp core becomes the stem **unchanged**, so the guide the design stage validated (seed-Tm, asymmetry, off-target) is exactly the guide the hairpin expresses: no re-shift, no re-validation. Emitted in sense-loop-antisense (SA) orientation as the bare hairpin insert (top and bottom DNA oligos, bottom = RC of top), with a Pol III poly-T terminator and **no vector-specific restriction overhangs**. Constants: `_POL3_TERMINATOR_RUN = "UUUU"` (a run of >= 4 U/T is read by Pol III as termination), `MIN_TERMINATOR_LENGTH = 4`, `DEFAULT_TERMINATOR_LENGTH = 6`. `ShrnaConstruct` fields include `fold`, `folds_as_designed`, `starts_with_g`, `premature_terminator_free`, `max_homopolymer_run`.
  - Literature parameters the scaffold rests on (`shrna_design_note.md:28-58`): canonical form is the **pSUPER** design (Brummelkamp 2002) - a Pol III promoter (U6 or H1) transcribing inverted repeat - loop - inverted repeat, closed by a poly-T terminator (5-6 T); Dicer removes the loop after nuclear export, regenerating an siRNA-like duplex. Stem 19-29 bp: 19 mimics a classical siRNA, 21-29 act as Dicer substrates and shift the guide by a few nt (Ge et al. 2009 is the key structure-activity study; stems <= 19 bp, "sshRNA", have distinct behaviour and can be extremely potent). Loop 4-11 nt, short (4-9) preferred; runs of T/U must be avoided because Pol III reads them as termination. Wang et al. 2013 found 8-nt loops best overall (`TTCTAGAA`, `TTGGCCAA`), a 4-nt loop dropping below 50% knockdown, and TTCAAGAGA beating CTCGAG. Orientation: both SA and AS are used, the optimum is target-dependent with no universal winner. U6 strongly prefers a **G at +1**.
  - The **AgoshRNA boundary**: a very short stem (~18-19 bp) with a small loop (3-5 nt) can be processed by Argonaute-2 instead of Dicer, yielding only the guide strand - eliminating passenger off-target entirely - but following different design rules. HUSH's 19-mer cores sit exactly on this boundary, so the module must choose deliberately (`shrna_design_note.md:60-68`).
  - **[STALE]** `docs/shrna_design_note.md` still carries `Status: scope fixed, references gathered, design decisions pending` and lists five open decisions, but the module is implemented and `README.md:27` says so. The five decisions have in fact been made in code: loop = selectable set with Jensen default; stem = 19-mer core as-is; orientation = SA only (enum-extensible); processing = Dicer; cloning = bare insert plus terminator, vector-agnostic. Worth fixing in the repo before anyone reads it alongside the wiki.

#### `service/` - the orchestration seam
- **`service.params`** (103 lines) - `DesignParams` / `AlleleTargetingSpec`: framework-agnostic, already-parsed requests. GC bounds are **fractions** here, not percent (the CLI/API convert).
- **`service.results`** (63 lines) - `InspectResult` / `DesignResult` / `ConvertResult`: lossless carriers, no formatting decisions.
- **`service.session`** (348 lines) - `HushSession`: one working instance of HUSH - a loaded target, an optional panel, and **the MSA alignment cached** so repeated design runs with different parameters do not re-align.

#### `app/` - the CLI layer
- **`app.cli`** (1215 lines, the largest module) - the one-shot commands and all rendering: `_print_criteria`, `_print_funnel`, the table, `_write_design_tsv`, `render_design_result`, `render_convert_result`. Also the **composition root** where the allele gate is ANDed with the selection verdict.
- **`app.shell`** (515 lines) - the REPL. Presentation only; all pipeline work delegated to the session. Renders through the very same functions as the one-shot CLI, so a design in the shell reads identically to `hush design`. Commands: `load`, `panel`, `set`, `design`, `convert`, `save`, `show`/`inspect`, `assess`, `limit`, `reset`, `help`, `quit`.
- **`app.parsing`** (92 lines) - shared flag-vocabulary parsing, so the one-shot CLI and the shell accept the same option words. Raises plain `ValueError` with a user-facing message; Typer re-wraps it as a usage error with a flag hint, the shell prints it as-is.

#### `api/` - the HTTP layer (optional extra)
- **`api.app`** (150 lines) - two ways in, one code path: stateless endpoints build a throwaway session per request, `/sessions` endpoints keep one alive. `HushError` and `ValueError` map to 400, unknown session id to 404.
- **`api.schemas`** (337 lines) - Pydantic v2 models defining the JSON contract and, through FastAPI, the OpenAPI schema the browser client generates its typed client from.
- **`api.serialize`** (292 lines) - the counterpart to the CLI's printers; kept separate from the schemas so wire shape and mapping evolve independently.
- **`api.store`** (55 lines) - in-memory session store: opaque ids, serialised access (FastAPI runs sync endpoints in a threadpool), capacity-bounded with **LRU eviction** so a long-lived server cannot grow without limit.
- **`api.inputs`** (35 lines) - request text to domain model.
- **`api.__main__`** (26 lines) - `python -m api` launcher, env-configured.

#### `tools/` - developer utilities
- **`tools/build_offtarget_db.py`** - downloads RefSeq, runs `bowtie-build`, writes the reproducibility manifest.
- **`tools/compare_sidirect.py`** - diffs a HUSH design TSV against a siDirect export, including the off-target columns.

#### Endpoints and commands, for a quick reference table
- CLI: `hush inspect`, `hush design`, `hush convert`, `hush shell`.
- API: `POST /design`, `POST /sessions`, `POST /sessions/{id}/design`, `POST /sessions/{id}/convert`, `DELETE /sessions/{id}`, plus `/docs` and `/openapi.json`.

---

## [LATER] Numbers parked for the Validation chapter

Collected so you do not have to re-mine them tomorrow.

- **`core.thermo`**: Tm agrees with Biopython's independent implementation of the same Xia table (`RNA_NN2`) to within **0.007 C**; the residual is entirely Biopython's gas constant rounded to 1.987 against the accurate 1.9872. `tests/test_thermo.py`.
- **`core.seed`**: computing the Tm for all 4^7 = **16384** heptamers reproduces siDirect's own published figures exactly - **4488 (27.4%)** below 21.5 C, spanning **-12.0 C to 59.6 C** (the paper reports "ranged from -12 C to 60 C"). One check that pins down every constant at once. `tests/test_seed.py`.
- **`core.rules` / Amarzguioui**: all **46** siRNAs of the paper's Table 3 reproduce exactly - every DT3, every motif, every total - plus the reported class means (2.6 functional / 0 non-functional) and selection statistic (14 selected, 11 functional). Fixture `tests/data/amarzguioui_table3.tsv`.
- **`core.asymmetry`**: on the real HLA CDS the sign of the thermodynamic differential agrees with Amarzguioui's A/U-count proxy in **95.2%** of candidates where the proxy is decisive.
- **siDirect parity (PDCD1 CDS, RefSeq NM_005018.3)**: with matched functional-rule settings HUSH kept a **perfect superset** - all of siDirect's picks plus **5 extra**. The extras were exactly the candidates siDirect removed for near-perfect off-targets; the decisive one at CDS **650-668** clears all three rules plus GC and homopolymer, so no rule difference could explain the disagreement. Enabling the homology gate drops them to zero.
- **Off-target gate calibration**: default `offtarget_min_mismatch = 2` because siDirect keeps candidates whose nearest off-target is >= 2 mismatches away and drops only the near-perfect (0/1) ones; a threshold of 3 would wrongly drop siDirect's own 2-mismatch picks. Search depth defaults to 3 (the bowtie `-v` cap) so the reported columns line up with siDirect's.
- **The seed filter's over-conservatism, measured** on `data/examples/hla_a01.fasta` (HLA-A*01:01 CDS, **1098 nt, 63.8% GC, 1080 candidates**): Ui-Tei functional = **85**; of those, seed-clean on both strands = **0**; seed-clean on both strands regardless of Ui-Tei = **2**. Mechanism: median guide seed Tm is 37.9 C over all candidates but 28.4 C among Ui-Tei passers (rule (iii) pulls it down, 16 of 85 clear 21.5 C), while median passenger seed Tm stays 38.2 C (rule (ii) requires G/C at the sense 5' terminus, so on a 64% GC gene **0 of 85** clear it). The two filters are in direct tension on GC-rich targets.
- **Freier vs Xia, measured across all 16384 heptamers**: Freier range -12.0 to 59.6 C, Xia -31.3 to 61.1 C; the Xia threshold selecting the same 27.4% is 22.31 C; Pearson r = **0.988**, Spearman rho = **0.990**; under quantile-matched thresholds only **2.3%** of heptamers change verdict. Conclusion: swapping the table is a near-no-op that would nonetheless invalidate the 21.5 C threshold. The shared systematic error is that both tables measure **free RNA duplexes in solution**, whereas the seed duplex forms inside Argonaute, which pre-organizes guide nucleotides g2-g8 in a near-helical conformation - seed pairing proceeds **50-200x faster** than annealing of the equivalent naked RNAs.
- **The anticorrelation siDirect surfaces but does not resolve** (its Figure 4): seeds with **lower** Tm are **more** frequently present in human 3' UTRs. So the filter selects for seeds that bind each site weakly while enriching for seeds that have more sites to bind. Proposed replacement: `expected off-target load ~ (seed-matched sites in EXPRESSED transcripts) x (repression per site)` - siDirect controls only the second factor. Predicted failure mode: a high-Tm seed that is rare in the expressed transcriptome is genuinely clean and wrongly rejected.
- **Test suite**: **548** test functions across 20 test modules (`tests/`), some parametrized so the case count is higher. Notable: a skipped-if-absent parity test asserts `BowtieSearcher` and `BruteForceSearcher` agree on a freshly built index.

## [LATER] Documentation chapter - what exists to point at
- `docs/HUSH_concept_dossier.md` - vision and prior art (85 lines).
- `docs/design_notes.md` - the long-form companion: rationale, formulae, validation records per stage (761 lines). Explicit division of labour: "the modules say *what* they do and the critical caller contract; the *why*, the derivations, and the checked-against-the-paper's-own-figures records live here so the code stays readable."
- `docs/design_options.md` - one-page reference to every option, its default and effect (210 lines), including a "What is deliberately not tunable" section.
- `docs/shrna_design_note.md`, `docs/offtarget_design_note.md`, `docs/seed_offtarget_design_note.md` - per-decision records.
- `docs/offtarget_module_issue.md` - a resolved hand-off spec, kept for provenance; a good example of the development process if the Documentation chapter wants one.
- Module docstrings carry the per-unit contract, scope and deferrals.
- `docker/README.md` - container and index operations (107 lines).

## [GAP] Things the wiki will want that the repo does not state
- A public repository URL, a CITATION file, and a citable DOI (Zenodo).
- Author/contributor list (`pyproject.toml` says only "HUSH contributors").
- The hosted instance's URL (the deployment posture is decided but no address is recorded).
- A changelog or release history; version is still `0.0.1`.
- Stage 5 "Output Layer" is marked planned but never specified anywhere - decide whether to present it as planned or drop it from the table.
- Citations for three things used but not cited: **GGGenome** (siDirect's search engine), **Gusfield's center-star construction**, and **IPD-IMGT/HLA** (the website is linked but the database's own paper is not cited). Dharmacon siDESIGN is also mentioned in passing without a reference.
