# Development

```bash
python -m venv .venv
# Windows:  .venv\Scripts\activate
# POSIX:    source .venv/bin/activate
pip install -e ".[dev,api]"
pytest
```

For the web client:

```bash
cd web
npm ci
npm run dev        # :5173, proxies the API paths to :8000
npm run typecheck
npm run lint
npm run build
```

Run `python -m api` alongside `npm run dev`. The client calls its own origin and
Vite proxies the API paths across, so dev and the container use identical
relative URLs and no API host is baked into the bundle.

## Tests

`pytest` over `tests/`, one module per domain module. A few things to know:

**The bowtie parity test skips without bowtie.** `tests/test_offtarget.py`
checks the bowtie backend against a brute-force search over a tiny synthetic
transcriptome. It skips when `bowtie` and `bowtie-build` are not on PATH, so a
green local run does not mean that path was exercised. CI installs bowtie, and
the container has it:

```bash
docker run --rm --user root -v "$PWD:/src" -w /src hush:local \
  sh -c 'pip install -q pytest httpx2 && python -m pytest -q'
```

**`httpx2`, not `httpx`.** Starlette's `TestClient` requires it; with `httpx`
alone every `tests/test_api.py` test fails at collection. The `dev` extra pins
the right one.

**Fixtures** are in `tests/conftest.py`: `target_fasta`, `panel_fasta`,
`hla_a01_fasta`, `demo_lowgc_fasta` point at `data/`, and `amarzguioui_table3`
loads the 46-siRNA validation table from `tests/data/`.

**Tests must not read ambient environment.** `HUSH_WEB_DIR` and
`HUSH_REFERENCE_DIR` change behaviour, and the container sets both; a test that
depends on either should clear it with `monkeypatch.delenv`.

## Conventions

**Selection is a filter, not a score.** Adding a gate means adding an
independent `bool | None` to `CandidateEvaluation`, not folding a number into a
total. See [design-rules.md](design-rules.md).

**Results are lossless carriers.** A result holds the computed state a renderer
needs and knows nothing about terminals or wire formats. If a front end has to
derive something from a result, that something belongs on the result. This is
why withheld shRNA constructs ride along rather than being dropped.

**Front ends do not reach into the domain.** `app` and `api` go through
`service.HushSession`. Duplicating domain logic in a renderer is how the same
overhang-stripping ended up in three places.

**Pinned values stay pinned.** The seed filter's internals and Reynolds'
internal criteria keep their published form; exposing them would let a run
silently invalidate a published calibration. Add a standalone tunable gate
instead, as `--gc` and `--hairpin` do.

Docstrings are short. Longer explanation belongs in these docs, not in the
source.

## CI

`.github/workflows/ci.yml` runs three jobs on every push and pull request:
Python tests with bowtie installed, the web typecheck/lint/build, and an image
build smoke-tested for `/health` and the served SPA. The Dockerfile is a shipped
artefact, so a build failure is a distribution failure.

`.github/workflows/release.yml` publishes the multi-arch image to GHCR on a
`v*` tag.

## The archive branch

The dated design notes and the full literature list are on
`archive/design-notes`, kept out of the release tree as provenance rather than
user documentation:

```bash
git show archive/design-notes:docs/design_notes.md
git show archive/design-notes:docs/wiki/references.md
git show archive/design-notes:docs/msa_ranking_note.md
```

Committed design-run outputs from before the overhang correction are there too,
under `results/`. They predate the geometry change and should not be read as
current output.
