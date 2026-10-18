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
Vite proxies the API paths to the API. Development and the container therefore use
identical relative URLs, and no API host is compiled into the bundle.

## Tests

Tests are run with `pytest` over `tests/`, with one test module per domain module.

**The bowtie parity test is skipped without bowtie.** `tests/test_offtarget.py`
compares the bowtie backend with a brute-force search over a small synthetic
transcriptome. It is skipped when `bowtie` and `bowtie-build` are not on PATH, so a
passing local run does not cover that path. CI installs bowtie, and the container
contains it:

```bash
docker run --rm --user root -v "$PWD:/src" -w /src hush:local \
  sh -c 'pip install -q pytest httpx2 && python -m pytest -q'
```

**Use `httpx2`, not `httpx`.** The Starlette `TestClient` requires `httpx2`. With
`httpx` alone, every test in `tests/test_api.py` fails at collection. The `dev`
extra installs the correct package.

**Fixtures** are in `tests/conftest.py`. `target_fasta`, `panel_fasta`,
`hla_a01_fasta` and `demo_lowgc_fasta` point to files in `data/`, and
`amarzguioui_table3` loads the 46-siRNA validation table from `tests/data/`.

**Tests must not depend on environment variables.** `HUSH_WEB_DIR` and
`HUSH_REFERENCE_DIR` change behaviour, and the container sets both. A test that
depends on either must clear it with `monkeypatch.delenv`.

## Conventions

**Selection is a filter, not a score.** A new gate is an independent
`bool | None` field on `CandidateEvaluation`. It is not a number added to a total.
See [design-rules.md](design-rules.md).

**Results hold the complete computed state.** A result contains everything a
renderer needs and has no dependency on terminals or wire formats. If a front end
has to derive a value from a result, add that value to the result. For this reason
withheld shRNA constructs stay on the result and are not dropped.

**Front ends do not call the domain packages directly.** `app` and `api` use
`service.HushSession`. Domain logic duplicated in a renderer previously produced
three copies of the overhang-stripping code.

**Published values are fixed.** The internals of the seed filter and the internal
criteria of Reynolds keep their published form. Making them adjustable would allow
a run to invalidate a published calibration without any indication. To make a
criterion tunable, add a separate gate, as `--gc` and `--hairpin` do.

Docstrings are short. Longer explanations belong in these docs and not in the
source.

## CI

`.github/workflows/ci.yml` runs three jobs on every push and pull request: the
Python tests with bowtie installed, the web typecheck, lint and build, and an image
build that is smoke-tested for `/health` and the served web client. The Dockerfile
is part of the release, so a failed image build blocks distribution.

`.github/workflows/release.yml` publishes the multi-arch image to GHCR on a `v*`
tag.

## The archive branch

The dated design notes and the full literature list are on `archive/design-notes`.
They are kept out of the release tree because they document provenance and are not
user documentation:

```bash
git show archive/design-notes:docs/design_notes.md
git show archive/design-notes:docs/wiki/references.md
git show archive/design-notes:docs/msa_ranking_note.md
```

Committed design-run outputs from before the overhang correction are also on that
branch, under `results/`. They predate the geometry change and are not current
output.
