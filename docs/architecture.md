# Architecture

HUSH is a five-stage pipeline with one orchestration seam. Three front ends (a
CLI, an HTTP API, a browser client) drive that seam; none of them contains
pipeline logic.

```
            core.inputs / core.sequence          stage 1  input
                        |
            core.sirna -> core.selection         stage 2  siRNA design
              (rules, thermo, seed, offtarget,
               asymmetry, hairpin)
                        |
            msa.*  (optional, with a panel)      stage 3  allele targeting
                        |
            shrna.convert                        stage 4  shRNA conversion
                        |
            app.cli / api.serialize / web        stage 5  output
```

## Layers

| Layer | Package | Knows about |
|---|---|---|
| Domain | `core`, `msa`, `shrna` | Sequences and biology. No I/O, no formatting, no CLI. |
| Orchestration | `service` | How to run a pipeline. Owns `HushSession` and the parameter and result objects. |
| Presentation | `app`, `api`, `web` | Terminals, HTTP and browsers. No pipeline logic. |

The rule is one-directional: `app` and `api` both depend on `service`, `service`
depends on the domain packages, and nothing depends back up. A front end that
reaches into `core` directly is a bug - that is what put the same overhang-
stripping logic in three places before it was consolidated.

## The orchestration seam

`service.session.HushSession` is the single place a design run is assembled. It
holds the loaded target and panel, caches the panel alignment (computed once,
reused across re-designs), and remembers the last design so `convert` can work
from it without a file round-trip.

```python
session = HushSession()
session.load_target("target.fasta")
session.load_panel("panel.fasta")          # optional: engages stage 3
result = session.design(DesignParams(...))  # DesignResult
constructs = session.convert(ShrnaOptions(...))  # ConvertResult
```

`service.params.DesignParams` is the full request; `service.results` holds the
lossless carriers (`DesignResult`, `ConvertResult`, `InspectResult`). A result
never knows about terminals, colours or wire formats - it carries the computed
state and a renderer decides how to show it. That is why `DesignResult` keeps
every candidate's per-criterion breakdown alongside the selected subset, and why
`ConvertResult` keeps withheld constructs rather than returning a shorter list.

## Stage 2: selection

`core.selection.evaluate` runs each enabled gate independently and ANDs the
results. There is no scoring function and no weighting; see
[design-rules.md](design-rules.md). Each gate contributes a `bool | None` -
`None` meaning "not enabled" - and `CandidateEvaluation` keeps all of them, so a
rejection is always attributable.

## Stage 3: allele targeting

Supplying a panel engages `msa`. The panel is star-aligned to the target
(`msa.align_star`), per-position variation is computed (`msa.variation`), and
each candidate-length window is evaluated (`msa.windows`) for whether it hits
the alleles to silence and discriminates against those to spare. Windows and
candidates share the target's coordinates one-for-one, so a candidate is looked
up by start position.

The allele gate is ANDed with every stage-2 gate. Ranking by discrimination
(`msa.discrimination`) orders survivors; it never admits or rejects one.

## Stage 4: conversion

`shrna.convert` grafts a candidate onto a stem-loop scaffold. It takes the
paired-core strands and the two overhangs separately, because it needs the
core/overhang boundary to compute its QC read-outs. See
[shrna-conversion.md](shrna-conversion.md).

## Where a change belongs

- A new **gate**: `core/` (or `msa/` if it needs the panel), then a field on
  `SelectionCriteria` and `DesignParams`, then the three front ends.
- A new **rendering**: `app/cli.py` for the terminal, `api/serialize.py` for
  JSON, `web/` for the GUI. Never recompute in a renderer - add it to the result.
- A new **pipeline step**: a package at the domain layer, wired in `HushSession`.
- Anything a front end needs to *derive* from a result is a sign the result
  should carry it instead.

## Configuration

There is no configuration framework. Everything is CLI flags and function
arguments. A handful of environment variables exist for deployment only:

| Variable | Default | Used by |
|---|---|---|
| `HUSH_REFERENCE_DIR` | `data/reference` | `core.reference`, to locate the off-target index |
| `HUSH_WEB_DIR` | unset | `api.app`, to serve a built web client |
| `HUSH_API_HOST` / `HUSH_API_PORT` | `127.0.0.1` / `8000` | `api.__main__` |

`HUSH_WEB_DIR` being unset means the API serves no client. That is deliberate: a
stale `web/dist` in a checkout must not start being served, nor silently drop
the CORS the Vite dev server needs.
