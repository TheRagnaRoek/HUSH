# Architecture

HUSH is a five-stage pipeline with a single orchestration layer (`service`). Three
front ends (a CLI, an HTTP API and a browser client) call this layer; none of them
contains pipeline logic.

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

| Layer | Package | Contents |
|---|---|---|
| Domain | `core`, `msa`, `shrna` | Sequences and biology. No I/O, no formatting, no CLI. |
| Orchestration | `service` | Running a pipeline. Contains `HushSession` and the parameter and result objects. |
| Presentation | `app`, `api`, `web` | Terminals, HTTP and browsers. No pipeline logic. |

Dependencies point one way: `app` and `api` depend on `service`, `service` depends
on the domain packages, and no package depends on a package above it. A front end
must not call `core` directly. Doing so previously duplicated the overhang-stripping
logic in three places.

## The orchestration layer

`service.session.HushSession` assembles a design run in one place. It holds the
loaded target and panel, caches the panel alignment (computed once and reused
across re-designs), and keeps the last design so that `convert` can use it without
writing a file.

```python
session = HushSession()
session.load_target("target.fasta")
session.load_panel("panel.fasta")          # optional: engages stage 3
result = session.design(DesignParams(...))  # DesignResult
constructs = session.convert(ShrnaOptions(...))  # ConvertResult
```

`service.params.DesignParams` is the full request. `service.results` contains the
result objects (`DesignResult`, `ConvertResult`, `InspectResult`). A result holds
the complete computed state and no formatting: it does not depend on terminals,
colours or wire formats, and a renderer decides how to display it. For this reason
`DesignResult` keeps the per-criterion breakdown of every candidate in addition to
the selected subset, and `ConvertResult` keeps withheld constructs instead of
returning a shorter list.

## Stage 2: selection

`core.selection.evaluate` runs each enabled gate independently and combines the
results with AND. There is no scoring function and no weighting; see
[design-rules.md](design-rules.md). Each gate contributes a `bool | None`, where
`None` means "not enabled". `CandidateEvaluation` keeps all of them, so every
rejection can be traced to a gate.

## Stage 3: allele targeting

Supplying a panel enables `msa`. The panel is star-aligned to the target
(`msa.align_star`), per-position variation is computed (`msa.variation`), and each
window of candidate length is evaluated (`msa.windows`) for whether it hits the
alleles to silence and discriminates against the alleles to spare. Windows and
candidates use the same coordinates as the target, so a candidate is matched to its
window by start position.

The allele gate is combined with every stage-2 gate by AND. Ranking by
discrimination (`msa.discrimination`) orders the survivors; it never admits or
rejects a candidate.

## Stage 4: conversion

`shrna.convert` grafts a candidate onto a stem-loop scaffold. It takes the
paired-core strands and the two overhangs separately, because it needs the boundary
between core and overhang to compute the QC read-outs. See
[shrna-conversion.md](shrna-conversion.md).

## Where a change belongs

- A new **gate**: `core/` (or `msa/` if it needs the panel), then a field on
  `SelectionCriteria` and `DesignParams`, then the three front ends.
- A new **rendering**: `app/cli.py` for the terminal, `api/serialize.py` for JSON,
  `web/` for the GUI. A renderer must not recompute values; add them to the result.
- A new **pipeline step**: a package at the domain layer, wired in `HushSession`.
- If a front end has to derive a value from a result, add that value to the result.

## Configuration

HUSH has no configuration files. Settings are CLI flags and function arguments. A
few environment variables exist for deployment:

| Variable | Default | Used by |
|---|---|---|
| `HUSH_REFERENCE_DIR` | `data/reference` | `core.reference`, to locate the off-target index |
| `HUSH_WEB_DIR` | unset | `api.app`, to serve a built web client |
| `HUSH_API_HOST` / `HUSH_API_PORT` | `127.0.0.1` / `8000` | `api.__main__` |

If `HUSH_WEB_DIR` is unset, the API serves no web client. This prevents a stale
`web/dist` in a checkout from being served and from removing the CORS settings that
the Vite dev server needs.
