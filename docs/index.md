# HUSH documentation

Start with the [README](../README.md) for installation and a first run. These
pages go deeper.

| Page | Read it when you want to |
|---|---|
| [architecture.md](architecture.md) | Understand how the stages and layers fit together, or decide where a change belongs |
| [design-rules.md](design-rules.md) | Know what each gate actually tests, why it exists, and what it is based on |
| [shrna-conversion.md](shrna-conversion.md) | Understand the hairpin HUSH emits: geometry, loops, QC read-outs |
| [cli.md](cli.md) | Look up a command or a flag |
| [api.md](api.md) | Call HUSH over HTTP |
| [development.md](development.md) | Run the tests, or work on the code |
| [../docker/README.md](../docker/README.md) | Run or build the container |

## Two things worth knowing up front

**Selection is a filter, not a score.** Every criterion is an independent gate,
and a candidate is selected only if it clears every enabled one. Nothing is
combined into a single number, so any verdict decomposes into the gates that
produced it. The one ordering HUSH applies - discrimination ranking for
allele-aware runs - reorders survivors but never admits or rejects one.

**The transcriptome off-target check is the safety floor.** It is on by default
and finds its own reference. A run without one still works, because you may want
to supply a different reference, but it is labelled as unreferenced everywhere
it surfaces. See [design-rules.md](design-rules.md#transcriptome-off-target).

## Source material

The literature behind these decisions, and the dated design notes recording how
each was made, live on the `archive/design-notes` branch. They are kept out of
the release tree because they are provenance, not user documentation:

```bash
git show archive/design-notes:docs/design_notes.md
git show archive/design-notes:docs/wiki/references.md
```
