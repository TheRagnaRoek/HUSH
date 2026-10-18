# HUSH documentation

Start with the [README](../README.md) for installation and a first run. These pages
give the details.

| Page | Read it to |
|---|---|
| [architecture.md](architecture.md) | Understand how the stages and layers fit together, or decide where a change belongs |
| [design-rules.md](design-rules.md) | Find out what each gate tests, why it exists and what it is based on |
| [shrna-conversion.md](shrna-conversion.md) | Understand the emitted hairpin: geometry, loops, QC read-outs |
| [cli.md](cli.md) | Look up a command or a flag |
| [api.md](api.md) | Call HUSH over HTTP |
| [development.md](development.md) | Run the tests, or work on the code |
| [../docker/README.md](../docker/README.md) | Run or build the container |
| [igem-mirror.md](igem-mirror.md) | Produce the frozen iGEM GitLab copy and see why it is pinned |

## Key principles

**Selection is a filter, not a score.** Every criterion is an independent gate, and
a candidate is selected only if it clears every enabled gate. Nothing is combined
into a single number, so every verdict can be traced to the gates that produced it.
The only ordering HUSH applies is the discrimination ranking in allele-aware runs.
It reorders the selected candidates and never admits or rejects one.

**The transcriptome off-target check is on by default.** It locates the reference
built by `hush reference build` automatically. A run without the built reference
still works, for example to use a different reference genome, but every output
labels it as unreferenced. See
[design-rules.md](design-rules.md#transcriptome-off-target).

## Source material

The literature behind these decisions, and the dated design notes that record how
each decision was made, are on the `archive/design-notes` branch. They are kept out
of the release tree because they document provenance and are not user
documentation:

```bash
git show archive/design-notes:docs/design_notes.md
git show archive/design-notes:docs/wiki/references.md
```
