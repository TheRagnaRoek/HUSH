# The iGEM GitLab mirror

Competition rules require the submitted software to live on the team's
repository on [iGEM's GitLab](https://gitlab.igem.org/); code hosted elsewhere
is not evaluated. Development happens on GitHub, so the GitLab repository is a
mirror, frozen at the submitted version.

This page records how to produce that mirror and exactly how it differs, so the
difference is a deliberate, repeatable edit rather than two copies drifting
apart.

## The frozen copy must stand alone

A judge clones the GitLab repository and runs it. Nothing outside that clone is
guaranteed to exist at that moment, so the documented path must not depend on
one.

That is why the quickstart leads with:

```bash
docker compose run --rm reference-build
docker compose up hush
```

Compose builds from `docker/Dockerfile` when the image is not already present,
so this path touches no registry at all. It is also what the rules ask for:
the source must produce the output, rather than a prebuilt binary being handed
over.

The published image is offered second and is **pinned to an exact version**,
never `latest`. A moving tag would mean the frozen repository's own documented
command pulls whatever was published most recently - so the code a judge reads
and the software a judge runs could silently differ. Every reference in the
repository names `ghcr.io/theragnaroek/hush:1.0.0`.

The one thing deliberately not shipped is the transcriptome reference index.
It is multi-gigabyte, it is built from public RefSeq data by
`hush reference build`, and the index manifest records exactly which RefSeq
release went into it. That is reproducible, not missing.

## Producing the mirror

```bash
git remote add igem https://gitlab.igem.org/2026/software-tools/<team>.git
git push igem main
git push igem v1.0.0
```

Push the tag as well: it is what ties the mirror to the published image and to
the version in `pyproject.toml` and `CITATION.cff`.

`.gitlab-ci.yml` runs there in place of `.github/workflows/`. The GitHub
workflows are harmless on GitLab - they are simply ignored - and keeping them
avoids a difference that would have to be re-applied on every sync.

## The differences to apply

Only the README changes. Add this directly under the logo, above the `# HUSH`
heading:

```markdown
> **iGEM 2026 competition submission.** This is the frozen copy of HUSH as
> submitted for judging: version 1.0.0, tag `v1.0.0`. Development continues at
> <https://github.com/TheRagnaRoek/HUSH>; later versions there will differ from
> what is recorded here.
```

And in `CITATION.cff`, add the mirror alongside the development repository:

```yaml
repository: https://gitlab.igem.org/2026/software-tools/<team>
```

That is the whole diff. Everything else is identical by design.

## Checking the mirror before you call it done

Do this from a clone of the **GitLab** repository, on a machine that has never
built HUSH, with no registry login:

```bash
docker compose up hush          # must build from source and serve :8000
pytest -q                       # after: pip install -r requirements.txt -r requirements-dev.txt
python tools/compare_sidirect.py \
  --sidirect data/testing/siDirect_PDCD1-CDS.ods \
  --hush data/testing/pdcd1_hush_relaxed.tsv
```

If the first command reaches <http://localhost:8000> without pulling anything,
the frozen copy is self-sufficient.

Also confirm, outside the clone:

- the repository is under 50 MB (`git count-objects -vH`),
- `LICENSE` is present and the license is OSI-approved (MIT),
- `CITATION.cff` names the real author list rather than a placeholder,
- the published image, if you reference it, is actually public.

## If you have to fix something after submitting

Fix it on GitHub, tag a new version, and push that tag to the mirror as well.
Do not edit the mirror directly: a change that exists only there is lost the
next time anyone syncs, and the two copies stop being the same software.
