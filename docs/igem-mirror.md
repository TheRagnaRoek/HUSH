# The iGEM GitLab mirror

The competition rules require the submitted software to be in the team repository
on [iGEM's GitLab](https://gitlab.igem.org/); code hosted elsewhere is not
evaluated. Development takes place on GitHub, so the GitLab repository is a mirror,
frozen at the submitted version.

This page describes how to produce the mirror and exactly how it differs from the
GitHub repository, so that the difference is a repeatable edit and the two copies
do not drift apart.

## The frozen copy must be self-contained

A judge clones the GitLab repository and runs it. Nothing outside that clone is
guaranteed to exist at that time, so the documented path must not depend on
anything outside it.

For this reason the quickstart starts with:

```bash
docker compose run --rm reference-build
docker compose up hush
```

Compose builds from `docker/Dockerfile` when the image is not already present, so
this path does not use any registry. The rules also require this: the source must
produce the output, and a prebuilt binary is not accepted.

The published image is offered second and is **pinned to an exact version**, never
`latest`. A moving tag would make the documented command pull the most recently
published image, so the code that is read and the software that runs could differ
without notice. Every reference in the repository names
`ghcr.io/theragnaroek/hush:1.0.0`.

The transcriptome reference index is intentionally not included. It is several
gigabytes, it is built from public RefSeq data by `hush reference build`, and its
manifest records the RefSeq release that was used. The index is reproducible and
its absence is not an omission.

## Producing the mirror

Build the mirror from the release tag, not from `main`. `main` keeps moving, and
the mirror has to stay frozen at the submitted version. Apply the two edits under
"The differences to apply" on a branch cut from the tag, merge in what iGEM has
already put in the project, then push that branch as the mirror's `main`:

```bash
git checkout -b gitlab-mirror v1.0.0
# apply the two edits below, then:
git commit -am "Mark this copy as the iGEM submission"
git remote add igem https://gitlab.igem.org/2026/software/gu-frankfurt/hush.git
git fetch igem
git merge --allow-unrelated-histories igem/main
# README.md and LICENSE conflict; keep this project's versions of both:
git checkout --ours -- README.md LICENSE
git add README.md LICENSE
git commit --no-edit
git push igem gitlab-mirror:main
git push igem v1.0.0
```

The iGEM project is not empty. It is created from a team template with its own
history: a README, an Apache-2.0 `LICENSE`, `CLAUDE.md` and
`.claude/RESPONSIBLE_AI_USE.md`. A plain push is rejected, and a forced push would
delete iGEM's files and history, so the template is merged in instead. Its two
AI-use files are kept unchanged. This project's README and MIT `LICENSE` replace
the template's, which iGEM allows for any OSI-approved license.

Push the tag as well. It connects the mirror to the published image and to the
version in `pyproject.toml` and `CITATION.cff`. iGEM creates the judging release
from `main` by itself at the Wiki Freeze, so push before then.

`.gitlab-ci.yml` runs on GitLab in place of `.github/workflows/`. GitLab ignores
the GitHub workflows. Keeping them avoids a difference that would have to be
re-applied at every sync.

## The differences to apply

Two files change. In `README.md`, add this directly under the logo, above the
`# HUSH` heading:

```markdown
> **iGEM 2026 competition submission.** This is the frozen copy of HUSH as
> submitted for judging: version 1.0.0, tag `v1.0.0`. Development continues at
> <https://github.com/TheRagnaRoek/HUSH>; later versions there will differ from
> what is recorded here.
```

In `CITATION.cff`, add the mirror next to the development repository:

```yaml
repository: https://gitlab.igem.org/2026/software/gu-frankfurt/hush
```

The merge with the template adds iGEM's `CLAUDE.md` and
`.claude/RESPONSIBLE_AI_USE.md` unchanged, so the mirror differs from the tag in
four files. No other file differs.

## Checking the mirror

Check from a clone of the **GitLab** repository, on a machine that has never built
HUSH and has no registry login:

```bash
docker compose up hush          # must build from source and serve :8000
pytest -q                       # after: pip install -r requirements.txt -r requirements-dev.txt
```

If the first command reaches <http://localhost:8000> without pulling any image, the
frozen copy is self-contained.

Also check, outside the clone, that:

- the repository is under 50 MB (`git count-objects -vH`),
- `git diff --stat v1.0.0 HEAD` lists only the four files named above,
- `LICENSE` is present and the license is OSI-approved (MIT),
- `CITATION.cff` contains the real author list and not a placeholder,
- the published image, if referenced, is public.

## Fixing something after submission

Fix it on GitHub, tag a new version, and repeat the steps above from the new tag:
cut the branch from it, apply the same two edits (with the new version and tag in
the banner), merge `igem/main` and push the branch and the tag. The mirror already
contains this project's history, so that merge needs no `--allow-unrelated-histories`
and has no conflicts. Do not edit the mirror directly: a change that exists only in
the mirror is lost at the next sync, and the two copies are no longer the same
software.
