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
without notice. Every reference in the repository names the image of the current
version, for example `ghcr.io/theragnaroek/hush:1.1.0`.

The transcriptome reference index is intentionally not included. It is several
gigabytes, it is built from public RefSeq data by `hush reference build`, and its
manifest records the RefSeq release that was used. The index is reproducible and
its absence is not an omission.

## Producing the mirror

Build the mirror from the release tag, not from `main`. `main` keeps moving, and
the mirror has to stay frozen at the submitted version. Apply the two edits under
"The differences to apply" on a branch cut from the tag, merge in what iGEM has
already put in the project, then push that branch as the mirror's `main`. The
commands below are the ones used to create the mirror at 1.0.0; for a later
version see "Fixing something after submission":

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
> submitted for judging: version 1.1.0, tag `v1.1.0`. Development continues at
> <https://github.com/TheRagnaRoek/HUSH>; later versions there will differ from
> what is recorded here.
```

In `CITATION.cff`, add the mirror next to the development repository:

```yaml
repository: https://gitlab.igem.org/2026/software/gu-frankfurt/hush
```

The merge with the template adds iGEM's `CLAUDE.md` and
`.claude/RESPONSIBLE_AI_USE.md` unchanged, so the mirror differs from the tag in
four files. For 1.0.0 it also differed in `LICENSE`: the copyright holder was
changed to the team name on GitLab before the software freeze. The same change is
on GitHub `main`, so every later release carries it and `LICENSE` no longer
differs. No other file differs.

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
- `git diff --stat v<version> HEAD`, for the version the mirror is at, lists only
  the files named above (four, or five for 1.0.0 with `LICENSE`) and any
  documentation files merged in as described under "Fixing something after
  submission",
- `LICENSE` is present and the license is OSI-approved (MIT),
- `CITATION.cff` contains the real author list and not a placeholder,
- the published image, if referenced, is public.

## Fixing something after submission

Fix it on GitHub, release it there (version bump in `pyproject.toml`,
`CITATION.cff` and `web/package.json`, the image pins, the changelog, then the
tag on the merged `main`), and bring the new tag to the mirror. Do not edit the
mirror directly: a change that exists only in the mirror is lost at the next
sync, and the two copies are no longer the same software.

For 1.1.0 this was:

```bash
git fetch origin --tags
git fetch igem
git checkout -b gitlab-mirror-1.1.0 v1.1.0
# apply the two edits above, with the new version and tag in the banner, then:
git commit -am "Mark this copy as the iGEM submission"
git merge igem/main
```

The mirror already contains this project's history, so the merge needs no
`--allow-unrelated-histories`. The one conflict to expect is the version line of
the banner in `README.md`, because `igem/main` still carries the old banner: keep
the new version. `CITATION.cff` already has the `repository:` line and merges
cleanly. Before pushing, check that `git diff --stat v1.1.0 HEAD` lists only the
four files named above, and that `igem/main` is an ancestor of the branch
(`git merge-base --is-ancestor igem/main HEAD`), so that the push is a
fast-forward and needs no force. Then push the branch and the tag:

```bash
git push igem gitlab-mirror-1.1.0:main
git push igem v1.1.0
```

Documentation-only changes need no new version. For 1.0.0 the README and the docs
were improved on GitHub after the tag and before the software freeze, and reached
the mirror by merging GitHub `main` into the mirror's `main`, not from a tag:

```bash
git fetch origin
git fetch igem
git checkout -b gitlab-docs-sync igem/main
git merge origin/main
```

Before pushing, check that `git diff --stat igem/main HEAD` lists only
documentation files. Then push the branch as the mirror's `main`:

```bash
git push igem gitlab-docs-sync:main
```

The banner still names the released version, which stays true because the
software is unchanged. `git diff --stat v<version> HEAD` then also lists the
documentation files that were merged in. A change to code, a lockfile or CI is not documentation-only
and needs a new version, as described above.
