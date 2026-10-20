# Changelog

All notable changes to HUSH are recorded here. Versions follow semantic
versioning. The version is declared in `pyproject.toml` and must match
`CITATION.cff` and `web/package.json` (checked by `tests/test_version.py`).

## Unreleased

- The web GUI can export all converted shRNA constructs at once (TSV or TXT,
  in the column layout of `hush convert -o`), from the construct step and the
  candidate list. Withheld constructs are left out, as in the CLI.
  A second export adds the siRNA candidate metrics (the `hush design` columns)
  in front of the construct columns, one row per construct.

## 1.0.0 - 2026-10-20

First release. Includes the allele-aware siRNA/shRNA design pipeline, the CLI,
the HTTP API and web GUI, the optional MCP server, and the Docker image with
bowtie1.

- `hush --version` and the API `/health` report the installed version.
- Docker Compose and the `docker/hush.*` launchers build and run `hush:local`
  from the checkout by default; set `HUSH_IMAGE` to run a published image.
- Release workflow refuses a `v*` tag that differs from the `pyproject.toml`
  version.
- `requirements-mcp.txt` is a lock on top of `requirements.txt`, and
  `pyproject.toml` has an `mcp` extra. The MCP server is a script
  (`mcp/mcp_server.py`) and is not installed as a package.
