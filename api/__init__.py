"""HUSH HTTP API (FastAPI).

A thin JSON frontend over the service layer, meant for the browser client. Like
the CLI, it drives the pipeline only through service.HushSession and never
reaches into core/msa/shrna directly; unlike the CLI it renders results as JSON
(api.schemas) rather than terminal tables.

Two ways to call it, sharing one code path (see api.app):

* stateless - POST the target (and panel) with each request; nothing is kept.
* session   - POST the inputs once to get an id, then design/convert by id; the
              server caches the MSA alignment so repeat designs do not re-align.

This package is an optional extra ('pip install -e ".[api]"'); the core library,
CLI and tests do not depend on it.
"""
