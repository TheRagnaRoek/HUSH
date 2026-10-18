# HTTP API

The pipeline as a JSON API (FastAPI). Interactive documentation is at `/docs` and
the OpenAPI schema at `/openapi.json`. A TypeScript client can be generated from
the schema with `npm run gen:api` in `web/`.

```bash
python -m api                      # or: uvicorn api.app:app --reload
```

The container serves the API and the web client together on port 8000.

## Endpoints

| Method | Path | Purpose |
|---|---|---|
| `GET` | `/health` | Liveness check; returns service name and version, and whether bowtie and the off-target reference are available |
| `POST` | `/inspect` | Preview a target and optional panel |
| `POST` | `/design` | One-shot design; the request contains all inputs |
| `POST` | `/sessions` | Create a session from a target and optional panel |
| `GET` | `/sessions/{id}` | Session state, including whether a design exists |
| `POST` | `/sessions/{id}/design` | Design with the cached inputs of the session |
| `POST` | `/sessions/{id}/convert` | Convert the last design of the session |
| `DELETE` | `/sessions/{id}` | Delete a session |

Target and panel are sent as FASTA **text** in the JSON body, not as file uploads.

## Stateless use

A single request carries all inputs:

```bash
curl -s localhost:8000/design -H 'content-type: application/json' -d '{
  "target_fasta": ">t\nACGT...",
  "params": {"use_seed_filter": false}
}'
```

## Sessions

Post the inputs once, then design repeatedly. The server caches the panel
alignment, so re-designing with new parameters does not repeat the alignment:

```bash
SID=$(curl -s localhost:8000/sessions -H 'content-type: application/json' \
  -d '{"target_fasta": ">t\nACGT...", "panel_fasta": ">a\nACGT..."}' | jq -r .id)

curl -s localhost:8000/sessions/$SID/design -H 'content-type: application/json' \
  -d '{"use_seed_filter": false, "sort": "discrimination"}'

curl -s localhost:8000/sessions/$SID/convert -H 'content-type: application/json' \
  -d '{"loop": "jensen"}'
```

Sessions are held in a size-limited in-memory cache, not in a database. They are
lost on restart, and the least recently used session is evicted when the limit is
reached.

## Parameters

`params` has the same fields as the CLI options, for example `use_seed_filter`,
`gc_min`, `offtarget_min_mismatch`, `targeting` and `sort`. See [cli.md](cli.md)
for every field and its default.

The API differs from the CLI in two fields: **`bowtie_path` and
`offtarget_self_fasta` are not accepted over HTTP.** The API has no
authentication, so it must not accept an executable path or a path of a file to
read. The CLI keeps both options. Resolution of the design's own transcripts still
runs in the API, because it reads the uploaded target and panel and not the file
system. `offtarget_self_accessions`, `resolve_offtarget_self` and
`offtarget_self_coverage` are accepted.

## Responses

A design response contains the criteria that were applied, the funnel (the number
of candidates that passed each gate), and the selected candidates with their
per-criterion breakdown. This is the same data that the CLI renders, so a client
does not need to recompute a verdict.

`offtarget_skipped` is true when the run had no reference. In that case an empty
off-target hit count means "not checked", not "no off-targets found". See
[design-rules.md](design-rules.md#transcriptome-off-target).

The criteria summary contains `offtarget_self_accessions` (the transcripts that the
gate counted as the design's own). `offtarget_self_declared` and
`offtarget_self_resolved` state how many were named and how many were found by
homology, and `offtarget_self_coverage` is the threshold that was used. An empty
list with an active filter means that every hit counted against every candidate.
For a gene family this indicates a faulty run, not a clean result; see
[What the design owns](design-rules.md#what-the-design-owns).

A convert response separates the constructs that passed from the constructs that
were withheld:

| Field | Meaning |
|---|---|
| `constructs` | Constructs that passed the conversion filters |
| `withheld` | Withheld constructs, with their QC flags |
| `construct_count` / `withheld_count` | Their counts |
| `converted_count` | The total built, before filtering |
| `stem_overhang` | The sense-overhang mode used, `keep` or `trim` (request field of the same name; default `keep`) |

`construct_count + withheld_count == converted_count` always holds, so every
converted construct is accounted for.

## Errors

Domain errors (`HushError`) and invalid requests (`ValueError`) return status
**400** with a `{"detail": "..."}` body. An unknown session id returns **404**.

## Security

The API has no authentication, no rate limiting and no TLS. It is intended to run
on the local machine. Do not expose the port beyond localhost without adding
authentication and a reverse proxy.

CORS is permissive only when the API does not serve a web client. When it serves
one, all calls are same-origin and the CORS middleware is not installed.
`create_app(allow_origins=[...])` overrides this.
