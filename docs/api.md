# HTTP API

The same pipeline as JSON, over FastAPI. Interactive docs at `/docs`, OpenAPI
schema at `/openapi.json` - a TypeScript client can be generated from it with
`npm run gen:api` in `web/`.

```bash
python -m api                      # or: uvicorn api.app:app --reload
```

The container serves this and the web client together on port 8000.

## Endpoints

| Method | Path | Does |
|---|---|---|
| `GET` | `/health` | Liveness, service name and version |
| `POST` | `/inspect` | Preview a target and optional panel |
| `POST` | `/design` | One-shot design; the request carries everything |
| `POST` | `/sessions` | Create a session from a target and optional panel |
| `GET` | `/sessions/{id}` | Session state, including whether a design exists |
| `POST` | `/sessions/{id}/design` | Design against the session's cached inputs |
| `POST` | `/sessions/{id}/convert` | Convert the session's last design |
| `DELETE` | `/sessions/{id}` | Drop a session |

Target and panel are sent as FASTA **text** in the JSON body, not as uploads.

## Stateless

Everything in one request:

```bash
curl -s localhost:8000/design -H 'content-type: application/json' -d '{
  "target_fasta": ">t\nACGT...",
  "params": {"use_seed_filter": false}
}'
```

## Sessions

Post the inputs once, then design repeatedly. The server caches the panel
alignment, so re-designing with new parameters does not re-align:

```bash
SID=$(curl -s localhost:8000/sessions -H 'content-type: application/json' \
  -d '{"target_fasta": ">t\nACGT...", "panel_fasta": ">a\nACGT..."}' | jq -r .id)

curl -s localhost:8000/sessions/$SID/design -H 'content-type: application/json' \
  -d '{"use_seed_filter": false, "sort": "discrimination"}'

curl -s localhost:8000/sessions/$SID/convert -H 'content-type: application/json' \
  -d '{"loop": "jensen"}'
```

Sessions are a bounded in-memory warm cache, not a database. They are lost on
restart and evicted under pressure.

## Parameters

`params` mirrors the CLI field-for-field - `use_seed_filter`, `gc_min`,
`offtarget_min_mismatch`, `targeting`, `sort` and so on. See
[cli.md](cli.md) for every field and its default.

Two deliberate differences: **`bowtie_path` and `offtarget_self_fasta` are not
accepted over HTTP.** The API has no authentication, so it must not take an
executable path or an arbitrary file to read from the wire. The CLI keeps both
flags. Resolution of the design's own transcripts still runs here - it reads the
uploaded target and panel, not the filesystem - and
`offtarget_self_accessions`, `resolve_offtarget_self` and
`offtarget_self_coverage` are accepted as usual.

## Responses

A design response carries the criteria actually applied, a funnel of how many
candidates cleared each gate, and the selected candidates with their full
per-criterion breakdown - the same lossless state the CLI renders, so a client
never has to recompute a verdict.

`offtarget_skipped` reports whether the run was unreferenced. Treat it as
load-bearing: an empty off-target hit count means "not checked", not "nothing
found". See [design-rules.md](design-rules.md#transcriptome-off-target).

The criteria summary carries `offtarget_self_accessions` (the transcripts the
gate counted as the design's own), with `offtarget_self_declared` /
`offtarget_self_resolved` saying how they were found and
`offtarget_self_coverage` the threshold used. An empty list on an active filter
means every hit counted against every candidate, which on a gene family is a
broken run rather than a clean one - see
[What the design owns](design-rules.md#what-the-design-owns).

A convert response separates what survived from what was withheld:

| Field | Meaning |
|---|---|
| `constructs` | Constructs that passed the conversion filters |
| `withheld` | Constructs withheld, with their QC flags intact |
| `construct_count` / `withheld_count` | Their sizes |
| `converted_count` | The total built, before filtering |

`construct_count + withheld_count == converted_count` always holds, so nothing
disappears silently between design and conversion.

## Errors

Domain errors (`HushError`) and incoherent requests (`ValueError`) become
**400** with a `{"detail": "..."}` body. An unknown session id is **404**.

## Running it safely

There is no authentication, no rate limiting and no TLS. It is built to run on
your own machine. Do not publish the port beyond localhost without putting your
own authentication and a reverse proxy in front of it.

CORS is permissive only when the API is not serving a web client; when it serves
one, every call is same-origin and the middleware is not installed at all.
`create_app(allow_origins=[...])` overrides this.
