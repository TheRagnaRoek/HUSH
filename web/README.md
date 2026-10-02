# HUSH Designer (web)

The browser GUI for HUSH: a guided, allele-aware siRNA/shRNA design workflow that
renders over the FastAPI backend. Built to match the Project HIDE wiki stack -
**React 19 + Vite + TypeScript + Bootstrap 5 / react-bootstrap + react-router** -
in the HIDE palette and typefaces, with Heidi as the design-wizard guide.

## Develop

```bash
cd web
npm install
npm run dev          # http://localhost:5173
```

The app opens on representative **sample data**, so it runs and is fully
clickable with no backend. To drive the real pipeline, start the API and point
the app at it:

```bash
# from the repo root, in a Python env:
pip install -e ".[api]"
python -m api                       # serves http://localhost:8000

# then, in web/ :
cp .env.example .env.local          # VITE_API_URL=http://localhost:8000
npm run dev
```

Paste a target (and optional panel) FASTA in step 1 and run a design; with a
reachable backend the Input -> Criteria(design) -> Construct(convert) flow hits
the real endpoints (via the `/sessions` flow, so the panel alignment is cached).
On any API failure it falls back to sample data and shows an offline banner.

## Scripts

| Script | Does |
|---|---|
| `npm run dev` | Vite dev server with HMR |
| `npm run build` | type-check then production build to `dist/` |
| `npm run preview` | serve the production build |
| `npm run typecheck` | `tsc --noEmit` |
| `npm run lint` | eslint (flat config, zero warnings) |
| `npm run format` | prettier |
| `npm run gen:api` | regenerate `src/api/schema.d.ts` from `$VITE_API_URL/openapi.json` |

## Layout

```
src/
  api/        typed client: types.ts (mirror of api/schemas.py), client.ts, hush.ts
  session/    SessionProvider + context + useSession (app state over the API)
  lib/        disc.ts (decode discrimination), alignment.ts (view model), sample.ts
  hooks/      useTheme (light default + persisted toggle)
  styles/     tokens.css (HIDE palette, light+dark), theme.css (Bootstrap vars),
              fonts.css, app.css (ported mockup component styles)
  components/ layout/ steps/ criteria/ candidates/ viewer/ common/
```

## Types / codegen

`src/api/types.ts` is a hand-maintained mirror of the backend models
(`api/schemas.py`) so the app type-checks without a running server. `npm run
gen:api` emits `src/api/schema.d.ts` from the live `/openapi.json` for drift
checks (and if you later switch to an openapi-fetch generated client). Keep
`types.ts` in step with `api/schemas.py` when the backend changes.

## Deferred (follow-ups, gated on backend additions)

- **fornac** real stem-loop fold - needs the dot-bracket added to
  `ConvertResponse` (the engine already computes the fold). Until then the fold
  is a schematic (`HairpinFold`).
- **Nightingale** panel-MSA + protein-position tab - needs the aligned panel
  windows surfaced. The alignment currently shows each spared allele's deciding
  mismatch (`buildAlignment`, marked `complete: false`); sample data shows the
  full view it will become.
- **SeqViz** richer construct/plasmid map.
- Off-target per-hit drill-down (which transcript, where) - needs the API to
  surface hit details it now collapses to counts.

Integration points are marked `TODO(viz)` / `TODO(api)` in the source.
