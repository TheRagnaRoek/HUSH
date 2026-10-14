"""The HUSH FastAPI application.

Two ways in, one code path: the stateless endpoints build a throwaway session
per request, while the /sessions endpoints keep one alive so repeat designs
reuse its cached alignment. Every handler delegates to service.HushSession and
renders the result with api.serialize; domain errors (HushError) and incoherent
requests (ValueError) become 400, an unknown session id 400 -> 404.

When the built web client is present it is served from this same app, so a
container exposes one port and the browser talks to its own origin.

Run it with:  uvicorn api.app:app --reload   (or  python -m api)
"""

from __future__ import annotations

import os
from importlib import metadata
from pathlib import Path

from fastapi import FastAPI, HTTPException, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, JSONResponse

from api import inputs, schemas, serialize
from api.store import SessionStore
from core.errors import HushError
from service.session import HushSession

try:
    _VERSION = metadata.version("hush")
except metadata.PackageNotFoundError:  # pragma: no cover - editable/uninstalled
    _VERSION = "0.0.0"

# Where the built SPA lives. The container points this at the dist it copied
# in; a source checkout falls back to web/dist, so 'npm run build' is enough.
WEB_DIR_ENV = "HUSH_WEB_DIR"


def _web_root(explicit: Path | None = None) -> Path | None:
    """The directory holding a built index.html, or None if there is none."""
    if explicit is not None:
        root = Path(explicit)
    elif os.environ.get(WEB_DIR_ENV):
        root = Path(os.environ[WEB_DIR_ENV])
    else:
        # Nothing implicit: a dist/ left in a checkout must not start being
        # served (and drop CORS) just because it exists.
        return None
    return root if (root / "index.html").is_file() else None


_DESCRIPTION = (
    "Allele-aware siRNA/shRNA design over HTTP. Post a target (and optional "
    "panel) to /design for a one-shot run, or create a /sessions instance to "
    "design and convert repeatedly against a cached alignment."
)


def create_app(
    store: SessionStore | None = None,
    allow_origins: list[str] | None = None,
    web_dir: Path | None = None,
) -> FastAPI:
    """Build the app.

    A caller may inject a store (tests), CORS origins, or the directory
    holding the built web client.
    """
    app = FastAPI(title="HUSH API", version=_VERSION, description=_DESCRIPTION)
    store = store if store is not None else SessionStore()
    web_root = _web_root(web_dir)

    # Serving the client ourselves makes every call same-origin, so CORS is
    # not needed at all. Without a client, stay permissive so the Vite dev
    # server can call in from its own port.
    if allow_origins is None:
        allow_origins = [] if web_root is not None else ["*"]
    if allow_origins:
        app.add_middleware(
            CORSMiddleware,
            allow_origins=allow_origins,
            allow_methods=["*"],
            allow_headers=["*"],
        )

    @app.exception_handler(HushError)
    async def _hush_error(request: Request, exc: HushError) -> JSONResponse:
        return JSONResponse(status_code=400, content={"detail": str(exc)})

    @app.exception_handler(ValueError)
    async def _value_error(request: Request, exc: ValueError) -> JSONResponse:
        return JSONResponse(status_code=400, content={"detail": str(exc)})

    def _build_session(target_fasta: str, panel_fasta: str | None) -> HushSession:
        session = HushSession()
        session.set_target(inputs.parse_target(target_fasta))
        if panel_fasta is not None:
            session.set_panel(inputs.parse_panel(panel_fasta))
        return session

    def _require(sid: str) -> HushSession:
        session = store.get(sid)
        if session is None:
            raise HTTPException(status_code=404, detail=f"no session {sid!r}")
        return session

    def _session_response(sid: str, session: HushSession) -> schemas.SessionResponse:
        return schemas.SessionResponse(
            id=sid,
            target=serialize.sequence_summary(session.target),
            panel=serialize.panel_summary(session.panel) if session.has_panel else None,
            has_last_design=session.last_design is not None,
        )

    # -- meta -------------------------------------------------------------- #

    @app.get("/health", response_model=schemas.HealthResponse, tags=["meta"])
    def health() -> schemas.HealthResponse:
        return schemas.HealthResponse(status="ok", service="HUSH API", version=_VERSION)

    # -- stateless --------------------------------------------------------- #

    @app.post("/inspect", response_model=schemas.InspectResponse, tags=["stateless"])
    def inspect(req: schemas.InspectRequest) -> schemas.InspectResponse:
        session = _build_session(req.target_fasta, req.panel_fasta)
        return serialize.inspect_response(session.design_input())

    @app.post("/design", response_model=schemas.DesignResponse, tags=["stateless"])
    def design(req: schemas.DesignRequest) -> schemas.DesignResponse:
        session = _build_session(req.target_fasta, req.panel_fasta)
        return serialize.design_response(session.design(req.params.to_domain()))

    # -- sessions ---------------------------------------------------------- #

    @app.post(
        "/sessions",
        response_model=schemas.SessionResponse,
        status_code=201,
        tags=["sessions"],
    )
    def create_session(req: schemas.CreateSessionRequest) -> schemas.SessionResponse:
        session = _build_session(req.target_fasta, req.panel_fasta)
        return _session_response(store.create(session), session)

    @app.get(
        "/sessions/{sid}", response_model=schemas.SessionResponse, tags=["sessions"]
    )
    def get_session(sid: str) -> schemas.SessionResponse:
        return _session_response(sid, _require(sid))

    @app.post(
        "/sessions/{sid}/design",
        response_model=schemas.DesignResponse,
        tags=["sessions"],
    )
    def session_design(
        sid: str, params: schemas.DesignParamsSchema | None = None
    ) -> schemas.DesignResponse:
        session = _require(sid)
        request = params or schemas.DesignParamsSchema()
        return serialize.design_response(session.design(request.to_domain()))

    @app.post(
        "/sessions/{sid}/convert",
        response_model=schemas.ConvertResponse,
        tags=["sessions"],
    )
    def session_convert(
        sid: str, options: schemas.ConvertOptionsSchema | None = None
    ) -> schemas.ConvertResponse:
        session = _require(sid)
        request = options or schemas.ConvertOptionsSchema()
        return serialize.convert_response(session.convert(request.to_domain()))

    @app.delete("/sessions/{sid}", status_code=204, tags=["sessions"])
    def delete_session(sid: str) -> None:
        if not store.delete(sid):
            raise HTTPException(status_code=404, detail=f"no session {sid!r}")

    if web_root is not None:
        _serve_web_client(app, web_root)

    return app


def _serve_web_client(app: FastAPI, web_root: Path) -> None:
    """Serve the built SPA, falling back to index.html for client routes.

    Registered after every API route, so those still win. The client uses
    BrowserRouter with a catch-all route, so an unknown path is its to
    resolve, not a 404 - hence the fallback rather than StaticFiles alone.
    """
    root = web_root.resolve()
    index = root / "index.html"

    @app.get("/{spa_path:path}", include_in_schema=False)
    def web_client(spa_path: str) -> FileResponse:
        try:
            candidate = (root / spa_path).resolve()
        except OSError:
            return FileResponse(index)
        # is_relative_to keeps a crafted path from escaping the web root.
        if candidate.is_file() and candidate.is_relative_to(root):
            return FileResponse(candidate)
        return FileResponse(index)


app = create_app()
