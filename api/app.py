"""The HUSH FastAPI application.

Handlers delegate to service.HushSession and render with api.serialize. The
built web client is served from the same app when HUSH_WEB_DIR is set.
Endpoints, errors and CORS behaviour: docs/api.md.
"""

from __future__ import annotations

import os
import shutil
from pathlib import Path

from fastapi import FastAPI, HTTPException, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, JSONResponse

from api import inputs, schemas, serialize
from api.store import SessionStore
from core.errors import HushError
from core.reference import default_index_prefix, index_exists
from core.version import __version__
from service.params import DesignParams
from service.session import HushSession

_VERSION = __version__

WEB_DIR_ENV = "HUSH_WEB_DIR"


def _web_root(explicit: Path | None = None) -> Path | None:
    """The directory with a built index.html, or None."""
    if explicit is not None:
        root = Path(explicit)
    elif os.environ.get(WEB_DIR_ENV):
        root = Path(os.environ[WEB_DIR_ENV])
    else:
        # No implicit lookup: a stale web/dist must not be served.
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
    """Build the app; store, CORS origins and web directory can be injected."""
    app = FastAPI(title="HUSH API", version=_VERSION, description=_DESCRIPTION)
    store = store if store is not None else SessionStore()
    web_root = _web_root(web_dir)

    # With a served client every call is same-origin and needs no CORS;
    # without one, allow all origins for the Vite dev server.
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

    @app.get("/health", response_model=schemas.HealthResponse, tags=["meta"])
    def health() -> schemas.HealthResponse:
        return schemas.HealthResponse(
            status="ok",
            service="HUSH API",
            version=_VERSION,
            bowtie_available=shutil.which(DesignParams().bowtie_path) is not None,
            offtarget_reference=index_exists(default_index_prefix()),
        )

    @app.post("/inspect", response_model=schemas.InspectResponse, tags=["stateless"])
    def inspect(req: schemas.InspectRequest) -> schemas.InspectResponse:
        session = _build_session(req.target_fasta, req.panel_fasta)
        return serialize.inspect_response(session.design_input())

    @app.post("/design", response_model=schemas.DesignResponse, tags=["stateless"])
    def design(req: schemas.DesignRequest) -> schemas.DesignResponse:
        session = _build_session(req.target_fasta, req.panel_fasta)
        return serialize.design_response(session.design(req.params.to_domain()))

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
    """Serve the built SPA, with index.html as fallback for client-side routes.

    Registered after the API routes, so those take precedence.
    """
    root = web_root.resolve()
    index = root / "index.html"

    @app.get("/{spa_path:path}", include_in_schema=False)
    def web_client(spa_path: str) -> FileResponse:
        try:
            candidate = (root / spa_path).resolve()
        except OSError:
            return FileResponse(index)
        # Reject paths that escape the web root.
        if candidate.is_file() and candidate.is_relative_to(root):
            return FileResponse(candidate)
        return FileResponse(index)


app = create_app()
