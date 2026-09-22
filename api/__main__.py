"""Run the HUSH API: ``python -m api`` (or ``uvicorn api.app:app --reload``).

A thin launcher for local use. Host/port/reload are read from the environment
(HUSH_API_HOST, HUSH_API_PORT, HUSH_API_RELOAD) so a container or a dev shell
can point it wherever without code changes.
"""

from __future__ import annotations

import os


def main() -> None:
    import uvicorn

    reload = os.environ.get("HUSH_API_RELOAD", "").lower() in ("1", "true", "yes")
    uvicorn.run(
        "api.app:app",
        host=os.environ.get("HUSH_API_HOST", "127.0.0.1"),
        port=int(os.environ.get("HUSH_API_PORT", "8000")),
        reload=reload,
    )


if __name__ == "__main__":
    main()
