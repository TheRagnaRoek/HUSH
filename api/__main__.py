"""Run the API: python -m api.

Reads HUSH_API_HOST, HUSH_API_PORT and HUSH_API_RELOAD from the environment.
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
