"""In-memory session store for the API's stateful mode.

Holds live HushSession instances behind opaque ids so a client can upload the
target and panel once and then design/convert by id, reusing the cached MSA
alignment. State lives only in this process: it is a warm cache, not a database,
and is lost on restart. Access is serialised (FastAPI runs sync endpoints in a
threadpool), and the store is capacity-bounded with least-recently-used
eviction so a long-lived server cannot grow without limit.
"""

from __future__ import annotations

import threading
import uuid
from collections import OrderedDict

from service.session import HushSession

DEFAULT_MAX_SESSIONS = 128


class SessionStore:
    """A bounded, thread-safe, LRU map of id -> HushSession."""

    def __init__(self, max_sessions: int = DEFAULT_MAX_SESSIONS) -> None:
        self._sessions: OrderedDict[str, HushSession] = OrderedDict()
        self._lock = threading.Lock()
        self._max = max_sessions

    def create(self, session: HushSession) -> str:
        """Store a session under a fresh id, evicting the oldest if full."""
        sid = uuid.uuid4().hex
        with self._lock:
            self._sessions[sid] = session
            self._sessions.move_to_end(sid)
            while len(self._sessions) > self._max:
                self._sessions.popitem(last=False)
        return sid

    def get(self, sid: str) -> HushSession | None:
        """Return the session, marking it most-recently-used, or None."""
        with self._lock:
            session = self._sessions.get(sid)
            if session is not None:
                self._sessions.move_to_end(sid)
            return session

    def delete(self, sid: str) -> bool:
        """Drop a session; return whether it existed."""
        with self._lock:
            return self._sessions.pop(sid, None) is not None

    def __len__(self) -> int:
        with self._lock:
            return len(self._sessions)
