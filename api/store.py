"""In-memory session store: HushSession instances by opaque id.

A warm cache, lost on restart. Access is locked, and the store evicts the
least recently used session beyond its capacity.
"""

from __future__ import annotations

import threading
import uuid
from collections import OrderedDict

from service.session import HushSession

DEFAULT_MAX_SESSIONS = 128


class SessionStore:
    """Bounded, thread-safe LRU map of id to HushSession."""

    def __init__(self, max_sessions: int = DEFAULT_MAX_SESSIONS) -> None:
        self._sessions: OrderedDict[str, HushSession] = OrderedDict()
        self._lock = threading.Lock()
        self._max = max_sessions

    def create(self, session: HushSession) -> str:
        """Store a session under a new id, evicting the oldest if full."""
        sid = uuid.uuid4().hex
        with self._lock:
            self._sessions[sid] = session
            self._sessions.move_to_end(sid)
            while len(self._sessions) > self._max:
                self._sessions.popitem(last=False)
        return sid

    def get(self, sid: str) -> HushSession | None:
        """The session (marked most recently used), or None."""
        with self._lock:
            session = self._sessions.get(sid)
            if session is not None:
                self._sessions.move_to_end(sid)
            return session

    def delete(self, sid: str) -> bool:
        """Drop a session; True if it existed."""
        with self._lock:
            return self._sessions.pop(sid, None) is not None

    def __len__(self) -> int:
        with self._lock:
            return len(self._sessions)
