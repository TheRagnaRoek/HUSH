"""Tests for the HTTP API (api.app).

Driven through Starlette's TestClient against a fresh app + store per test, so
sessions never leak between cases. The pipeline itself is covered by the service
and CLI suites; these cover the HTTP contract: status codes, JSON shape, session
lifecycle, and that domain errors become 400 (not 500).
"""

import pytest
from fastapi.testclient import TestClient

from api.app import create_app
from api.store import SessionStore

_ALLELE_1 = "HLA-DEMO*01:01"
# A short, unambiguous amino-acid sequence (auto-detected as protein).
_PROTEIN = "MKWVTFISLLFLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFE"


@pytest.fixture
def client():
    return TestClient(create_app(store=SessionStore()))


@pytest.fixture
def target(target_fasta):
    return target_fasta.read_text(encoding="ascii")


@pytest.fixture
def panel(panel_fasta):
    return panel_fasta.read_text(encoding="ascii")


@pytest.fixture
def demo_target(demo_lowgc_fasta):
    return demo_lowgc_fasta.read_text(encoding="ascii")


# --------------------------------------------------------------------------- #
# meta + inspect
# --------------------------------------------------------------------------- #


def test_health(client):
    r = client.get("/health")
    assert r.status_code == 200
    assert r.json()["status"] == "ok"


def test_openapi_is_served(client):
    r = client.get("/openapi.json")
    assert r.status_code == 200
    assert "/design" in r.json()["paths"]


def test_inspect_target_only(client, target):
    r = client.post("/inspect", json={"target_fasta": target})
    assert r.status_code == 200
    body = r.json()
    assert body["target_ready"] is True
    assert body["panel"] is None


def test_inspect_with_panel(client, target, panel):
    r = client.post("/inspect", json={"target_fasta": target, "panel_fasta": panel})
    assert r.status_code == 200
    assert r.json()["panel"]["count"] == 3


def test_inspect_protein_target_flags_not_ready(client):
    r = client.post("/inspect", json={"target_fasta": _PROTEIN})
    assert r.status_code == 200
    body = r.json()
    assert body["target_ready"] is False
    assert body["note"]


def test_inspect_accepts_bare_sequence(client):
    # A raw sequence (no FASTA header) is accepted as the target.
    r = client.post("/inspect", json={"target_fasta": "ACGUACGUACGUACGUACGU"})
    assert r.status_code == 200
    assert r.json()["target_ready"] is True


# --------------------------------------------------------------------------- #
# stateless design
# --------------------------------------------------------------------------- #


def test_design_bare_default_selects(client, demo_target):
    r = client.post("/design", json={"target_fasta": demo_target})
    assert r.status_code == 200
    body = r.json()
    assert body["candidate_count"] > 0
    assert body["funnel"]["selected"] == len(body["selected"])
    assert len(body["selected"]) > 0
    assert body["offtarget_skipped"] is True


def test_design_allele_specific_shape(client, target, panel):
    r = client.post(
        "/design",
        json={
            "target_fasta": target,
            "panel_fasta": panel,
            "params": {
                "use_seed_filter": False,
                "targeting": {"allele_specific": _ALLELE_1},
            },
        },
    )
    assert r.status_code == 200
    body = r.json()
    assert body["targeting"]["mode"] == "allele_specific"
    assert body["targeting"]["discriminating"] is True
    assert "allele" in body["funnel"]
    # A discriminating run carries the per-candidate allele block, now with the
    # full per-spared-allele breakdown, not just the weakest one.
    allele = body["selected"][0]["allele"]
    assert allele["spare_min_mismatch"] is not None
    assert allele["spare_count"] == 2
    assert len(allele["per_allele"]) == 2
    assert allele["discrimination"]["strongest_position"] in (True, False)
    # Placement and count are reported side by side, never combined: the client
    # gets every difference, not only the one the rank rests on.
    positions = allele["discrimination"]["all_guide_positions"]
    assert positions and positions[0] == allele["discrimination"]["guide_position"]


def test_design_target_specific_shape(client, target, panel):
    r = client.post(
        "/design",
        json={
            "target_fasta": target,
            "panel_fasta": panel,
            "params": {
                "use_seed_filter": False,
                "targeting": {"target_specific": True},
            },
        },
    )
    assert r.status_code == 200
    body = r.json()
    assert body["targeting"]["mode"] == "target_specific"
    assert body["targeting"]["hit"] == []
    assert body["targeting"]["discriminating"] is True


def test_design_inverted_gc_window_is_400(client, demo_target):
    r = client.post(
        "/design",
        json={
            "target_fasta": demo_target,
            "params": {"use_gc_filter": True, "gc_min": 0.6, "gc_max": 0.3},
        },
    )
    assert r.status_code == 400
    assert "detail" in r.json()


def test_design_unknown_allele_is_400(client, target, panel):
    r = client.post(
        "/design",
        json={
            "target_fasta": target,
            "panel_fasta": panel,
            "params": {"targeting": {"allele_specific": "NOT_IN_PANEL"}},
        },
    )
    assert r.status_code == 400


def test_design_unknown_overhang_is_400(client, demo_target):
    r = client.post(
        "/design",
        json={"target_fasta": demo_target, "params": {"overhang": "bogus"}},
    )
    assert r.status_code == 400


def test_design_protein_target_is_400(client):
    r = client.post("/design", json={"target_fasta": _PROTEIN})
    assert r.status_code == 400


def test_design_empty_target_is_400(client):
    r = client.post("/design", json={"target_fasta": "   "})
    assert r.status_code == 400


# --------------------------------------------------------------------------- #
# sessions
# --------------------------------------------------------------------------- #


def test_session_lifecycle(client, demo_target):
    created = client.post("/sessions", json={"target_fasta": demo_target})
    assert created.status_code == 201
    sid = created.json()["id"]
    assert created.json()["has_last_design"] is False

    got = client.get(f"/sessions/{sid}")
    assert got.status_code == 200

    designed = client.post(f"/sessions/{sid}/design", json={})
    assert designed.status_code == 200
    assert len(designed.json()["selected"]) > 0

    # After a design the session reports it, and convert works off it.
    assert client.get(f"/sessions/{sid}").json()["has_last_design"] is True
    converted = client.post(f"/sessions/{sid}/convert", json={"loop": "jensen"})
    assert converted.status_code == 200
    body = converted.json()
    # Conversion accounts for every selected candidate, kept or withheld.
    assert body["converted_count"] == len(designed.json()["selected"])
    assert body["construct_count"] + body["withheld_count"] == body["converted_count"]
    assert len(body["constructs"]) == body["construct_count"]
    assert len(body["withheld"]) == body["withheld_count"]

    assert client.delete(f"/sessions/{sid}").status_code == 204
    assert client.get(f"/sessions/{sid}").status_code == 404


def test_session_design_reuses_and_reparams(client, target, panel):
    sid = client.post(
        "/sessions", json={"target_fasta": target, "panel_fasta": panel}
    ).json()["id"]
    first = client.post(f"/sessions/{sid}/design", json={"use_seed_filter": False})
    second = client.post(
        f"/sessions/{sid}/design", json={"use_seed_filter": False, "length": 21}
    )
    assert first.status_code == 200 and second.status_code == 200
    assert second.json()["length"] == 21


def test_session_convert_without_design_is_400(client, demo_target):
    sid = client.post("/sessions", json={"target_fasta": demo_target}).json()["id"]
    r = client.post(f"/sessions/{sid}/convert", json={})
    assert r.status_code == 400


def test_unknown_session_is_404(client):
    assert client.get("/sessions/deadbeef").status_code == 404
    assert client.post("/sessions/deadbeef/design", json={}).status_code == 404
    assert client.delete("/sessions/deadbeef").status_code == 404


def test_store_evicts_least_recently_used():
    store = SessionStore(max_sessions=2)
    from service.session import HushSession

    a = store.create(HushSession())
    b = store.create(HushSession())
    store.get(a)  # touch a, so b is now least-recently-used
    c = store.create(HushSession())  # evicts b
    assert store.get(a) is not None
    assert store.get(c) is not None
    assert store.get(b) is None
    assert len(store) == 2


# --------------------------------------------------------------------------- #
# serving the built web client from the API's own origin
# --------------------------------------------------------------------------- #


@pytest.fixture
def fake_dist(tmp_path):
    """A minimal built-SPA layout: an index and one asset."""
    dist = tmp_path / "dist"
    (dist / "assets").mkdir(parents=True)
    (dist / "index.html").write_text("<!doctype html><div id=root></div>", encoding="ascii")
    (dist / "assets" / "app.js").write_text("console.log(1);", encoding="ascii")
    return dist


def test_no_web_dir_serves_api_only(monkeypatch):
    # A dist left lying around in a checkout must not start being served.
    # HUSH_WEB_DIR is cleared because the image sets it, and this asserts the
    # behaviour when nothing configures a client at all.
    monkeypatch.delenv("HUSH_WEB_DIR", raising=False)
    client = TestClient(create_app())
    assert client.get("/health").status_code == 200
    assert client.get("/").status_code == 404


def test_web_client_is_served_at_the_root(fake_dist):
    client = TestClient(create_app(web_dir=fake_dist))
    res = client.get("/")
    assert res.status_code == 200
    assert "<div id=root>" in res.text
    assert client.get("/assets/app.js").status_code == 200


def test_api_routes_still_win_over_the_web_client(fake_dist):
    client = TestClient(create_app(web_dir=fake_dist))
    res = client.get("/health")
    assert res.status_code == 200
    assert res.json()["status"] == "ok"


def test_client_routes_fall_back_to_the_index(fake_dist):
    # BrowserRouter owns unknown paths, so they are the client's to resolve.
    client = TestClient(create_app(web_dir=fake_dist))
    res = client.get("/some/client/route")
    assert res.status_code == 200
    assert "<div id=root>" in res.text


def test_web_client_does_not_serve_files_outside_its_root(fake_dist):
    secret = fake_dist.parent / "secret.txt"
    secret.write_text("do-not-serve", encoding="ascii")
    client = TestClient(create_app(web_dir=fake_dist))
    for probe in ("../secret.txt", "..%2f..%2fsecret.txt", "/etc/passwd"):
        res = client.get(f"/{probe}")
        assert "do-not-serve" not in res.text


def test_design_params_schema_carries_the_self_set():
    # The wire field is a list: one accession never covers a gene family.
    from api.schemas import DesignParamsSchema

    params = DesignParamsSchema(
        offtarget_self_accessions=["NM_002116", "NR_001434"],
        offtarget_self_coverage=0.25,
    ).to_domain()
    assert params.offtarget_self_accessions == ("NM_002116", "NR_001434")
    assert params.offtarget_self_coverage == 0.25
    assert params.resolve_offtarget_self is True
    # Defaults stay best-practice: resolution on, nothing named by hand.
    assert DesignParamsSchema().to_domain().offtarget_self_accessions == ()


def test_criteria_summary_reports_the_self_set():
    from api.serialize import criteria_summary
    from core.offtarget import BruteForceSearcher, SelfTranscripts
    from core.selection import SelectionCriteria

    own = SelfTranscripts.of(
        declared=["NM_002116.8"], resolved=["xm_041680767"], min_coverage=0.1
    )
    criteria = SelectionCriteria(
        use_offtarget_filter=True,
        offtarget_searcher=BruteForceSearcher([("NM_1", "ACGT")], max_mismatch=3),
        offtarget_self_transcripts=own,
    )
    summary = criteria_summary(criteria)
    assert summary.offtarget_self_accessions == ["NM_002116", "XM_041680767"]
    assert summary.offtarget_self_declared == 1
    assert summary.offtarget_self_resolved == 1
    assert summary.offtarget_self_coverage == 0.1
