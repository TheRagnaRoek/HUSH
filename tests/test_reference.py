"""Tests for the off-target reference helpers (core.reference)."""

import os
import shutil

import pytest

from core import reference
from core.errors import ReferenceBuildError


class TestLocatingTheIndex:
    def test_reference_dir_follows_the_environment(self, tmp_path, monkeypatch):
        monkeypatch.setenv(reference.REFERENCE_DIR_ENV, str(tmp_path))
        assert reference.reference_dir() == tmp_path
        assert reference.default_index_prefix() == tmp_path / reference.INDEX_NAME

    def test_reference_dir_falls_back_to_the_default(self, monkeypatch):
        monkeypatch.delenv(reference.REFERENCE_DIR_ENV, raising=False)
        assert reference.reference_dir() == reference.DEFAULT_REFERENCE_DIR

    def test_index_exists_keys_on_the_first_ebwt_shard(self, tmp_path):
        prefix = tmp_path / "idx"
        assert not reference.index_exists(prefix)
        (tmp_path / "idx.1.ebwt").write_bytes(b"")
        assert reference.index_exists(prefix)

    def test_describe_index_returns_none_without_a_manifest(self, tmp_path):
        assert reference.describe_index(tmp_path / "idx") is None

    def test_describe_index_reads_the_manifest(self, tmp_path):
        prefix = tmp_path / "idx"
        reference.manifest_path(prefix).write_text(
            '{"refseq_release": "226"}', encoding="ascii"
        )
        assert reference.describe_index(prefix)["refseq_release"] == "226"

    def test_describe_index_survives_a_corrupt_manifest(self, tmp_path):
        prefix = tmp_path / "idx"
        reference.manifest_path(prefix).write_text("{not json", encoding="ascii")
        assert reference.describe_index(prefix) is None


class TestBuildPreconditions:
    def test_missing_bowtie_build_is_reported_before_downloading(
        self, tmp_path, monkeypatch
    ):
        monkeypatch.setattr(shutil, "which", lambda name: None)
        with pytest.raises(ReferenceBuildError) as exc:
            reference.build_index(tmp_path)
        assert "bowtie-build" in str(exc.value)

    def test_an_unwritable_outdir_is_reported_with_the_fix(
        self, tmp_path, monkeypatch
    ):
        # A reference volume created by an older image carries that image's uid.
        # Reading still works, so this only surfaces on a rebuild - and without
        # this check it surfaces as a bare EACCES partway through a multi-GB
        # download.
        monkeypatch.setattr(shutil, "which", lambda name: "/usr/bin/bowtie-build")
        monkeypatch.setattr(os, "access", lambda path, mode: False)
        with pytest.raises(ReferenceBuildError) as exc:
            reference.build_index(tmp_path)
        message = str(exc.value)
        assert "not writable" in message
        assert "chown" in message

    def test_a_writable_outdir_passes_the_precondition(self, tmp_path, monkeypatch):
        # Guards the check itself: it must not reject a directory that is fine.
        monkeypatch.setattr(shutil, "which", lambda name: "/usr/bin/bowtie-build")
        monkeypatch.setattr(
            reference, "_read_release_number", lambda: "test"
        )

        def _no_network(base_url, workdir):
            raise ReferenceBuildError("reached the download step")

        monkeypatch.setattr(reference, "_fetch_transcripts", _no_network)
        with pytest.raises(ReferenceBuildError) as exc:
            reference.build_index(tmp_path)
        assert "reached the download step" in str(exc.value)
