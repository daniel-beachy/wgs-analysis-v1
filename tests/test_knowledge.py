from pathlib import Path

import duckdb

from wgs import knowledge
from wgs.config import Config
from wgs.knowledge import Built, Source, Store, Upstream


class Fake(Source):
    id = "fake"
    title = "Fake source"
    upstream_version = "v1"

    def probe(self, pin=None):
        return [Upstream(url="https://example.org/f", etag=pin or self.upstream_version)]

    def build(self, upstream, scratch: Path, out: Path) -> Built:
        p = out / "t.parquet"
        duckdb.connect().execute(f"COPY (SELECT 1 AS x) TO '{p}' (FORMAT parquet)")
        return Built(version=upstream[0].etag, tables={"t": p}, upstream=upstream)


def test_refresh_versions_skips_unchanged_and_prunes(tmp_path, monkeypatch):
    cfg = Config(workspace=tmp_path / "ws", search=[tmp_path])
    fake = Fake()
    monkeypatch.setattr("wgs.knowledge.registry.SOURCES", [fake])
    assert knowledge.refresh(cfg, tmp_path / "scratch")["fake"] == "updated (v1)"
    assert knowledge.refresh(cfg, tmp_path / "scratch")["fake"] == "up to date (v1)"
    store = Store(cfg)
    assert store.versions() == {"fake": "v1"} and store.table("fake", "t").exists()
    for v in ("v2", "v3"):
        fake.upstream_version = v
        knowledge.refresh(cfg, tmp_path / "scratch", keep=2)
    assert store.versions() == {"fake": "v3"}
    kept = [v["version"] for v in store.index()["sources"]["fake"]["versions"]]
    assert kept == ["v3", "v2"]
    assert not (store.root / "fake" / "v1").exists()
    # pinning an older version reproduces it and makes it current
    knowledge.refresh(cfg, tmp_path / "scratch", pin={"fake": "v0"}, keep=5)
    assert store.versions() == {"fake": "v0"}
    knowledge.refresh(cfg, tmp_path / "scratch", pin={"fake": "a0"}, keep=1)
    assert store.versions() == {"fake": "a0"} and store.table("fake", "t").exists()


def test_unpinning_reuses_stored_version(tmp_path, monkeypatch):
    cfg = Config(workspace=tmp_path / "ws", search=[tmp_path])
    fake = Fake()
    monkeypatch.setattr("wgs.knowledge.registry.SOURCES", [fake])
    knowledge.refresh(cfg, tmp_path / "scratch")
    knowledge.refresh(cfg, tmp_path / "scratch", pin={"fake": "v0"})
    built = []
    monkeypatch.setattr(fake, "build", lambda *a: built.append(1))
    assert knowledge.refresh(cfg, tmp_path / "scratch")["fake"] == "updated (v1)"
    assert not built and Store(cfg).versions() == {"fake": "v1"}
