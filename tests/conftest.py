import pytest


@pytest.fixture(autouse=True)
def isolate_from_local_config(tmp_path_factory, monkeypatch):
    """Never let the developer's wgs.local.toml or env leak into tests."""
    fake_repo = tmp_path_factory.mktemp("repo")
    monkeypatch.setattr("wgs.config.REPO_ROOT", fake_repo)
    for var in ("WGS_CONFIG", "WGS_DATA", "WGS_WORKSPACE"):
        monkeypatch.delenv(var, raising=False)
