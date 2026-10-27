import os
from pathlib import Path
import pytest

from paths import get_data_path, get_repository_path, get_visuals_path


def test_repository_paths_use_configured_root(monkeypatch, tmp_path):
    monkeypatch.setenv("LUXETENEBRAE_PATH", str(tmp_path))
    assert get_repository_path() == tmp_path
    assert get_data_path() == tmp_path / "data"
    assert get_visuals_path() == tmp_path / "visuals"


@pytest.mark.parametrize("value", [None, "", " ", "\t"])
def test_repository_path_defaults_to_checkout(monkeypatch, value):
    if value is None:
        monkeypatch.delenv("LUXETENEBRAE_PATH", raising=False)
    else:
        monkeypatch.setenv("LUXETENEBRAE_PATH", value)
    expected = Path(__file__).resolve().parents[1]
    assert get_repository_path() == expected
    assert os.environ["LUXETENEBRAE_PATH"] == str(expected)


@pytest.mark.parametrize("value", ["~/nested", "nested", "outer/../nested"])
def test_repository_path_normalizes_configured_root(monkeypatch, tmp_path, value):
    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("LUXETENEBRAE_PATH", value)
    assert get_repository_path() == tmp_path / "nested"
    assert get_data_path() == tmp_path / "nested" / "data"
    assert get_visuals_path() == tmp_path / "nested" / "visuals"
    assert os.environ["LUXETENEBRAE_PATH"] == value
    assert not (tmp_path / "nested").exists()