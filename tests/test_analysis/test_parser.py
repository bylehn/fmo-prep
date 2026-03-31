"""Tests for GAMESS log file parser."""

import pytest
import pandas as pd


@pytest.mark.skip(reason="parse_gamout not yet implemented")
def test_parse_gamout_pieda(tmp_path):
    from fmo_prep.analysis.parser import parse_gamout

    # Minimal GAMESS log snippet would go in fixtures/sample.log
    log = tmp_path / "sample.log"
    # TODO: write minimal log fixture
    df = parse_gamout(log, pieda=True)
    assert isinstance(df, pd.DataFrame)
    assert set(df.columns) >= {"I", "J", "COMPONENT", "ENERGY"}


@pytest.mark.skip(reason="parse_gamout not yet implemented")
def test_parse_gamout_nopieda(tmp_path):
    from fmo_prep.analysis.parser import parse_gamout

    log = tmp_path / "sample.log"
    df = parse_gamout(log, pieda=False)
    assert "COMPONENT" not in df.columns
    assert "ENERGY" in df.columns
