"""Tests for GAMESS log file parser."""

import pytest
import pandas as pd
from pathlib import Path

FIXTURES = Path(__file__).parent.parent / "fixtures"


def test_parse_gamout_pieda_shape():
    from fmo_prep.analysis.parser import parse_gamout

    df = parse_gamout(FIXTURES / "sample_pieda.log", pieda=True)
    assert isinstance(df, pd.DataFrame)
    # 3 pairs × 5 components
    assert df.shape == (15, 9)
    assert set(df.columns) == {"I", "IFRG", "J", "JFRG", "R", "Q", "COMPONENT", "ENERGY", "TOTAL"}


def test_parse_gamout_pieda_components():
    from fmo_prep.analysis.parser import parse_gamout

    df = parse_gamout(FIXTURES / "sample_pieda.log", pieda=True)
    assert set(df["COMPONENT"].unique()) == {"ES", "EX", "CT", "DI", "SOL"}


def test_parse_gamout_pieda_fragment_names():
    from fmo_prep.analysis.parser import parse_gamout

    df = parse_gamout(FIXTURES / "sample_pieda.log", pieda=True)
    assert "ALA001" in df["IFRG"].values or "ALA001" in df["JFRG"].values
    assert "LIG001" in df["IFRG"].values or "LIG001" in df["JFRG"].values


def test_parse_gamout_pieda_energy_numeric():
    from fmo_prep.analysis.parser import parse_gamout

    df = parse_gamout(FIXTURES / "sample_pieda.log", pieda=True)
    assert pd.api.types.is_float_dtype(df["ENERGY"])
    assert pd.api.types.is_float_dtype(df["TOTAL"])


def test_parse_gamout_pieda_known_values():
    from fmo_prep.analysis.parser import parse_gamout

    df = parse_gamout(FIXTURES / "sample_pieda.log", pieda=True)
    # Pair (2,1) ES component should be -9235.440
    row = df[(df["I"] == 2) & (df["J"] == 1) & (df["COMPONENT"] == "ES")]
    assert len(row) == 1
    assert abs(row["ENERGY"].iloc[0] - (-9235.440)) < 0.001


def test_parse_gamout_nopieda_shape():
    from fmo_prep.analysis.parser import parse_gamout

    df = parse_gamout(FIXTURES / "sample_nopieda.log", pieda=False)
    assert df.shape == (3, 7)
    assert "COMPONENT" not in df.columns
    assert set(df.columns) == {"I", "IFRG", "J", "JFRG", "R", "Q", "ENERGY"}


def test_parse_gamout_nopieda_energy_numeric():
    from fmo_prep.analysis.parser import parse_gamout

    df = parse_gamout(FIXTURES / "sample_nopieda.log", pieda=False)
    assert pd.api.types.is_float_dtype(df["ENERGY"])


def test_parse_gamout_nopieda_known_values():
    from fmo_prep.analysis.parser import parse_gamout

    df = parse_gamout(FIXTURES / "sample_nopieda.log", pieda=False)
    # Pair (2,1) ENERGY = -0.033
    row = df[(df["I"] == 2) & (df["J"] == 1)]
    assert len(row) == 1
    assert abs(row["ENERGY"].iloc[0] - (-0.033)) < 0.001


def test_parse_gamout_empty_log(tmp_path):
    from fmo_prep.analysis.parser import parse_gamout

    empty = tmp_path / "empty.log"
    empty.write_text("no gamess output here\n")
    df = parse_gamout(empty, pieda=True)
    assert df.empty
    assert list(df.columns) == ["I", "IFRG", "J", "JFRG", "R", "Q", "COMPONENT", "ENERGY", "TOTAL"]
