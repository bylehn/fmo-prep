"""Tests for PIEDA plotting utilities."""

import pytest
import pandas as pd


@pytest.fixture
def sample_pieda_df():
    return pd.DataFrame({
        "I": [1, 1, 1, 2, 2],
        "IFRG": ["LIG001", "LIG001", "LIG001", "ALA001", "ALA001"],
        "J": [2, 3, 4, 3, 4],
        "JFRG": ["ALA001", "GLY002", "VAL003", "GLY002", "VAL003"],
        "R": [3.5, 4.1, 5.0, 6.2, 7.0],
        "Q": [0.0, 0.0, 0.0, 0.0, 0.0],
        "COMPONENT": ["ES", "ES", "ES", "ES", "ES"],
        "ENERGY": [-5.2, -1.1, -0.3, 0.8, -0.2],
        "TOTAL": [-5.2, -1.1, -0.3, 0.8, -0.2],
    })


@pytest.mark.skip(reason="process_csv not yet implemented")
def test_detect_ligand_fragments(sample_pieda_df):
    from fmo_prep.analysis.plots import detect_ligand_fragments

    id_to_name = {1: "LIG001", 2: "ALA001", 3: "GLY002", 4: "VAL003"}
    frag_ids, labels = detect_ligand_fragments(sample_pieda_df, id_to_name)
    assert frag_ids == {1}
    assert "LIG001" in labels
