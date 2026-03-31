"""Tests for PIEDA plotting utilities."""

import pytest
import pandas as pd
from pathlib import Path

FIXTURES = Path(__file__).parent.parent / "fixtures"


@pytest.fixture
def pieda_df():
    """Minimal PIEDA DataFrame: ligand frag 3 interacting with protein frags 1 and 2."""
    return pd.DataFrame({
        "I":         [3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
        "IFRG":      ["LIG001"] * 10,
        "J":         [1, 1, 1, 1, 1, 2, 2, 2, 2, 2],
        "JFRG":      ["ALA001"] * 5 + ["GLY002"] * 5,
        "R":         [3.5] * 5 + [4.1] * 5,
        "Q":         [0.12] * 10,
        "COMPONENT": ["ES", "EX", "CT", "DI", "SOL"] * 2,
        "ENERGY":    [-3.1, 0.8, -2.1, -1.9, 0.0, -1.5, 0.3, -1.2, -1.3, 0.0],
        "TOTAL":     [-6.3] * 5 + [-3.7] * 5,
    })


@pytest.fixture
def frag_info():
    return {
        1: {"chain": "A", "plot_label": "ALA1",  "name": "ALA001"},
        2: {"chain": "A", "plot_label": "GLY2",  "name": "GLY002"},
        3: {"chain": "A", "plot_label": "LIG276", "name": "LIG001"},
    }


def test_detect_ligand_fragments_auto(pieda_df, frag_info):
    from fmo_prep.analysis.plots import build_fragment_maps, detect_ligand_fragments

    id_to_name, _ = build_fragment_maps(pieda_df, frag_info)
    # Only LIG001 is non-standard
    frag_ids, labels = detect_ligand_fragments(pieda_df, id_to_name)
    assert frag_ids == {3}
    assert "LIG001" in labels


def test_detect_ligand_fragments_hint(pieda_df, frag_info):
    from fmo_prep.analysis.plots import build_fragment_maps, detect_ligand_fragments

    id_to_name, _ = build_fragment_maps(pieda_df, frag_info)
    frag_ids, labels = detect_ligand_fragments(pieda_df, id_to_name, ligand_hint="LIG")
    assert frag_ids == {3}


def test_detect_ligand_fragments_hint_not_found(pieda_df, frag_info):
    from fmo_prep.analysis.plots import build_fragment_maps, detect_ligand_fragments

    id_to_name, _ = build_fragment_maps(pieda_df, frag_info)
    with pytest.raises(ValueError, match="No fragment label matched"):
        detect_ligand_fragments(pieda_df, id_to_name, ligand_hint="ZZZ")


def test_build_fragment_maps_from_df(pieda_df):
    from fmo_prep.analysis.plots import build_fragment_maps

    id_to_name, id_to_label = build_fragment_maps(pieda_df, {})
    assert id_to_name[3] == "LIG001"
    assert id_to_name[1] == "ALA001"


def test_process_csv_ligand_mode(tmp_path, pieda_df, frag_info):
    from fmo_prep.analysis.plots import process_csv

    csv = tmp_path / "test.csv"
    pieda_df.to_csv(csv, index=False)

    pivot, hm, hm_all, name_map, lig_ids, mode, partner = process_csv(
        csv, frag_info, ligand_hint="LIG"
    )
    assert mode == "ligand"
    assert lig_ids == {3}
    assert "LIG001" in partner
    assert isinstance(pivot, pd.DataFrame)
    # Should have protein frags 1 and 2 as rows
    assert set(pivot.index) == {1, 2}


def test_process_csv_cov_threshold_zeroes_large(tmp_path, frag_info):
    from fmo_prep.analysis.plots import process_csv

    df = pd.DataFrame({
        "I": [3, 3], "IFRG": ["LIG001", "LIG001"],
        "J": [1, 2], "JFRG": ["ALA001", "GLY002"],
        "R": [3.5, 4.1], "Q": [0.1, 0.1],
        "COMPONENT": ["TOTAL", "TOTAL"],
        "ENERGY": [500.0, -2.5],   # 500 should be zeroed by cov_threshold=150
        "TOTAL": [500.0, -2.5],
    })
    csv = tmp_path / "test.csv"
    df.to_csv(csv, index=False)

    pivot, _, _, _, _, _, _ = process_csv(csv, frag_info, cov_threshold=150.0, ligand_hint="LIG")
    # Frag 1 energy should be 0 (clamped), frag 2 should be -2.5
    assert pivot.loc[1, "TOTAL"] == 0.0
    assert abs(pivot.loc[2, "TOTAL"] - (-2.5)) < 0.01


def test_run_analysis_creates_files(tmp_path, pieda_df, frag_info):
    from fmo_prep.analysis.plots import run_analysis

    csv = tmp_path / "interactions.csv"
    pieda_df.to_csv(csv, index=False)
    out = tmp_path / "out"

    run_analysis(
        csv_path=csv,
        compare_csv=None,
        frag_info=frag_info,
        ligand_hint="LIG",
        interaction_mode="ligand",
        significant_threshold=0.0,
        cov_threshold=150.0,
        output_dir=out,
    )

    assert (out / "interactions_stacked_bar.png").exists()
    assert (out / "interactions_heatmap_all.png").exists()
    assert (out / "interactions_summary.txt").exists()
