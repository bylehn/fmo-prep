"""Interaction energy visualisation.

Unified port of:
- fmo-poc/scripts/plot_pieda.py       (single-system plots)
- fmo-poc/scripts/plot_pieda_delta.py (delta / comparative plots)

Main entry point for the CLI::

    run_analysis(csv_path, compare_csv, frag_info, ...)
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from fmo_prep.analysis.reports import write_summary

# Canonical amino acids plus common histidine protonation variants.
STANDARD_RESIDUES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "HID", "HIE", "HIP",
}

_COMPONENT_COLORS = {
    "ES":    "#4C72B0",
    "EX":    "#C44E52",
    "CT":    "#55A868",
    "DI":    "#4122A8",
    "SOL":   "#64B5CD",
    "TOTAL": "#808080",
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def extract_residue_code(fragment_label) -> str:
    if fragment_label is None:
        return ""
    token = str(fragment_label).strip().upper()
    match = re.match(r"^([A-Z]+)", token)
    return match.group(1) if match else token


def is_standard_residue(fragment_label) -> bool:
    return extract_residue_code(fragment_label) in STANDARD_RESIDUES


def _has_chain_ab_pairs(df: pd.DataFrame) -> bool:
    mask_ab = (df["I_chain"] == "A") & (df["J_chain"] == "B")
    mask_ba = (df["I_chain"] == "B") & (df["J_chain"] == "A")
    return bool((mask_ab | mask_ba).any())


def map_original_residues(
    frag_info: dict,
    pdb_renum: str | Path,
    pdb_orig: str | Path,
) -> dict:
    """Map renumbered residue labels back to original sequence numbers via coordinates."""
    orig_coords: dict[tuple, str] = {}
    with open(pdb_orig) as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                res_name = line[17:20].strip()
                res_seq  = line[22:26].strip()
                x, y, z  = float(line[30:38]), float(line[38:46]), float(line[46:54])
                orig_coords[(round(x, 3), round(y, 3), round(z, 3))] = f"{res_name}{res_seq}"

    renum_to_orig: dict[str, str] = {}
    with open(pdb_renum) as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                res_name  = line[17:20].strip()
                res_seq   = line[22:26].strip()
                x, y, z   = float(line[30:38]), float(line[38:46]), float(line[46:54])
                coord_key = (round(x, 3), round(y, 3), round(z, 3))
                renum_key = f"{res_name}{res_seq}"
                if renum_key not in renum_to_orig and coord_key in orig_coords:
                    renum_to_orig[renum_key] = orig_coords[coord_key]

    for info in frag_info.values():
        label = info["plot_label"]
        if label in renum_to_orig:
            info["plot_label"] = renum_to_orig[label]

    return frag_info


# ---------------------------------------------------------------------------
# Core data processing
# ---------------------------------------------------------------------------

def build_fragment_maps(df: pd.DataFrame, frag_info: dict):
    """Build fragment ID → name/plot-label maps.

    CSV columns (IFRG/JFRG) take priority; frag_info fills gaps.

    Returns:
        (id_to_name, id_to_plot_label) — both dicts keyed by int fragment ID.
    """
    id_to_name: dict[int, str] = {}
    id_to_plot_label: dict[int, str] = {}

    if "IFRG" in df.columns:
        for frag_id, frag_name in df[["I", "IFRG"]].dropna().drop_duplicates().itertuples(index=False):
            id_to_name[int(frag_id)] = str(frag_name).strip()
    if "JFRG" in df.columns:
        for frag_id, frag_name in df[["J", "JFRG"]].dropna().drop_duplicates().itertuples(index=False):
            id_to_name[int(frag_id)] = str(frag_name).strip()

    for frag_id, info in frag_info.items():
        id_to_name.setdefault(int(frag_id), str(info.get("name", info.get("plot_label", frag_id))).strip())
        id_to_plot_label[int(frag_id)] = str(info.get("plot_label", info.get("name", frag_id))).strip()

    for frag_id, frag_name in id_to_name.items():
        id_to_plot_label.setdefault(int(frag_id), frag_name)

    return id_to_name, id_to_plot_label


def detect_ligand_fragments(
    df: pd.DataFrame,
    id_to_name: dict,
    ligand_hint: Optional[str] = None,
):
    """Auto-detect ligand fragment IDs from the interaction DataFrame.

    If *ligand_hint* is given, match fragment names that start with it.
    Otherwise, pick non-standard-residue fragments; if multiple exist,
    choose the one with the largest total |ENERGY| against standard residues.

    Returns:
        (ligand_frag_ids: set[int], ligand_labels: list[str])
    """
    if ligand_hint:
        hint = ligand_hint.strip().upper()
        hinted = [fid for fid, name in id_to_name.items() if str(name).upper().startswith(hint)]
        if not hinted:
            candidates = ", ".join(sorted({n for n in id_to_name.values() if not is_standard_residue(n)}))
            raise ValueError(
                f"No fragment label matched --ligand '{ligand_hint}'. "
                f"Available non-standard labels: {candidates}"
            )
        return set(hinted), [id_to_name[fid] for fid in hinted]

    candidate_ids = [fid for fid, name in id_to_name.items() if not is_standard_residue(name)]
    if not candidate_ids:
        raise ValueError(
            "Could not auto-detect a ligand fragment. Provide --ligand explicitly or "
            "ensure IFRG/JFRG labels are present in the CSV."
        )
    if len(candidate_ids) == 1:
        fid = candidate_ids[0]
        return {fid}, [id_to_name[fid]]

    standard_ids = {fid for fid, name in id_to_name.items() if is_standard_residue(name)}
    scores: dict[int, float] = {}
    for fid in candidate_ids:
        mask = ((df["I"] == fid) & df["J"].isin(standard_ids)) | \
               ((df["J"] == fid) & df["I"].isin(standard_ids))
        scores[fid] = df.loc[mask, "ENERGY"].abs().sum()

    best = max(scores, key=scores.get)
    return {best}, [id_to_name[best]]


def process_csv(
    csv_file: str | Path,
    frag_info: dict,
    cov_threshold: float = 150.0,
    ligand_frag_ids=None,
    ligand_hint: Optional[str] = None,
    interaction_mode: str = "auto",
):
    """Process interaction CSV and return pivot DataFrames ready for plotting.

    Port of plot_pieda.py:process_csv().

    Returns:
        (pivot_df, heatmap_pivot, heatmap_all_pivot, name_map,
         ligand_frag_ids, mode_used, partner_label_str)
    """
    df = pd.read_csv(csv_file)

    if "COMPONENT" not in df.columns:
        df["COMPONENT"] = "TOTAL"
    df["ENERGY"] = pd.to_numeric(df["ENERGY"])
    df["I_chain"] = df["I"].map(lambda x: frag_info.get(x, {}).get("chain", "U"))
    df["J_chain"] = df["J"].map(lambda x: frag_info.get(x, {}).get("chain", "U"))

    id_to_name, id_to_plot_label = build_fragment_maps(df, frag_info)

    df_all = df.copy()
    df_all["I_label"] = df_all["I"].map(lambda x: id_to_plot_label.get(int(x), f"Frag {x}"))
    df_all["J_label"] = df_all["J"].map(lambda x: id_to_plot_label.get(int(x), f"Frag {x}"))

    # All-vs-all heatmap (symmetrised, covalent-filtered)
    hm_all = df_all.groupby(["I_label", "J_label"])["ENERGY"].sum().reset_index()
    hm_all_sym = hm_all.rename(columns={"I_label": "J_label", "J_label": "I_label"})
    hm_all_combined = pd.concat([hm_all, hm_all_sym]).drop_duplicates(subset=["I_label", "J_label"])
    hm_all_combined.loc[hm_all_combined["ENERGY"].abs() >= cov_threshold, "ENERGY"] = 0.0
    heatmap_all_pivot = hm_all_combined.pivot(index="I_label", columns="J_label", values="ENERGY").fillna(0)

    # Determine mode
    use_chain_mode = (
        interaction_mode == "chain"
        or (interaction_mode == "auto" and ligand_hint is None and _has_chain_ab_pairs(df))
    )

    if use_chain_mode:
        mask_ab = (df["I_chain"] == "A") & (df["J_chain"] == "B")
        mask_ba = (df["I_chain"] == "B") & (df["J_chain"] == "A")

        df_target = pd.concat([
            df[mask_ab].assign(
                Protein_Frag=df.loc[mask_ab, "I"],
                Partner_Frag=df.loc[mask_ab, "J"],
            ),
            df[mask_ba].assign(
                Protein_Frag=df.loc[mask_ba, "J"],
                Partner_Frag=df.loc[mask_ba, "I"],
            ),
        ])
        if df_target.empty:
            raise ValueError("No chain A/B protein-peptide interaction rows found.")

        partner_label_str = "Chain B (peptide)"
        mode_used = "chain"
        ligand_frag_ids = None
    else:
        if ligand_frag_ids is None:
            ligand_frag_ids, ligand_labels = detect_ligand_fragments(df, id_to_name, ligand_hint=ligand_hint)
        else:
            ligand_frag_ids = set(ligand_frag_ids)
            ligand_labels = sorted({id_to_name.get(fid, f"Frag {fid}") for fid in ligand_frag_ids})

        mask_i = df["I"].isin(ligand_frag_ids) & ~df["J"].isin(ligand_frag_ids)
        mask_j = df["J"].isin(ligand_frag_ids) & ~df["I"].isin(ligand_frag_ids)

        df_target = pd.concat([
            df[mask_i].assign(Protein_Frag=df.loc[mask_i, "J"], Partner_Frag=df.loc[mask_i, "I"]),
            df[mask_j].assign(Protein_Frag=df.loc[mask_j, "I"], Partner_Frag=df.loc[mask_j, "J"]),
        ])

        if df_target.empty:
            names = ", ".join(sorted({id_to_name.get(fid, str(fid)) for fid in ligand_frag_ids}))
            raise ValueError(f"No protein-ligand interaction rows found for: {names}")

        partner_label_str = ", ".join(ligand_labels)
        mode_used = "ligand"

    df_target.loc[df_target["ENERGY"].abs() >= cov_threshold, "ENERGY"] = 0.0
    df_target["Protein_Label"] = df_target["Protein_Frag"].map(
        lambda x: id_to_plot_label.get(int(x), f"Frag {x}")
    )
    df_target["Partner_Label"] = df_target["Partner_Frag"].map(
        lambda x: id_to_plot_label.get(int(x), f"Frag {x}")
    )

    agg_df = df_target.groupby(
        ["Protein_Frag", "Protein_Label", "COMPONENT"], as_index=False
    )["ENERGY"].sum()
    pivot_df = (
        agg_df.pivot(index="Protein_Frag", columns="COMPONENT", values="ENERGY")
        .fillna(0)
        .sort_index()
    )

    hm_df = df_target.groupby(["Protein_Label", "Partner_Label"])["ENERGY"].sum().reset_index()
    heatmap_pivot = hm_df.pivot(index="Protein_Label", columns="Partner_Label", values="ENERGY").fillna(0)

    name_map = {int(row["Protein_Frag"]): row["Protein_Label"] for _, row in agg_df.iterrows()}
    return pivot_df, heatmap_pivot, heatmap_all_pivot, name_map, ligand_frag_ids, mode_used, partner_label_str


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def _plot_heatmap_all(pivot: pd.DataFrame, title_prefix: str, output_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(14, 12))
    sns.heatmap(pivot, cmap="coolwarm", center=0, ax=ax, square=True)
    ax.collections[0].colorbar.set_label("kcal/mol", size=16)
    plt.title(f"{title_prefix}All-vs-All Interaction Heatmap (Non-Covalent)", fontsize=18, pad=15)
    plt.xlabel("Fragment J", fontsize=16)
    plt.ylabel("Fragment I", fontsize=16)
    plt.xticks(rotation=90, fontsize=12)
    plt.yticks(rotation=0, fontsize=12)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close(fig)


def _plot_heatmap_focused(
    pivot: pd.DataFrame,
    title_prefix: str,
    partner_label_str: str,
    mode_used: str,
    output_path: Path,
) -> None:
    if mode_used != "chain":
        return  # focused heatmap only meaningful for chain A/B case
    fig, ax = plt.subplots(figsize=(12, 10))
    sns.heatmap(pivot, cmap="coolwarm", center=0, ax=ax, square=True)
    ax.collections[0].colorbar.set_label("kcal/mol", size=16)
    plt.title(
        f"{title_prefix}Protein-Peptide Interaction Heatmap ({partner_label_str})",
        fontsize=18, pad=15,
    )
    plt.xlabel("Peptide Residue", fontsize=16)
    plt.ylabel("Protein Residue", fontsize=16)
    plt.xticks(rotation=90, fontsize=12)
    plt.yticks(rotation=0, fontsize=12)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close(fig)


def _plot_bar(
    target_df: pd.DataFrame,
    net_totals: pd.Series,
    title_prefix: str,
    net_energy_label: str,
    partner_label_str: str,
    mode_used: str,
    output_path: Path,
) -> None:
    plot_colors = [_COMPONENT_COLORS.get(str(c).strip(), "#888888") for c in target_df.columns]
    target_df = target_df.astype(float).fillna(0)

    fig, ax = plt.subplots(figsize=(14, 7))
    if not target_df.empty and len(target_df.columns) > 0:
        target_df.plot(kind="bar", stacked=True, color=plot_colors, ax=ax, edgecolor="black", linewidth=0.5)

    ax.scatter(range(len(net_totals)), net_totals.values, color="black", marker="D", s=40, zorder=5, label=net_energy_label)
    ax.axhline(0, color="black", linewidth=1, zorder=0)

    if mode_used == "chain":
        plt.title(f"{title_prefix}PIEDA Interactions: Protein vs Peptide ({partner_label_str})", fontsize=18, pad=15)
    else:
        plt.title(f"{title_prefix}PIEDA Interactions: Protein vs Ligand ({partner_label_str})", fontsize=18, pad=15)
    plt.xlabel("Protein Residue", fontsize=16)
    plt.ylabel(f"{title_prefix}Interaction Energy (kcal/mol)", fontsize=16)
    plt.xticks(rotation=45, ha="right", fontsize=13)

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, title="Component", bbox_to_anchor=(1.01, 1), loc="upper left")
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Top-level runner
# ---------------------------------------------------------------------------

def run_analysis(
    csv_path: Path,
    compare_csv: Optional[Path],
    frag_info: dict,
    ligand_hint: Optional[str],
    interaction_mode: str,
    significant_threshold: float,
    cov_threshold: float,
    output_dir: Path,
) -> None:
    """Top-level analysis runner called by the CLI.

    Generates stacked bar plot, heatmaps, and summary text file.

    Args:
        csv_path: Primary interaction CSV (from parse_gamout).
        compare_csv: Optional second CSV for delta analysis (csv2 − csv1).
        frag_info: Fragment info dict from parse_frag_map_file.
        ligand_hint: Optional ligand label hint for auto-detection.
        interaction_mode: 'auto', 'chain', or 'ligand'.
        significant_threshold: Min |net energy| for bar plot inclusion.
        cov_threshold: Covalent outlier suppression threshold.
        output_dir: Where to write PNG and TXT outputs.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    pivot1, heat1, heat1_all, name_map1, ligand_frag_ids, mode_used, partner_label_str = process_csv(
        csv_path,
        frag_info,
        cov_threshold=cov_threshold,
        ligand_hint=ligand_hint,
        interaction_mode=interaction_mode,
    )

    if mode_used == "chain":
        print(f"[+] Using legacy chain mode: protein chain A vs peptide chain B")
    else:
        print(f"[+] Using ligand fragment(s): {partner_label_str}")

    is_diff = compare_csv is not None

    if is_diff:
        pivot2, heat2, heat2_all, _, _, _, _ = process_csv(
            compare_csv,
            frag_info,
            cov_threshold=cov_threshold,
            ligand_frag_ids=ligand_frag_ids,
            interaction_mode=mode_used,
        )
        target_df    = pivot2.subtract(pivot1, fill_value=0)
        target_heat  = heat2.subtract(heat1, fill_value=0)
        target_heat_all = heat2_all.subtract(heat1_all, fill_value=0)
        title_prefix     = r"$\Delta$ "
        net_energy_label = r"Net $\Delta$H"
        base_name        = "delta_H"
    else:
        target_df    = pivot1
        target_heat  = heat1
        target_heat_all = heat1_all
        title_prefix     = ""
        net_energy_label = "Net Total"
        base_name        = Path(csv_path).stem

    # Heatmaps
    _plot_heatmap_all(target_heat_all, title_prefix, output_dir / f"{base_name}_heatmap_all.png")
    _plot_heatmap_focused(target_heat, title_prefix, partner_label_str, mode_used, output_dir / f"{base_name}_heatmap.png")

    # Bar plot
    if is_diff:
        net_totals = target_df.sum(axis=1)
        target_df_bar = target_df.copy()
    else:
        net_totals    = target_df.sum(axis=1)
        target_df_bar = target_df.copy()

    # Relabel index with residue names
    target_df_bar.index = [name_map1.get(i, f"Frag {i}") for i in target_df_bar.index]
    net_totals.index    = target_df_bar.index

    # Write summary before filtering
    write_summary(
        net_totals,
        output_dir,
        mode=mode_used,
        partner_label=partner_label_str,
        is_delta=is_diff,
    )

    # Filter for significant interactions in bar plot
    energy_mask   = net_totals.abs() > significant_threshold
    target_df_bar = target_df_bar[energy_mask]
    net_totals    = net_totals[energy_mask]

    _plot_bar(
        target_df_bar, net_totals, title_prefix, net_energy_label,
        partner_label_str, mode_used,
        output_dir / f"{base_name}_stacked_bar.png",
    )

    print(f"[+] Outputs written to: {output_dir}")
