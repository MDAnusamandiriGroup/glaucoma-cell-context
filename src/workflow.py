"""Small, auditable integration of published glaucoma and retinal expression data."""
from __future__ import annotations

import hashlib
import json
import os
import platform
import urllib.request
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


PRIMARY_CLUSTERS = {
    "Rod": ["C0 rod PR", "C1 Rod PR", "C2 Rod PR", "C3 Rod PR", "C4 Rod PR", "C7 Rod PR"],
    "Cone": ["C10 Cone PR"],
    "Bipolar": ["C6 Bipolar", "C8 Bipolar", "C11 Bipolar"],
    "Amacrine": ["C17 Amacrine"],
    "Horizontal": ["C15 Horizontal"],
    "RGC": ["C12 RGC"],
    "Muller glia": ["C9 Muller glia"],
    "Astrocyte": ["C16 Astrocytes"],
    "Microglia": ["C13 Microglia"],
}
CCA_CLUSTERS = {
    "Rod": ["CCA0 Rod PR", "CCA1 Rod PR"],
    "Cone": ["CCA5 Cone PR"],
    "Bipolar": ["CCA2 Bipolar", "CCA6 Bipolar", "CCA7 Bipolar", "CCA8 Bipolar"],
    "Amacrine": ["CCA9 Amacrine"],
    "RGC": ["CCA12 RGC"],
    "Muller glia": ["CCA4 Muller glia"],
    "Microglia": ["CCA11 Microglia"],
}
COL_QC = "Pass_QC_QTL_FDR05_P1E04_GWAS_P1E05"
SHEETS = {
    "eCAVIAR": ("SuppData7_eCAVIAR_POAG_CrossAnc", "CLPP"),
    "enloc": ("SupplData10_enloc_POAG_CrossAnc", "RCP"),
}


def sha256(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def write_new(path: Path, content: str | bytes) -> None:
    """Create a new file, never silently replace an existing result."""
    path.parent.mkdir(parents=True, exist_ok=True)
    mode = "xb" if isinstance(content, bytes) else "x"
    handle = path.open(mode) if isinstance(content, bytes) else path.open(mode, encoding="utf-8", newline="")
    with handle:
        handle.write(content)


def save_csv(frame: pd.DataFrame, path: Path, index: bool = False) -> None:
    write_new(path, frame.to_csv(index=index, float_format="%.12g"))


def acquire(root: Path, config: dict, offline: bool) -> dict:
    records = {}
    for source in config["sources"]:
        path = root / "data" / "raw" / source["file"]
        if path.exists():
            if sha256(path) != source["sha256"]:
                raise ValueError(f"Input checksum mismatch: {path}. Preserve this file and inspect it.")
            print(f"Reusing verified input: {path.name}", flush=True)
        else:
            if offline:
                raise FileNotFoundError(f"Offline input missing: {path}. Run once without --offline.")
            path.parent.mkdir(parents=True, exist_ok=True)
            print(f"Downloading: {path.name}", flush=True)
            request = urllib.request.Request(source["url"], headers={"User-Agent": "glaucoma-cell-context/0.1"})
            with urllib.request.urlopen(request, timeout=60) as response:
                payload = response.read()
            if hashlib.sha256(payload).hexdigest() != source["sha256"]:
                raise ValueError(f"Remote content changed: {path.name}; no output was overwritten.")
            temporary = path.with_suffix(".xlsx.partial")
            write_new(temporary, payload)
            if not zipfile.is_zipfile(temporary):
                raise ValueError(f"Downloaded input is not a valid XLSX: {temporary}")
            # Hard link provides no-clobber publication on the same filesystem.
            os.link(temporary, path)
            temporary.unlink()
        if not zipfile.is_zipfile(path):
            raise ValueError(f"Invalid XLSX input: {path}")
        records[source["file"]] = {"sha256": sha256(path), "bytes": path.stat().st_size, "url": source["url"], "doi": source["doi"]}
    return records


def load_colocalization(path: Path, config: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    workbook = pd.ExcelFile(path)
    tables = []
    audit = []
    for method, (sheet, score) in SHEETS.items():
        original = pd.read_excel(workbook, sheet_name=sheet, header=1)
        required = ["GWAS Trait", "GWAS variant", "rsID_GWAS", "gene_id", "Gene symbol", "QTL Type", "Tissue", COL_QC, score]
        missing = set(required) - set(original.columns)
        if missing:
            raise ValueError(f"Missing columns in {sheet}: {sorted(missing)}")
        # The workbook contains a blank separator row. Never manufacture a gene for it.
        data = original.loc[original["GWAS Trait"].eq("POAG Cross-ancestry")].copy()
        if data[required].isna().any().any():
            raise ValueError(f"Incomplete genuine association record in {sheet}")
        qc = pd.to_numeric(data[COL_QC], errors="raise")
        posterior = pd.to_numeric(data[score], errors="raise")
        if not qc.isin([0, 1]).all() or not posterior.between(0, 1).all():
            raise ValueError(f"Invalid QC flag/posterior in {sheet}")
        table = pd.DataFrame({
            "locus": data["GWAS variant"],
            "rsid": data["rsID_GWAS"],
            "gene_id": data["gene_id"].str.replace(r"\.\d+$", "", regex=True),
            "gene": data["Gene symbol"],
            "qtl_type": data["QTL Type"],
            "tissue": data["Tissue"],
            "method": method,
            "score": posterior,
            "qc_pass": qc.eq(1),
        })
        table["baseline_pass"] = table["qc_pass"] & table["score"].gt(config["thresholds"][method])
        table["strict_pass"] = table["qc_pass"] & table["score"].gt(config["sensitivity_thresholds"][method])
        audit.append({"method": method, "records": len(table), "qc_pass_records": int(table["qc_pass"].sum()), "baseline_pass_records": int(table["baseline_pass"].sum()), "strict_pass_records": int(table["strict_pass"].sum()), "retina_records": int(table["tissue"].eq("Retina").sum())})
        tables.append(table)
    return pd.concat(tables, ignore_index=True), pd.DataFrame(audit)


def evidence_matrix(records: pd.DataFrame) -> pd.DataFrame:
    keys = ["locus", "gene_id", "gene", "qtl_type", "tissue"]
    result = records.groupby(keys + ["method"], as_index=False).agg(score=("score", "max"), baseline_pass=("baseline_pass", "any"), strict_pass=("strict_pass", "any"))
    wide = result.pivot(index=keys, columns="method", values=["score", "baseline_pass", "strict_pass"])
    wide.columns = [f"{method}_{value}" for value, method in wide.columns]
    wide = wide.reset_index()
    # Missing method scores stay NA: absence from a thresholded result table is not a fitted zero probability.
    for method in SHEETS:
        for flag in ["baseline_pass", "strict_pass"]:
            wide[f"{method}_{flag}"] = wide[f"{method}_{flag}"].eq(True)
    wide["baseline_union"] = wide["eCAVIAR_baseline_pass"] | wide["enloc_baseline_pass"]
    wide["strict_union"] = wide["eCAVIAR_strict_pass"] | wide["enloc_strict_pass"]
    wide["baseline_same_context_agreement"] = wide["eCAVIAR_baseline_pass"] & wide["enloc_baseline_pass"]
    wide["strict_same_context_agreement"] = wide["eCAVIAR_strict_pass"] & wide["enloc_strict_pass"]
    return wide.sort_values(keys).reset_index(drop=True)


def evidence_counts(wide: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for rule in ["baseline_union", "strict_union", "baseline_same_context_agreement", "strict_same_context_agreement"]:
        subset = wide.loc[wide[rule]]
        rows.append({"rule": rule, "loci": subset["locus"].nunique(), "genes": subset["gene_id"].nunique(), "locus_gene_pairs": len(subset[["locus", "gene_id"]].drop_duplicates()), "qtl_tissue_contexts": len(subset)})
    return pd.DataFrame(rows)


def load_expression(path: Path, sheet: str, cluster_map: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    table = pd.read_excel(path, sheet_name=sheet, index_col=0)
    # The source XLSX has gene symbols converted by Excel to dates. Some are
    # ambiguous duplicates. Exclude every non-text label; never guess a repair.
    valid_label = np.array([isinstance(label, str) and bool(label.strip()) for label in table.index])
    excluded = pd.DataFrame({
        "source_file": path.name,
        "excel_row": np.arange(2, len(table) + 2)[~valid_label],
        "source_label": [str(label) for label in table.index[~valid_label]],
        "reason": "Non-text gene label, including Excel date conversion; ambiguous identity, excluded without guessed repair.",
    })
    table = table.loc[valid_label].copy()
    if not table.index.is_unique or table.index.hasnans:
        raise ValueError(f"Non-unique or missing atlas gene symbols in {path.name}")
    selected = [cluster for values in cluster_map.values() for cluster in values]
    if set(selected) - set(table.columns):
        raise ValueError(f"Atlas cluster labels changed: {path.name}")
    values = table[selected].apply(pd.to_numeric, errors="raise")
    if not np.isfinite(values.to_numpy()).all() or (values < 0).any().any():
        raise ValueError(f"Invalid expression averages in {path.name}")
    # Give each broad cell class one column; six rod clusters must not get six votes.
    expression = pd.DataFrame({cell: values[clusters].mean(axis=1) for cell, clusters in cluster_map.items()}).rename_axis("gene")
    return expression, excluded


def cell_shares(expression: pd.DataFrame) -> pd.DataFrame:
    denominator = expression.sum(axis=1).replace(0, np.nan)
    return expression.div(denominator, axis=0)


def gene_priorities(wide: pd.DataFrame, expression: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    supported = wide.loc[wide["baseline_union"]].copy()
    supported["retina_eQTL"] = supported["tissue"].eq("Retina") & supported["qtl_type"].eq("eQTL")
    supported["strict_retina_eQTL"] = supported["retina_eQTL"] & supported["strict_union"]
    grouped = supported.groupby(["gene_id", "gene"], as_index=False).agg(n_loci=("locus", "nunique"), n_qtl_tissue_contexts=("locus", "size"), retina_eQTL=("retina_eQTL", "any"), strict_retina_eQTL=("strict_retina_eQTL", "any"), same_context_agreement=("baseline_same_context_agreement", "any"), strict_same_context_agreement=("strict_same_context_agreement", "any"), loci=("locus", lambda x: ";".join(sorted(set(x)))))
    # Exact symbols only. Unmapped lncRNA aliases remain visibly unmatched.
    joined = grouped.merge(expression, left_on="gene", right_index=True, how="left", validate="many_to_one")
    joined["atlas_matched"] = joined[list(expression.columns)].notna().all(axis=1)
    shares = cell_shares(expression)
    joined = joined.merge(shares.add_suffix("_share"), left_on="gene", right_index=True, how="left", validate="many_to_one")
    matched = joined["atlas_matched"] & joined[[f"{cell}_share" for cell in expression]].notna().all(axis=1)
    joined["highest_mean_expression_cell"] = pd.Series(pd.NA, index=joined.index, dtype="object")
    joined.loc[matched, "highest_mean_expression_cell"] = joined.loc[matched, list(expression.columns)].idxmax(axis=1)
    joined["RGC_vs_other_classes_log2"] = np.log2((joined["RGC"] + 0.01) / (joined[[c for c in expression if c != "RGC"]].mean(axis=1) + 0.01))
    joined = joined.sort_values(["retina_eQTL", "strict_retina_eQTL", "strict_same_context_agreement", "same_context_agreement", "gene"], ascending=[False, False, False, False, True]).reset_index(drop=True)
    usable = joined["atlas_matched"] & joined[[f"{cell}_share" for cell in expression]].notna().all(axis=1)
    cases = sorted(joined.loc[joined["retina_eQTL"] & usable, "gene"].drop_duplicates())
    if not cases:
        raise ValueError("No retina-supported genes map to the expression atlas")
    return joined, cases


def bh_adjust(p_values: np.ndarray) -> np.ndarray:
    values = np.asarray(p_values, dtype=float)
    order = np.argsort(values)
    adjusted = np.minimum.accumulate((values[order] * len(values) / np.arange(1, len(values) + 1))[::-1])[::-1]
    result = np.empty(len(values))
    result[order] = np.minimum(adjusted, 1)
    return result


def matched_reference(expression: pd.DataFrame, cases: list[str], config: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Exploratory gene-level reference, NOT LD-aware GWAS cell-type enrichment."""
    expression = expression.loc[expression.sum(axis=1).gt(0)].copy()
    bins = pd.qcut(expression.mean(axis=1).rank(method="first"), q=config["expression_match_bins"], labels=False)
    shares = cell_shares(expression)
    observed = shares.loc[cases].mean(axis=0).to_numpy()
    null_sum = np.zeros((config["n_permutations"], expression.shape[1]))
    rng = np.random.default_rng(config["seed"])
    audit = []
    for bin_number in range(config["expression_match_bins"]):
        n_cases = int(bins.loc[cases].eq(bin_number).sum())
        controls = bins.index[bins.eq(bin_number) & ~bins.index.isin(cases)]
        audit.append({"expression_decile": bin_number + 1, "candidate_genes": n_cases, "available_control_genes": len(controls)})
        if n_cases == 0:
            continue
        if len(controls) < n_cases:
            raise ValueError(f"Insufficient matched controls in expression decile {bin_number + 1}")
        control_values = shares.loc[controls].to_numpy()
        for draw in range(config["n_permutations"]):
            selected = rng.choice(len(controls), size=n_cases, replace=False)
            null_sum[draw] += control_values[selected].sum(axis=0)
    null_mean = null_sum / len(cases)
    p_values = (1 + (null_mean >= observed).sum(axis=0)) / (1 + len(null_mean))
    result = pd.DataFrame({"cell_class": expression.columns, "n_candidate_genes": len(cases), "observed_mean_class_share": observed, "matched_null_mean": null_mean.mean(axis=0), "matched_null_2.5_percentile": np.quantile(null_mean, 0.025, axis=0), "matched_null_97.5_percentile": np.quantile(null_mean, 0.975, axis=0), "one_sided_empirical_p": p_values, "BH_q_across_9_classes": bh_adjust(p_values)})
    result["observed_minus_null_mean"] = result["observed_mean_class_share"] - result["matched_null_mean"]
    return result, pd.DataFrame(audit)


def cca_sensitivity(primary: pd.DataFrame, cca: pd.DataFrame, cases: list[str]) -> tuple[pd.DataFrame, dict]:
    common_classes = list(CCA_CLUSTERS)
    common_genes = sorted(set(cases) & set(cca.index))
    first = cell_shares(primary.loc[common_genes, common_classes])
    second = cell_shares(cca.loc[common_genes, common_classes])
    valid = first.notna().all(axis=1) & second.notna().all(axis=1)
    first = first.loc[valid]
    second = second.loc[valid]
    table = pd.DataFrame({"gene": first.index, "primary_RGC_share_common_7_classes": first["RGC"], "CCA_RGC_share_common_7_classes": second["RGC"], "primary_highest_class_common_7": first.idxmax(axis=1), "CCA_highest_class_common_7": second.idxmax(axis=1)})
    table["highest_class_agrees"] = table["primary_highest_class_common_7"].eq(table["CCA_highest_class_common_7"])
    rho = float(spearmanr(first["RGC"], second["RGC"]).statistic)
    summary = {"n_genes": len(table), "n_highest_class_agrees": int(table["highest_class_agrees"].sum()), "RGC_share_spearman": rho, "note": "Same atlas donors; processing sensitivity, not independent replication. Both shares use the same seven classes."}
    return table.reset_index(drop=True), summary


def run_analysis(root: Path, output: Path, config: dict) -> dict:
    raw = root / "data" / "raw"
    records, audit = load_colocalization(raw / "hamel2024_supplement.xlsx", config)
    evidence = evidence_matrix(records)
    counts = evidence_counts(evidence)
    atlas, primary_label_audit = load_expression(raw / "lukowski2019_EV1.xlsx", "Dataset EV1", PRIMARY_CLUSTERS)
    cca, cca_label_audit = load_expression(raw / "lukowski2019_EV2.xlsx", "Dataset EV2", CCA_CLUSTERS)
    priorities, cases = gene_priorities(evidence, atlas)
    reference, match_audit = matched_reference(atlas, cases, config)
    cca_table, cca_summary = cca_sensitivity(atlas, cca, cases)
    tables = {
        "source_record_audit.csv": audit,
        "colocalization_records.csv": records,
        "evidence_by_locus_gene_qtl_tissue.csv": evidence,
        "evidence_rule_counts.csv": counts,
        "gene_evidence_and_retinal_context.csv": priorities,
        "retina_supported_genes.csv": priorities.loc[priorities["retina_eQTL"]],
        "exploratory_expression_matched_reference.csv": reference,
        "expression_matching_audit.csv": match_audit,
        "cca_processing_sensitivity.csv": cca_table,
        "atlas_excluded_nontext_gene_labels.csv": pd.concat([primary_label_audit, cca_label_audit], ignore_index=True),
    }
    for name, table in tables.items():
        save_csv(table, output / "tables" / name)
    save_csv(atlas, output / "tables" / "retinal_class_average_expression.csv", index=True)
    summary = {
        "analysis_type": "Secondary-data portfolio pilot; published GWAS-QTL results integrated with published single-cell cluster averages.",
        "GWAS_loci_in_original_POAG_cross_ancestry_study": 127,
        "atlas_cells_in_original_publication": 20009,
        "atlas_donors_in_original_publication": 3,
        "atlas_genes": len(atlas),
        "atlas_primary_excluded_nontext_gene_labels": len(primary_label_audit),
        "atlas_CCA_excluded_nontext_gene_labels": len(cca_label_audit),
        "atlas_classes": len(atlas.columns),
        "evidence_counts": counts.to_dict("records"),
        "candidate_genes_atlas_matched": int(priorities["atlas_matched"].sum()),
        "candidate_genes_atlas_unmatched": int((~priorities["atlas_matched"]).sum()),
        "retina_supported_genes": int(priorities["retina_eQTL"].sum()),
        "strict_retina_supported_genes": int(priorities["strict_retina_eQTL"].sum()),
        "retina_supported_unique_genes_used_in_reference": len(cases),
        "retina_genes_used_in_reference": cases,
        "exploratory_classes_with_BH_q_lt_0.05": reference.loc[reference["BH_q_across_9_classes"].lt(0.05), "cell_class"].tolist(),
        "CCA_processing_sensitivity": cca_summary,
        "retina_enloc_results_available": False,
        "retina_enloc_note": "The source result table contains no retina enloc rows. Method agreement here is assessable for listed GTEx contexts, not retina-specific validation.",
    }
    write_new(output / "summary.json", json.dumps(summary, indent=2, allow_nan=False) + "\n")
    return summary


def software_versions() -> dict:
    import scipy
    import matplotlib
    import openpyxl
    return {"python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__, "matplotlib": matplotlib.__version__, "openpyxl": openpyxl.__version__}
