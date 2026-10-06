"""Figures and a short report from saved analysis tables only."""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.workflow import write_new


def render(analysis_dir: Path, figure_dir: Path) -> None:
    figure_dir.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.spines.top": False, "axes.spines.right": False, "savefig.dpi": 180})
    counts = pd.read_csv(analysis_dir / "tables/evidence_rule_counts.csv")
    labels = ["Source thresholds\nUnion", "Stricter thresholds\nUnion", "Same QTL + tissue\nTwo-method support", "Stricter + same context\nTwo-method support"]
    fig, ax = plt.subplots(figsize=(10.5, 4.7), constrained_layout=True)
    x = np.arange(len(counts))
    for shift, column, color in [(-0.2, "genes", "#247BA0"), (0.2, "loci", "#E09028")]:
        bars = ax.bar(x + shift, counts[column], width=0.36, label=column.capitalize(), color=color)
        ax.bar_label(bars, padding=3)
    ax.set_xticks(x, labels)
    ax.set_ylim(0, counts[["genes", "loci"]].to_numpy().max() * 1.18)
    ax.set_ylabel("Count in POAG cross-ancestry result tables")
    ax.set_title("Glaucoma gene nomination depends on the evidence rule", loc="left", pad=16, fontweight="bold")
    ax.legend(frameon=False)
    fig.savefig(figure_dir / "01_evidence_rule_sensitivity.png")
    plt.close(fig)

    genes = pd.read_csv(analysis_dir / "tables/retina_supported_genes.csv")
    genes = genes.loc[genes["atlas_matched"]].sort_values(["strict_retina_eQTL", "gene"], ascending=[False, True])
    cells = ["RGC", "Muller glia", "Astrocyte", "Microglia", "Amacrine", "Bipolar", "Horizontal", "Cone", "Rod"]
    values = genes[[f"{cell}_share" for cell in cells]].to_numpy()
    fig, ax = plt.subplots(figsize=(10.1, max(5, 0.32 * len(genes) + 2)), constrained_layout=True)
    plot = ax.imshow(values, aspect="auto", cmap="viridis", vmin=0, vmax=max(0.5, np.nanmax(values)))
    labels = [row.gene + (" *" if row.strict_retina_eQTL else "") for row in genes.itertuples()]
    ax.set_yticks(np.arange(len(genes)), labels)
    ax.set_xticks(np.arange(len(cells)), cells, rotation=35, ha="right")
    ax.set_title("Retinal context of genes with published retina eQTL support\n* retina CLPP > 0.5; class averages give each broad class equal weight", loc="left", pad=14, fontsize=11)
    fig.colorbar(plot, ax=ax, shrink=0.8, label="Class mean / sum of nine class means")
    fig.savefig(figure_dir / "02_retinal_cell_context.png")
    plt.close(fig)

    reference = pd.read_csv(analysis_dir / "tables/exploratory_expression_matched_reference.csv")
    reference = reference.set_index("cell_class").loc[cells].reset_index()
    x = np.arange(len(reference))
    fig, ax = plt.subplots(figsize=(10.4, 5.1), constrained_layout=True)
    mean = reference["matched_null_mean"].to_numpy()
    lower = reference["matched_null_2.5_percentile"].to_numpy()
    upper = reference["matched_null_97.5_percentile"].to_numpy()
    ax.errorbar(x, mean, yerr=[mean - lower, upper - mean], fmt="o", color="#7B8490", capsize=5, label="Expression-matched reference: mean and 95% null interval")
    ax.scatter(x, reference["observed_mean_class_share"], color="#B65064", marker="D", s=45, label="Retina eQTL-supported gene set", zorder=4)
    ax.set_xticks(x, cells, rotation=35, ha="right")
    ax.set_ylabel("Mean class-expression share across genes")
    ax.set_title("Exploratory expression context against matched reference genes\n10,000 draws; this is not LD-aware GWAS cell-type enrichment", loc="left", pad=14, fontsize=11)
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    fig.savefig(figure_dir / "03_expression_matched_reference.png")
    plt.close(fig)

    cca = pd.read_csv(analysis_dir / "tables/cca_processing_sensitivity.csv")
    fig, ax = plt.subplots(figsize=(6.7, 5.8), constrained_layout=True)
    ax.scatter(cca["primary_RGC_share_common_7_classes"], cca["CCA_RGC_share_common_7_classes"], color="#247BA0", s=40)
    limit = max(0.2, cca[["primary_RGC_share_common_7_classes", "CCA_RGC_share_common_7_classes"]].to_numpy().max() * 1.15)
    ax.plot([0, limit], [0, limit], "--", color="#9A9A9A", linewidth=1)
    label_offsets = {"RP11-466F5.8": (-72, 24), "COL8A2": (0, 45), "SLC2A12": (8, 18), "DGKG": (-12, 18)}
    for row in cca.itertuples():
        offset = label_offsets.get(row.gene, (5, 5))
        ax.annotate(row.gene, (row.primary_RGC_share_common_7_classes, row.CCA_RGC_share_common_7_classes), xytext=offset, textcoords="offset points", fontsize=8, arrowprops={"arrowstyle": "-", "color": "#7B8490", "lw": 0.6} if row.gene in label_offsets else None)
    ax.set(xlim=(0, limit), ylim=(0, limit), xlabel="RGC expression share before CCA", ylabel="RGC expression share after CCA")
    ax.set_title("Atlas-processing sensitivity in the same donors\nBoth calculations use the seven shared cell classes", loc="left", pad=14, fontsize=11)
    fig.savefig(figure_dir / "04_cca_processing_sensitivity.png")
    plt.close(fig)


def make_report(analysis_dir: Path, figure_dir: Path) -> str:
    summary = json.loads((analysis_dir / "summary.json").read_text())
    counts = pd.read_csv(analysis_dir / "tables/evidence_rule_counts.csv")
    cca = summary["CCA_processing_sensitivity"]
    rows = ["| Evidence rule | Genes | Loci | Locus–gene pairs |", "|---|---:|---:|---:|"]
    for row in counts.itertuples():
        rows.append(f"| {row.rule.replace('_', ' ')} | {row.genes} | {row.loci} | {row.locus_gene_pairs} |")
    reference = pd.read_csv(analysis_dir / "tables/exploratory_expression_matched_reference.csv")
    reference_rows = ["| Cell class | Observed share | Matched reference mean | Empirical P | BH q |", "|---|---:|---:|---:|---:|"]
    for row in reference.itertuples():
        reference_rows.append(f"| {row.cell_class} | {row.observed_mean_class_share:.3f} | {row.matched_null_mean:.3f} | {row.one_sided_empirical_p:.4f} | {row.BH_q_across_9_classes:.4f} |")
    relative_figures = "../" + figure_dir.name
    significant = summary["exploratory_classes_with_BH_q_lt_0.05"]
    null_statement = ("No class crossed BH q < 0.05 in the exploratory matched-reference calculation." if not significant else "Classes crossing BH q < 0.05 in the exploratory matched-reference calculation: " + ", ".join(significant) + ". This is not a genome-wide enrichment result.")
    return f"""# Glaucoma evidence and retinal cell context: completed pilot

This is a secondary-data analysis of published results, completed on 5 October 2026. It integrates human POAG GWAS–QTL evidence with a separate healthy-human retinal single-cell expression reference. It does not fit a new GWAS, colocalization model, MR model, raw single-cell model, or spatial model.

## Evidence sensitivity

{chr(10).join(rows)}

The baseline uses the source thresholds (eCAVIAR CLPP > 0.01 or enloc RCP > 0.1) after the source QC filter. The stricter sensitivity uses > 0.5 separately for each method. These scores have different definitions and are not averaged or compared as a common scale. Two-method support requires the same locus, Ensembl gene, QTL type and tissue. It uses overlapping source data, so it is corroboration across methods rather than independent replication.

![Evidence sensitivity]({relative_figures}/01_evidence_rule_sensitivity.png)

## Independent expression reference

The reference publication profiled 20,009 cells from **three donors**, not 20,009 independent people. After excluding {summary['atlas_primary_excluded_nontext_gene_labels']} non-text gene labels, its averages contain {summary['atlas_genes']:,} usable genes. The source XLSX includes gene labels converted by Excel to dates, including ambiguous duplicate dates. We preserve a row-level exclusion audit and do not guess corrected symbols. We collapse 16 annotated clusters into nine broad cell classes using an unweighted mean within each class. This estimates a class-balanced expression profile, not a donor-level pseudobulk profile or a cell-count-weighted population mean. Two unassigned clusters are excluded.

Of the nominated genes, {summary['candidate_genes_atlas_matched']} map by exact symbol and {summary['candidate_genes_atlas_unmatched']} do not. Unmatched genes remain in the evidence table. There are {summary['retina_supported_genes']} gene-ID/symbol entries with retina eQTL support, {summary['strict_retina_supported_genes']} retaining retina CLPP > 0.5; {summary['retina_supported_unique_genes_used_in_reference']} unique genes map to the atlas for the reference comparison.

There are **no retina enloc rows** in the source result table. Missing retina method evidence is not a negative finding. A strict retina eCAVIAR result plus GTEx method agreement still does not demonstrate retina-specific two-method agreement.

![Retinal context]({relative_figures}/02_retinal_cell_context.png)

## Exploratory expression-matched reference

We match random control gene sets on deciles of mean retinal expression, exclude candidate genes from the control pool, and sample without replacement within each draw. We compare the mean of class-expression shares using 10,000 draws with a fixed seed. Empirical one-sided P = (1 + number of null means at least as large as observed)/(10,001); BH adjustment covers nine cell classes. The plotted null interval is **not** a confidence interval for a disease effect.

{null_statement}

{chr(10).join(reference_rows)}

![Matched reference]({relative_figures}/03_expression_matched_reference.png)

This reference does not control LD, gene length, QTL ascertainment, or correlation between genes from the same locus. Cell classes are compositional. Treat the calculation as an exploratory context check, not confirmation that glaucoma risk acts in a particular cell class.

## Processing sensitivity

The same publication also provides CCA-adjusted cluster averages. Using the seven classes available in both versions, {cca['n_highest_class_agrees']}/{cca['n_genes']} genes retain their highest-expression class; RGC-share Spearman correlation is {cca['RGC_share_spearman']:.3f}. The low-quality CCA10 rod cluster and unassigned CCA3 cluster are excluded according to the source description. Horizontal cells and astrocytes are unavailable after CCA and are omitted from **both** versions for this comparison.

![CCA sensitivity]({relative_figures}/04_cca_processing_sensitivity.png)

This is a sensitivity analysis within the same three donors, not independent replication. The atlas is healthy tissue; it cannot show glaucoma differential expression or protective response to injury.

## Practical interpretation

The completed result is an auditable evidence table for deciding what to investigate next. Retina regulatory evidence and retinal cell expression provide different information. A nearest gene or a highly expressed gene is not automatically the causal gene. Further work would require allele harmonization, locus-level fine-mapping/colocalization from full summary statistics, donor-aware disease or injury datasets, and experimental validation before proposing a therapeutic direction.

Sources: [Hamel et al. 2024](https://doi.org/10.1038/s41467-023-44380-y), Supplementary Data 7 and 10; [Lukowski et al. 2019](https://doi.org/10.15252/embj.2018100811), Dataset EV1 and EV2. The underlying POAG GWAS is [Gharahkhani et al. 2021](https://doi.org/10.1038/s41467-020-20851-4). All primary-data collection and original model fitting belong to those authors.
"""
