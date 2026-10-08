"""Render figures and a report only from completed, saved analysis stages."""
from pathlib import Path
import html
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D

from .common import publish, save_json


COLORS = {'Rod': '#3b6fb6', 'Cone': '#e09232', 'Bipolar': '#54a586',
          'Amacrine': '#a584b9', 'Horizontal': '#907048', 'RGC': '#cf5365',
          'Muller glia': '#35a0ac', 'Astrocyte': '#dd87af', 'Microglia': '#828c39',
          'Unassigned': '#aaaaaa'}


def read_summary(stage):
    return json.loads((Path(stage) / 'summary.json').read_text())


def save_figure(directory, stem, fig):
    publish(Path(directory) / (stem + '.png'), lambda path: fig.savefig(path, dpi=180, bbox_inches='tight', facecolor='white'))
    publish(Path(directory) / (stem + '.pdf'), lambda path: fig.savefig(path, bbox_inches='tight', facecolor='white'))
    plt.close(fig)


def qc_figure(stages, directory, config):
    frame = pd.read_csv(Path(stages['preprocess']) / 'cell_qc_audit.csv')
    donors = sorted(frame['donor'].unique())
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.9), constrained_layout=True)
    for ax, column, title in zip(axes, ['total_counts', 'n_genes_by_counts', 'pct_counts_mt'], ['UMIs per cell', 'Detected genes per cell', 'Mitochondrial UMIs (%)']):
        values = [frame.loc[frame['donor'].eq(d), column].to_numpy() for d in donors]
        boxes = ax.boxplot(values, tick_labels=donors, showfliers=False, patch_artist=True, widths=0.6)
        for box, color in zip(boxes['boxes'], ['#4178b9', '#d78c40', '#59a486']):
            box.set_facecolor(color)
            box.set_alpha(0.75)
        ax.set_title(title)
        if column != 'pct_counts_mt':
            ax.set_yscale('log')
        if column == 'n_genes_by_counts':
            ax.axhline(config['qc']['primary']['min_genes'], color='#bd4343', linestyle='--', label='Primary threshold')
            ax.axhline(config['qc']['strict']['min_genes'], color='#777777', linestyle=':', label='Strict threshold')
        if column == 'pct_counts_mt':
            ax.axhline(config['qc']['primary']['max_mito_percent'], color='#bd4343', linestyle='--')
            ax.axhline(config['qc']['strict']['max_mito_percent'], color='#777777', linestyle=':')
        ax.grid(axis='y', alpha=0.15)
    fig.suptitle('Quality control of the published filtered retinal UMI matrix', fontsize=13)
    save_figure(directory, '01_quality_control', fig)


def umap_figure(stages, directory):
    frame = pd.read_csv(Path(stages['embedding']) / 'cell_embeddings.csv')
    fig, axes = plt.subplots(1, 3, figsize=(16.8, 5.5), constrained_layout=True)
    for ax, column, title in zip(axes, ['donor', 'published_cell_type', 'leiden_0.5'], ['Donor', 'Published cell labels', 'New Leiden clusters (resolution 0.5)']):
        groups = sorted(frame[column].astype(str).unique(), key=lambda x: int(x) if x.isdigit() else x)
        palette = plt.get_cmap('tab20')
        for i, label in enumerate(groups):
            selected = frame[column].astype(str).eq(label)
            color = COLORS[label] if column == 'published_cell_type' else palette(i % 20)
            ax.scatter(frame.loc[selected, 'UMAP1'], frame.loc[selected, 'UMAP2'], s=1.6, c=[color], alpha=0.65, rasterized=True, linewidths=0)
        handles = [Line2D([0], [0], marker='o', linestyle='', markersize=4, color=COLORS[label] if column == 'published_cell_type' else palette(i % 20), label=label) for i, label in enumerate(groups)]
        ax.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, -0.10), ncol=3 if len(groups) > 10 else 2, fontsize=7, frameon=False, columnspacing=1.0)
        ax.set(title=title, xlabel='UMAP 1', ylabel='UMAP 2')
        ax.set_xticks([])
        ax.set_yticks([])
        ax.spines[['top', 'right']].set_visible(False)
    fig.suptitle('Count-level retinal single-cell analysis | 19,865 cells, three healthy donors', fontsize=14)
    save_figure(directory, '02_retina_umap', fig)


def marker_figure(stages, directory, config):
    frame = pd.read_csv(Path(stages['aggregate']) / 'marker_expression_published_types.csv')
    genes = list(dict.fromkeys(g for panel in config['markers'].values() for g in panel))
    types = config['cell_types'] + ['Unassigned']
    means = frame.pivot(index='group', columns='gene', values='mean_log1p_CP10k').reindex(index=types, columns=genes)
    detect = frame.pivot(index='group', columns='gene', values='pct_detected').reindex(index=types, columns=genes)
    relative = means / means.max(axis=0).replace(0, np.nan)
    fig, ax = plt.subplots(figsize=(13.5, 5.2), constrained_layout=True)
    x, y = np.meshgrid(np.arange(len(genes)), np.arange(len(types)))
    dots = ax.scatter(x.ravel(), y.ravel(), s=5 + detect.to_numpy().ravel() * 1.4, c=relative.to_numpy().ravel(), cmap='viridis', vmin=0, vmax=1, edgecolor='#555555', linewidths=0.2)
    n = frame.groupby('group')['n_cells'].first().reindex(types)
    ax.set_yticks(np.arange(len(types)), [f'{t} (n={int(n[t]):,})' for t in types])
    ax.set_xticks(np.arange(len(genes)), genes, rotation=60, ha='right')
    ax.invert_yaxis()
    ax.set_xlim(-0.6, len(genes) - 0.4)
    ax.set_title('Canonical-marker audit of published cell labels', pad=12)
    fig.colorbar(dots, ax=ax, fraction=0.025, pad=0.02, label='Mean log expression / gene maximum')
    ax.text(0, -0.40, 'Dot area: percentage of cells with >=1 UMI. Pooled descriptive audit; these are published labels.', transform=ax.transAxes, fontsize=9)
    ax.spines[['top', 'right']].set_visible(False)
    save_figure(directory, '03_canonical_marker_audit', fig)


def candidate_figure(stages, directory, config):
    frame = pd.read_csv(Path(stages['aggregate']) / 'balanced_gene_expression.csv')
    primary = frame.loc[frame['qc'].eq('primary')]
    types = config['cell_types']
    genes = config['target_genes']
    scores = primary.pivot(index='cell_type', columns='gene', values='donor_mean_log1p_pseudobulk_CPM').reindex(index=types, columns=genes)
    eligible = primary.pivot(index='cell_type', columns='gene', values='eligible_donors').reindex(index=types, columns=genes)
    scores = scores.where(eligible.eq(3))
    relative = scores / scores.max(axis=0).replace(0, np.nan)
    counts = pd.read_csv(Path(stages['aggregate']) / 'donor_cell_type_counts.csv')
    count_matrix = counts.loc[counts['qc'].eq('primary')].pivot(index='cell_type', columns='donor', values='n_cells').reindex(index=types)
    fig, axes = plt.subplots(1, 2, figsize=(12.4, 5.5), gridspec_kw={'width_ratios': [1.35, 1]}, constrained_layout=True)
    cmap = plt.get_cmap('viridis').copy()
    cmap.set_bad('#ededed')
    im = axes[0].imshow(relative, vmin=0, vmax=1, cmap=cmap, aspect='auto')
    axes[0].set_xticks(range(len(genes)), genes, rotation=45, ha='right')
    axes[0].set_yticks(range(len(types)), types)
    axes[0].set_title('Candidate expression: three-donor comparisons')
    for i, t in enumerate(types):
        for j, g in enumerate(genes):
            if not np.isfinite(relative.loc[t, g]):
                axes[0].text(j, i, 'NA', ha='center', va='center', color='#888888', fontsize=8)
    fig.colorbar(im, ax=axes[0], fraction=0.035, pad=0.02, label='Within-gene relative donor mean log1p CPM')
    axes[1].imshow(np.log1p(count_matrix), cmap='Blues', aspect='auto')
    axes[1].set_xticks(range(len(count_matrix.columns)), count_matrix.columns, rotation=30, ha='right')
    axes[1].set_yticks(range(len(types)), types)
    axes[1].set_title('Retained cells per donor and cell type')
    for i, t in enumerate(types):
        for j, d in enumerate(count_matrix.columns):
            n = int(count_matrix.loc[t, d])
            axes[1].text(j, i, f'{n:,}' + ('*' if n < 20 else ''), ha='center', va='center', fontsize=9, color='white' if n > 1000 else '#333333')
    fig.suptitle('Six genes nominated by the earlier glaucoma colocalisation pilot', fontsize=13)
    axes[0].text(0, -0.35, 'NA: not eligible in all three donors or gene absent. * Fewer than 20 cells; descriptive only.', transform=axes[0].transAxes, fontsize=9)
    save_figure(directory, '04_candidates_and_donor_coverage', fig)


def robustness_figure(stages, directory, config):
    summary = pd.read_csv(Path(stages['aggregate']) / 'robustness_summary.csv').set_index('gene')
    folds = pd.read_csv(Path(stages['aggregate']) / 'leave_one_donor_out.csv')
    donors = sorted(folds['omitted_donor'].unique())
    genes = config['target_genes']
    table = []
    for gene in genes:
        row = [summary.loc[gene, 'primary_top_type']]
        row += [folds.loc[folds['gene'].eq(gene) & folds['omitted_donor'].eq(donor), 'fold_top_type'].iloc[0] for donor in donors]
        table.append(row)
    labels = ['All three donors'] + ['Omit ' + donor for donor in donors]
    types = config['cell_types']
    color_index = {t: i for i, t in enumerate(types)}
    values = np.array([[color_index.get(str(value), np.nan) for value in row] for row in table])
    from matplotlib.colors import ListedColormap
    cmap = ListedColormap([COLORS[t] for t in types])
    cmap.set_bad('#eeeeee')
    fig, ax = plt.subplots(figsize=(11.2, 4.6), constrained_layout=True)
    ax.imshow(values, cmap=cmap, vmin=-0.5, vmax=len(types)-0.5, aspect='auto', alpha=0.35)
    ax.set_xticks(range(len(labels)), labels)
    ax.set_yticks(range(len(genes)), genes)
    for i, row in enumerate(table):
        for j, value in enumerate(row):
            ax.text(j, i, str(value) if pd.notna(value) and value else 'Not evaluable', ha='center', va='center', fontsize=11)
    ax.set_title('Leave-one-donor-out stability of the highest expression context', pad=15)
    ax.text(0, -0.15, 'Equal-donor mean log1p pseudobulk CPM; fixed eligible cell types. Saved pseudobulks are reaggregated; clusters are not refitted.', transform=ax.transAxes, fontsize=9)
    save_figure(directory, '05_leave_one_donor_out', fig)


def report_text(stages, config):
    ingest = read_summary(stages['ingest'])
    qc = read_summary(stages['preprocess'])
    pca = read_summary(stages['pca'])
    graph = read_summary(stages['graph'])
    aggregate = read_summary(stages['aggregate'])
    robustness = pd.read_csv(Path(stages['aggregate']) / 'robustness_summary.csv')
    cluster_audit = pd.read_csv(Path(stages['aggregate']) / 'cluster_annotation_audit.csv')
    expression = pd.read_csv(Path(stages['aggregate']) / 'donor_gene_expression.csv')
    umi_totals = expression.loc[expression['qc'].eq('primary')].groupby('gene')['sum_gene_UMIs'].sum()
    sparse_counts = ', '.join(f'{g}: {int(umi_totals[g])} UMIs' for g in ['COL8A2', 'DGKG', 'SLC2A12'])
    matched = int(cluster_audit.loc[cluster_audit['resolution'].eq(0.5), 'marker_panel_matches_published_dominant'].sum())
    table = ['| Gene | Highest context (eligible types) | LODO folds agreeing | Strict-QC comparison |', '|---|---|---:|---|']
    for row in robustness.to_dict('records'):
        qc_result = 'Same' if row['QC_top_agrees'] else ('Changed' if row['QC_evaluable'] else 'Not evaluable')
        table.append(f"| {row['gene']} | {row['primary_top_type'] if pd.notna(row['primary_top_type']) else 'Not evaluable'} | {row['LODO_agree_folds']}/3 | {qc_result} |")
    stable = ', '.join(aggregate['LODO_stable_genes']) or 'None'
    qc_stable = ', '.join(aggregate['QC_stable_genes']) or 'None'
    common = ', '.join(aggregate['primary_shared_cell_types'])
    qc_common = ', '.join(aggregate['QC_comparable_cell_types'])
    min10 = ', '.join(aggregate['min10_sensitivity_shared_cell_types'])
    text = f"""# Count-level single-cell extension of the glaucoma pilot

Author portfolio: Mila Desi Anasanti. Analysis run: 6 October 2026. Repository: https://github.com/MDAnusamandiriGroup/glaucoma-cell-context

## What was actually executed

This module reanalyses the published filtered retinal UMI count matrix of Lukowski et al. (2019), rather than relying only on the earlier pilot's published cell-class averages. The earlier 228-gene colocalisation evidence table is frozen and used solely to define the six retina-supported targets. No GWAS association, colocalisation or Mendelian-randomisation model is refitted here. Source files were checked against the published MD5 checksums. All {ingest['cells']:,} cell barcodes, count totals and detected-gene counts match the published metadata exactly.

The input contains {ingest['cells']:,} cells, {ingest['genes']:,} genes, three healthy donors and five libraries. Primary QC retained {qc['primary_cells']:,} cells ({100 * qc['primary_cells'] / ingest['cells']:.1f}%); stricter QC retained {qc['strict_cells']:,}. Total counts were scaled to 10,000 and log1p transformed, with integer UMIs retained separately. Library-aware selection produced {pca['highly_variable_genes']:,} highly variable genes; {pca['PCs']} principal components, a 15-neighbour graph, Leiden clustering and UMAP were computed. Leiden resolutions 0.5 and 1.0 produced {graph['clusters']['0.5']} and {graph['clusters']['1.0']} clusters respectively. No integrated or batch-corrected embedding was fitted.

## Main results

All {len(aggregate['target_genes_in_count_matrix'])} nominated targets map by exact symbol to the count matrix; {aggregate['published_candidate_genes_exactly_mapped']}/{aggregate['published_candidate_genes_total']} genes from the complete earlier evidence table map exactly. Published original cell labels remain the primary labels for expression summaries. Independently fitted clusters were audited against canonical marker panels and published-label composition: {matched}/{graph['clusters']['0.5']} clusters at resolution 0.5 had the same provisional marker-panel label as their dominant published label. This audit is descriptive and does not replace expert annotation.

With the predefined minimum of 20 cells per donor and cell type, the cell types eligible in all three donors are: {common}. Rare cell types with inadequate coverage are reported as unavailable for the primary three-donor comparison. Their measured expression is retained in the full donor table for descriptive inspection; missing groups are never treated as zero expression.

The primary ranking uses the equal-donor mean of log1p pseudobulk CPM. UMIs are summed within donor and published cell type, divided by all-gene UMIs in that group, and transformed with log1p. Each donor has equal weight. The highest context below means the highest score among eligible cell types; it is not a causal target-cell assignment.

{chr(10).join(table)}

The highest context remains the same in all three leave-one-donor-out (LODO) reaggregations for {len(aggregate['LODO_stable_genes'])}/6 genes: {stable}. LODO removes one donor from the expression summary and reaggregates the two retained donor pseudobulks. It does not rerun QC, feature selection, clustering or UMAP. The set of eligible cell types remains fixed across folds.

Several target genes have sparse counts across the nine assigned cell types: {sparse_counts}. These totals exclude Unassigned cells. COL8A2 has no detected UMIs among the eligible comparison types after strict QC, so its strict-QC ranking is not evaluable. LODO rank agreement for a low-count gene does not imply a precisely estimated expression level or strong biological evidence.

Strict-QC comparisons use only cell types eligible in all three donors under both QC rules: {qc_common}. The highest context agrees for {len(aggregate['QC_stable_genes'])}/6 genes: {qc_stable}. The separately reported minimum-10-cell sensitivity includes: {min10}. These sensitivity outputs are retained even where they differ from the primary result; thresholds were not chosen to maximise agreement.

## Interpretation and limitations

This is a reproducible technical pilot using real single-cell counts and donor-level descriptive sensitivity analyses. It demonstrates count-matrix ingestion, QC, normalisation, feature selection, dimensionality reduction, graph clustering, marker auditing and provenance-aware reporting. It is not a new glaucoma case-control single-cell study, an ageing association study, or a completed spatial-transcriptomics analysis. The donors are healthy, donor number is three, age is confounded with donor, and cell composition is unequal. Cells are not independent biological replicates. No differential-expression p-values, disease/age estimates, spatial mapping, or new doublet classifier are produced.

The published source matrix was already filtered by its authors. This workflow starts from their filtered UMI matrix, not FASTQ files or unfiltered droplets. Published original labels are retained; provisional labels based on marker panels are only an audit aid. Library-aware feature selection does not remove donor or technical batch effects. Candidate expression and LODO stability do not establish causality, independent replication, therapeutic efficacy, or publication readiness for a Q1 journal. Changes from the earlier atlas-average ranking can reflect donor balancing, different estimands and exclusion of underrepresented classes.

## Reproduce and resume

Install Python 3.12 and the packages in `single_cell/requirements.txt`, then run `python single_cell/run_single_cell.py`. After source files have been cached, `python single_cell/run_single_cell.py --offline` validates and reuses completed stages. Stage fingerprints include inputs, relevant code, parameters and software versions; output SHA256 values are validated before reuse. Partial runs are preserved in separate attempts, existing results are never overwritten, and a process lock prevents concurrent writers. Frozen computed matrices are intentionally excluded from ordinary GitHub upload packages.

## Sources

Lukowski SW et al. A single-cell transcriptome atlas of the adult human retina. The EMBO Journal (2019). https://doi.org/10.15252/embj.2018100811

Original filtered UMI count matrix and metadata: https://zenodo.org/records/5515631

Scanpy official preprocessing and clustering tutorial: https://scanpy.readthedocs.io/en/stable/tutorials/basics/clustering.html

The earlier candidate table and its original scientific sources remain documented in the parent repository.
"""
    return text


def render_pdf(directory, text):
    from reportlab.lib import colors
    from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
    from reportlab.lib.enums import TA_LEFT
    from reportlab.lib.pagesizes import A4
    from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, PageBreak, Table, TableStyle
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name='BodySmall', fontName='Helvetica', fontSize=9.2, leading=13, spaceAfter=7))
    styles.add(ParagraphStyle(name='CellSmall', fontName='Helvetica', fontSize=8, leading=10))
    story = []
    lines = text.splitlines()
    index = 0
    while index < len(lines):
        line = lines[index]
        if line.startswith('|'):
            rows = []
            while index < len(lines) and lines[index].startswith('|'):
                if not lines[index].startswith('|---'):
                    rows.append([Paragraph(html.escape(cell.strip()), styles['CellSmall']) for cell in lines[index].strip('|').split('|')])
                index += 1
            table = Table(rows, colWidths=[80, 165, 83, 145], repeatRows=1)
            table.setStyle(TableStyle([('BACKGROUND', (0,0), (-1,0), colors.HexColor('#e5edf5')), ('VALIGN', (0,0), (-1,-1), 'TOP'), ('GRID', (0,0), (-1,-1), 0.3, colors.HexColor('#cccccc')), ('TOPPADDING', (0,0), (-1,-1), 6), ('BOTTOMPADDING', (0,0), (-1,-1), 6)]))
            story.extend([table, Spacer(1,10)])
            continue
        if line.startswith('# '):
            story.append(Paragraph(html.escape(line[2:]), styles['Title']))
        elif line.startswith('## '):
            story.append(Paragraph(html.escape(line[3:]), styles['Heading2']))
        elif line:
            story.append(Paragraph(html.escape(line).replace('`', ''), styles['BodySmall']))
        index += 1
    for stem, caption in [
        ('01_quality_control', 'Quality control by donor. Whiskers show the central distribution; outliers are omitted from display.'),
        ('02_retina_umap', 'The same computed UMAP coloured by donor, published label and new Leiden cluster. Donor-specific structure may remain.'),
        ('03_canonical_marker_audit', 'Pooled canonical-marker audit. These labels are retained from the original authors.'),
        ('04_candidates_and_donor_coverage', 'Primary candidate comparison and donor coverage. Undersized groups are marked explicitly.'),
        ('05_leave_one_donor_out', 'Descriptive LODO expression reaggregation. No clusters or embedding are refitted.')
    ]:
        story.append(PageBreak())
        story.append(Paragraph(stem.replace('_', ' ').title(), styles['Heading2']))
        from PIL import Image as PILImage
        with PILImage.open(Path(directory) / (stem + '.png')) as im:
            width, height = im.size
        scale = min(473 / width, 610 / height)
        story.append(Image(str(Path(directory) / (stem + '.png')), width=width * scale, height=height * scale))
        story.append(Spacer(1,14))
        story.append(Paragraph(caption, styles['BodySmall']))
    def page_number(canvas, doc):
        canvas.setFont('Helvetica', 8)
        canvas.drawString(36, 22, 'Mila Desi Anasanti | Glaucoma single-cell pilot | 6 October 2026')
        canvas.drawRightString(A4[0]-36, 22, str(doc.page))
    def writer(path):
        doc = SimpleDocTemplate(str(path), pagesize=A4, leftMargin=36, rightMargin=36, topMargin=35, bottomMargin=35, title='Glaucoma single-cell pilot report', author='Mila Desi Anasanti')
        doc.build(story, onFirstPage=page_number, onLaterPages=page_number)
    publish(Path(directory) / 'Glaucoma_Single_Cell_Pilot_Report.pdf', writer)


def run(stages, directory, config):
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'axes.titleweight': 'semibold'})
    qc_figure(stages, directory, config)
    umap_figure(stages, directory)
    marker_figure(stages, directory, config)
    candidate_figure(stages, directory, config)
    robustness_figure(stages, directory, config)
    text = report_text(stages, config)
    publish(Path(directory) / 'REPORT.md', lambda path: path.write_text(text, encoding='utf-8'))
    render_pdf(directory, text)
    save_json(Path(directory) / 'figure_summary.json', {'figures': 5, 'formats': ['PNG', 'PDF'], 'report': 'Glaucoma_Single_Cell_Pilot_Report.pdf', 'analysis_refitted_for_rendering': False})
