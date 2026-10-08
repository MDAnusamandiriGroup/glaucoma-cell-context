"""Create figures and a report from validated fit checkpoints; never refit."""
import argparse
import hashlib
import json
import shutil
from pathlib import Path
from xml.sax.saxutils import escape

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import findfont, FontProperties
import numpy as np
import pandas as pd
from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, Image, PageBreak

from validate_results import validate, DEFAULT_INDEX

ROOT = Path(__file__).resolve().parent
GENES = ['NPC2', 'LTBP2', 'YLPM1']
STUDIES = ['GCST90011766', 'GCST90011770']
LABELS = {'GCST90011766': 'European primary', 'GCST90011770': 'Combined-ancestry sensitivity'}
PALETTE = ['#ccd5df', '#a8b7c9', '#91bcc4', '#dc985e', '#216675']


def build_report(index_path=DEFAULT_INDEX):
    verification = validate(index_path)
    index = json.loads((ROOT / index_path).read_text())
    fit, prep = ROOT / index['fit'], ROOT / index['preprocess']
    primary = pd.read_csv(fit / 'primary_settings_posteriors.csv')
    sensitivity = pd.read_csv(fit / 'all_sensitivity_posteriors.csv')
    coverage = pd.read_csv(prep / 'shared_coverage_audit.csv')
    audit = pd.read_csv(prep / 'eqtl_source_audit.csv')
    source_qtl = pd.read_csv(ROOT / 'data/processed/eyegex_target_nominal_raw.tsv', sep='\t')
    plan = json.loads((ROOT / 'PLAN.json').read_text())
    versions = {'matplotlib': matplotlib.__version__, 'reportlab': __import__('reportlab').Version,
                'openpyxl': __import__('openpyxl').__version__}
    contract = {'fit_manifest_sha256': index['fit_manifest_sha256'],
                'preprocess_manifest_sha256': index['preprocess_manifest_sha256'],
                'report_code_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                'validation_code_sha256': hashlib.sha256((ROOT / 'validate_results.py').read_bytes()).hexdigest(),
                'method_sources_sha256': hashlib.sha256((ROOT / 'data/method_sources.json').read_bytes()).hexdigest(),
                'report_software': versions}
    fingerprint = hashlib.sha256(json.dumps(contract, sort_keys=True).encode()).hexdigest()[:16]
    out = ROOT / 'published' / ('report_' + fingerprint)
    manifest_path = out / 'manifest.json'
    if out.exists():
        if not manifest_path.is_file():
            raise RuntimeError('Preserve incomplete report: ' + str(out))
        manifest = json.loads(manifest_path.read_text())
        assert manifest['contract'] == contract
        for name, digest in manifest['outputs'].items():
            assert hashlib.sha256((out / name).read_bytes()).hexdigest() == digest, name
        print('REUSE report ' + str(out))
        return out
    out.mkdir(parents=True, exist_ok=False)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                         'axes.spines.top': False, 'axes.spines.right': False,
                         'figure.facecolor': 'white', 'axes.titleweight': 'bold'})

    fig, ax = plt.subplots(figsize=(9.4, 4.5))
    ordered = primary.set_index(['study', 'gene']).loc[[(s, g) for s in STUDIES for g in GENES]].reset_index()
    y = np.arange(6)
    left = np.zeros(6)
    for i in range(5):
        values = ordered['PP_H' + str(i)].to_numpy()
        ax.barh(y, values, left=left, color=PALETTE[i], label=['H0: neither', 'H1: GWAS only', 'H2: eQTL only', 'H3: distinct', 'H4: shared'][i], height=0.66)
        left += values
    ax.set_yticks(y, [g + (' | European' if s == STUDIES[0] else ' | Combined') for s, g in zip(ordered.study, ordered.gene)])
    ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ax.set_xlabel('Posterior probability under the single-causal-variant model')
    ax.set_title('Candidate gene comparison at the prespecified settings', loc='left', pad=17)
    for i, row in ordered.iterrows():
        ax.text(1.018, i, f'H4 = {row.PP_H4:.3g}', va='center', fontsize=10)
    ax.legend(ncol=3, loc='upper left', bbox_to_anchor=(0, -0.18), frameon=False, fontsize=9)
    ax.grid(axis='x', alpha=0.15)
    fig.subplots_adjust(left=0.24, right=0.82, bottom=0.27, top=0.85)
    fig.savefig(out / 'Coloc_candidate_posteriors.png', dpi=190)
    plt.close(fig)

    fig, axes = plt.subplots(4, 1, figsize=(8.6, 7.2), sharex=True)
    raw_gwas = pd.read_csv(ROOT / 'data/processed/gwas_GCST90011766_NPC2_region.tsv', sep='\t')
    axes[0].scatter(raw_gwas.base_pair_location / 1e6, -np.log10(raw_gwas.p_value), s=4, color='#344a60', alpha=0.55, rasterized=True)
    axes[0].set_title('Regional association patterns | GRCh38 chromosome 14', loc='left', pad=12)
    axes[0].text(0.012, 0.85, 'European POAG GWAS', transform=axes[0].transAxes, fontsize=9)
    for axis, gene in zip(axes[1:], GENES):
        raw = source_qtl.loc[(source_qtl.V1 == plan['genes'][gene]) & source_qtl.V10.between(plan['index_position']-1000000, plan['index_position']+1000000)].copy()
        matched = pd.read_csv(prep / (STUDIES[0] + '_' + gene + '.tsv'), sep='\t')
        present = raw.V10.isin(matched.position_GRCh38)
        axis.scatter(raw.loc[~present, 'V10']/1e6, -np.log10(raw.loc[~present, 'V12']), s=5, color='#bbbbbb', alpha=0.6, rasterized=True)
        axis.scatter(matched.position_GRCh38/1e6, -np.log10(matched.eqtl_p), s=5, color='#216675', alpha=0.6, rasterized=True)
        axis.text(0.012, 0.85, gene + ' | bulk-retina eQTL', transform=axis.transAxes, fontsize=9)
        if gene == 'NPC2':
            lead = raw.loc[raw.V12.idxmin()]
            axis.scatter([lead.V10/1e6], [-np.log10(lead.V12)], marker='x', color='#c44734', s=48, zorder=5)
            axis.annotate('Source lead indel absent in GWAS', xy=(lead.V10/1e6, -np.log10(lead.V12)), xytext=(0.36, 0.85), textcoords='axes fraction', fontsize=8, color='#9e3425', arrowprops={'arrowstyle': '-', 'color': '#9e3425', 'lw': 0.7})
    for axis in axes:
        axis.axvline(plan['index_position']/1e6, color='#b68848', ls='--', lw=1)
        axis.set_ylabel('−log10(P)', fontsize=9)
        axis.grid(alpha=0.1)
        axis.set_xlim((plan['index_position']-1000000)/1e6, (plan['index_position']+1000000)/1e6)
        axis.set_ylim(bottom=0)
    axes[-1].set_xlabel('Chromosome 14 position (Mb) | dashed line: rs754458')
    fig.text(0.12, 0.015, 'Teal: matched eQTL variants. Gray: source variants absent from the matched set. No LD colouring.', fontsize=8)
    fig.tight_layout(rect=(0, 0.045, 1, 1), h_pad=0.75)
    fig.savefig(out / 'Coloc_regional_associations.png', dpi=170)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4), sharey=True)
    gene_colors = {'NPC2': '#216675', 'LTBP2': '#aa6a32', 'YLPM1': '#775c91'}
    for axis, study in zip(axes, STUDIES):
        for gene in GENES:
            selected = sensitivity.loc[(sensitivity.study == study) & (sensitivity.gene == gene)]
            grouped = selected.groupby('p12').PP_H4.agg(['min', 'max'])
            working = selected.loc[(selected.half_window_bp == 1000000) & (selected.allele_policy == 'exact_all') & (selected.eqtl_prior_sd == 0.15)].sort_values('p12')
            axis.fill_between(grouped.index.to_numpy(), grouped['min'].to_numpy(), grouped['max'].to_numpy(), color=gene_colors[gene], alpha=0.12)
            axis.plot(working.p12, working.PP_H4, 'o-', color=gene_colors[gene], label=gene, lw=1.7, ms=4)
        axis.axhline(0.8, color='#777777', ls='--', lw=0.9)
        axis.axvline(1e-5, color='#777777', ls=':', lw=0.9)
        axis.set_xscale('log')
        axis.set_xlabel('Per-variant shared-association prior p12')
        axis.set_title(LABELS[study], fontsize=11, loc='left')
        axis.set_ylim(-0.02, 1.02)
        axis.grid(alpha=0.15)
    axes[0].set_ylabel('PP.H4')
    axes[1].legend(frameon=False, loc='lower right')
    fig.text(0.09, 0.01, 'Lines: working prior SD 0.15, ±1 Mb, all exact matches. Shading: other frozen window / allele / effect-prior settings.', fontsize=8)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    fig.savefig(out / 'Coloc_prior_sensitivity.png', dpi=190)
    plt.close(fig)

    for source, name in [(fit/'primary_settings_posteriors.csv', 'Coloc_primary_results.csv'),
                         (fit/'all_sensitivity_posteriors.csv', 'Coloc_sensitivity_results.csv'),
                         (prep/'shared_coverage_audit.csv', 'Coloc_variant_coverage.csv')]:
        shutil.copy2(source, out/name)
    with pd.ExcelWriter(out / 'Coloc_results_and_QC.xlsx', engine='openpyxl') as writer:
        primary.to_excel(writer, sheet_name='Primary settings', index=False)
        sensitivity.to_excel(writer, sheet_name='Sensitivity 216', index=False)
        coverage.to_excel(writer, sheet_name='Shared variant QC', index=False)
        audit.to_excel(writer, sheet_name='Full nominal source QC', index=False)
        pd.DataFrame([{'item': k, 'value': str(v)} for k, v in plan.items()]).to_excel(writer, sheet_name='Frozen analysis plan', index=False)
        for sheet in writer.book.worksheets:
            sheet.freeze_panes = 'C2'
            sheet.auto_filter.ref = sheet.dimensions
            for cells in sheet.columns:
                maximum = max(len(str(c.value or '')) for c in cells)
                sheet.column_dimensions[cells[0].column_letter].width = min(45, max(12, maximum+2))

    pdfmetrics.registerFont(TTFont('Pilot', findfont(FontProperties(family='DejaVu Sans'))))
    pdfmetrics.registerFont(TTFont('PilotBold', findfont(FontProperties(family='DejaVu Sans', weight='bold'))))
    pdfmetrics.registerFontFamily('Pilot', normal='Pilot', bold='PilotBold', italic='Pilot', boldItalic='PilotBold')
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name='PilotBody', fontName='Pilot', fontSize=9, leading=13, spaceAfter=8))
    styles.add(ParagraphStyle(name='PilotSmall', fontName='Pilot', fontSize=7.5, leading=10, spaceAfter=6))
    styles.add(ParagraphStyle(name='PilotTitle', fontName='PilotBold', fontSize=20, leading=25, textColor=colors.HexColor('#163f50'), spaceAfter=10))
    styles.add(ParagraphStyle(name='PilotHeading', fontName='PilotBold', fontSize=12, leading=16, spaceAfter=9, textColor=colors.HexColor('#163f50')))
    story = []
    def para(text, style='PilotBody'):
        story.append(Paragraph(text, styles[style]))
    def table(data, widths):
        rows = [[Paragraph(escape(str(cell)), styles['PilotSmall']) for cell in row] for row in data]
        item = Table(rows, colWidths=widths, repeatRows=1, hAlign=TA_LEFT)
        item.setStyle(TableStyle([('BACKGROUND', (0,0),(-1,0), colors.HexColor('#e4eef1')), ('VALIGN',(0,0),(-1,-1),'TOP'),('LINEBELOW',(0,0),(-1,0),0.5,colors.HexColor('#a7bcc4')),('BOTTOMPADDING',(0,0),(-1,-1),7),('TOPPADDING',(0,0),(-1,-1),6),('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#f5f7f8')])]))
        story.append(item)
        story.append(Spacer(1, 10))

    para('Targeted glaucoma–retina colocalization', 'PilotTitle')
    para('Completed portfolio extension • 6 October 2026<br/>Project owner: Mila Desi Anasanti', 'PilotSmall')
    para('Three candidate genes, full regional statistics, explicit uncertainty', 'PilotHeading')
    para('A new targeted Python ABF colocalization analysis compared NPC2, LTBP2 and YLPM1 around rs754458. The data are public EyeGEx bulk-retina eQTL statistics and two POAG GWAS releases. Six gene–study comparisons and 216 fixed sensitivity settings completed.')
    table([['Gene','European PP.H4','Combined PP.H4','Primary shared variants'],
           *[[g, f"{primary.loc[(primary.gene==g)&primary.primary_GWAS,'PP_H4'].iloc[0]:.6f}", f"{primary.loc[(primary.gene==g)&~primary.primary_GWAS,'PP_H4'].iloc[0]:.6f}", f"{primary.loc[(primary.gene==g)&primary.primary_GWAS,'shared_variants'].iloc[0]:,}"] for g in GENES]], [70,105,110,185])
    story.append(Image(str(out/'Coloc_candidate_posteriors.png'), width=470, height=225))
    para('No gene exceeds the prespecified PP.H4 ≥ 0.80 reference in the European primary analysis. LTBP2 reaches 0.809 in the overlapping combined-ancestry analysis, but none of the genes has H4 ≥ 0.80 across all frozen sensitivity settings. YLPM1 favours H3 (different associated variants), with European PP.H3 = 0.957.', 'PilotBody')
    para('PP.H4 is posterior support for a shared association signal under the chosen model and priors. It is not a probability that a gene is causal, a therapeutic target, or an application outcome. The combined GWAS includes the European analysis and is not independent replication.', 'PilotSmall')
    story.append(PageBreak())

    para('Inputs, harmonization and validation', 'PilotTitle')
    table([['Input','Cohort / source','Role'],['European POAG','16,677 cases; 199,580 controls\nGCST90011766','Primary GWAS'],['Combined POAG','34,179 cases; 349,321 controls\nGCST90011770','Overlapping ancestry sensitivity'],['EyeGEx retina','406 individuals; nominal df = 404\nRatnapriya et al. 2019 / author Zenodo replacement','Bulk-retina expression QTL']], [90,225,155])
    para('Complete nominal statistics were extracted for the three genes: 5,995 NPC2, 5,870 LTBP2 and 5,642 YLPM1 rows. Each count equals the published cis-variant count. Many rows have P &gt; 0.05; no significance filter was applied to the colocalization input. GWAS regional extracts contain 9,241 European and 9,193 combined-study variants before matching.')
    table([['Gene','Full nominal rows','Matched European / combined','Eligible-region matching'],
           *[[g, f"{int(audit.loc[audit.gene==g,'nominal_source_rows'].iloc[0]):,}", ' / '.join(f"{int(coverage.loc[(coverage.gene==g)&(coverage.study==s),'shared_literal_variants'].iloc[0]):,}" for s in STUDIES), f"{coverage.loc[(coverage.gene==g)&(coverage.study==STUDIES[0]),'qtl_matching_fraction'].iloc[0]:.1%}"] for g in GENES]], [55,100,175,140])
    para('<b>Build and variant identity.</b> Explicit coordinates are GRCh38. Positions embedded inside some EyeGEx variant IDs are GRCh37 and were not used. Matching required chromosome, position and the exact unordered allele pair; direct/swapped effects were aligned. No strand rescue, LD proxies, or guessed indel normalization was performed. Exact palindromic matches were retained for the unsigned ABF analysis, flagged for direction ambiguity, and excluded in a separate sensitivity.')
    para('<b>Reconstructed eQTL standard errors.</b> The actual nominal file has beta and P, but no native SE. Approximate SE was reconstructed from the rounded published values as abs(beta) / t.isf(P/2, 404), using nominal dof1 rather than optimized permutation dof2. The authors used a normal-quantile transform of log2 CPM for eQTL mapping; the residual expression SD is unavailable. The working eQTL prior SD of 0.15 is therefore treated as a scale assumption, with 0.075 and 0.30 sensitivities. GWAS SE is the published native value.')
    para('<b>Verification.</b> Eight tests passed: five numerical/scientific tests and three checkpoint tests. All input/output hashes, posterior sums, aligned effects and conditional H4 sets passed. A separate long-double arithmetic calculation agrees with all six primary posteriors to a maximum absolute error of 1.89 × 10⁻¹⁵. A repeat entry-point call reused both stages with unchanged manifest hashes. These tests verify computation; they do not independently validate the reconstructed SE against unavailable native values.')
    para('Implementation: independent Python equations for Wakefield ABFs and the coloc uniform-prior hypothesis sums. The R coloc package was not run or used for an empirical parity comparison. The single-causal-variant model was used; SuSiE-coloc and MR were not run.', 'PilotSmall')
    story.append(PageBreak())

    para('Regional evidence and incomplete coverage', 'PilotTitle')
    story.append(Image(str(out/'Coloc_regional_associations.png'), width=465, height=389))
    para('The strongest source NPC2 eQTL is rs34848569, an indel at GRCh38 chr14:74469519 (P = 1.65534 × 10⁻⁶). That position is absent in both GWAS regional extracts. It contributes no Bayes factor to the fitted shared-variant analysis. No substitute was invented. This limits the completeness of the NPC2 comparison even though approximately 88–89% of eligible regional eQTL rows match the GWAS.')
    para('The European GWAS association at rs754458 has P = 5.939 × 10⁻⁶, versus 3.083 × 10⁻¹⁰ in the combined study. The primary study therefore provides weaker local disease association evidence. The LTBP2 nominal lead has P = 2.11758 × 10⁻⁴ and its source permutation-adjusted gene-level P is approximately 0.209; this is an additional reason to avoid interpreting its near-threshold H4 as confirmed regulation.')
    para('These association plots show source/matched coverage and competing regional patterns. They do not show LD, conditional signals, or cell-specific regulatory effects. The conditional 95% SNP sets in the output tables assume H4 is true; they are not separate fine-mapping evidence or validated causal-variant sets.', 'PilotSmall')
    story.append(PageBreak())

    para('Prior sensitivity and the next defensible step', 'PilotTitle')
    story.append(Image(str(out/'Coloc_prior_sensitivity.png'), width=470, height=205))
    table([['Gene','European H4 range','Combined H4 range','H4 ≥ 0.80 throughout?'],
           *[[g, '–'.join(f'{v:.4g}' for v in primary.loc[(primary.gene==g)&primary.primary_GWAS,['H4_sensitivity_min','H4_sensitivity_max']].iloc[0]), '–'.join(f'{v:.4g}' for v in primary.loc[(primary.gene==g)&~primary.primary_GWAS,['H4_sensitivity_min','H4_sensitivity_max']].iloc[0]), 'No'] for g in GENES]], [65,140,140,125])
    para('Sensitivity design: three genes × two GWAS × two region half-windows (1 Mb and 0.5 Mb) × two allele policies × three eQTL effect-prior SDs × three p12 values = 216 settings. Primary p1 = p2 = 10⁻⁴; primary p12 = 10⁻⁵, varied to 10⁻⁶ and 10⁻⁴. The GWAS prior SD was 0.20 on the log-odds scale. The analysis plan was frozen before posterior estimates, but was not externally registered.')
    para('<b>Scientific contribution.</b> This extension adds a newly fitted, auditable genetic comparison to the existing expression-context pilot. It shows why selecting NPC2 from a cell-expression pattern is insufficient by itself, and why LTBP2 and YLPM1 must remain explicit competing candidates. These are targeted exploratory results from established source cohorts, not independent genetic replication or a genome-wide novel discovery.')
    para('<b>Follow-up.</b> Obtain harmonized, ancestry-matched LD and full multi-signal QTL information for SuSiE-coloc, assess the absent NPC2 indel through an appropriate source dataset, and replicate using independent retinal or relevant cell-type QTL data. Only after adequate instruments and outcome coverage are established would MR be a defensible additional module. A microglia or Müller-glia expression profile alone does not supply cell-specific eQTL evidence.')
    para('Reproducibility and attribution', 'PilotHeading')
    para('Python 3.12.14; NumPy 2.3.5; pandas 2.2.3; SciPy 1.17.0; pysam 0.23.3. Code, filtered input data, stage manifests, tables and reporting scripts accompany this report. New computations and documentation were prepared with OpenAI Codex assistance; original cohort collection and source statistics belong to the cited authors.', 'PilotSmall')
    for label, url in [('Project', 'https://github.com/MDAnusamandiriGroup/glaucoma-cell-context'), ('EyeGEx author data', 'https://doi.org/10.5281/zenodo.18331341'), ('EyeGEx study', 'https://doi.org/10.1038/s41588-019-0351-9'), ('POAG GWAS', 'https://doi.org/10.1038/s41467-020-20851-4'), ('Colocalization method', 'https://doi.org/10.1371/journal.pgen.1004383')]:
        para(escape(label + ': ' + url), 'PilotSmall')

    def footer(canvas, document):
        canvas.saveState()
        canvas.setFont('Pilot', 7)
        canvas.setFillColor(colors.HexColor('#66757e'))
        canvas.drawString(45, 26, 'Glaucoma cell-context | Targeted genetic extension | 6 October 2026')
        canvas.drawRightString(A4[0]-45, 26, str(document.page))
        canvas.restoreState()
    document = SimpleDocTemplate(str(out/'Glaucoma_coloc_pilot_report.pdf'), pagesize=A4,
                                 rightMargin=45, leftMargin=45, topMargin=38, bottomMargin=40,
                                 title='Targeted glaucoma–retina colocalization', author='Mila Desi Anasanti; prepared with OpenAI Codex assistance')
    document.build(story, onFirstPage=footer, onLaterPages=footer)
    (out/'validation.json').write_text(json.dumps(verification, indent=2)+'\n')
    manifest = {'status': 'complete', 'contract': contract,
                'outputs': {str(p.relative_to(out)): hashlib.sha256(p.read_bytes()).hexdigest()
                            for p in sorted(out.iterdir()) if p.is_file()}}
    manifest_path.write_text(json.dumps(manifest, indent=2)+'\n')
    print('COMPLETE report ' + str(out))
    return out


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--index', default=DEFAULT_INDEX)
    args = parser.parse_args()
    build_report(args.index)
