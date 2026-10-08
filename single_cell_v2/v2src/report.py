"""Render scientific figures and a concise report from saved V2 estimates."""
from pathlib import Path
import gzip
import json
import re
import html
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, Table, TableStyle, PageBreak
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.units import inch
from .cache import save_json


CASES = ['COL8A2', 'DGKG', 'NPC2', 'PSMC3', 'RP11-466F5.8', 'SLC2A12']
PROGRAM_LABELS = ['Oxidative phosphorylation', 'Cholesterol homeostasis', 'Glycolysis',
                  'Fatty acid metabolism', 'Unfolded protein response', 'Reactive oxygen species']
PALETTE = {'Rod': '#686dcc', 'Cone': '#d99b30', 'Bipolar': '#49a886', 'RGC': '#397cbd',
    'Muller glia': '#c95881', 'Microglia': '#995eb6', 'Horizontal': '#a98655', '': '#dddddd'}


def figure_save(fig, directory, stem):
    for extension in ['png', 'pdf']:
        target = directory / (stem + '.' + extension)
        temporary = directory / (stem + '.partial.' + extension)
        with temporary.open('xb') as handle:
            fig.savefig(handle, format=extension, dpi=220, bbox_inches='tight', facecolor='white')
            handle.flush(); os.fsync(handle.fileno())
        if temporary.stat().st_size == 0:
            raise ValueError('Empty scientific figure')
        temporary.rename(target)
    plt.close(fig)


def run(root, stages, directory, config):
    figures = directory / 'figures'; figures.mkdir()
    summaries = {name: json.loads((path / 'summary.json').read_text()) for name, path in stages.items()}
    expanded = summaries['expanded']; replication = summaries['replication']; regulatory = summaries['regulatory']; tri = summaries['triangulation']
    evidence = pd.read_csv(stages['triangulation'] / 'candidate_evidence_matrix.csv')
    cases = evidence.set_index('gene').reindex(CASES)
    sampling = pd.read_csv(stages['sampling'] / 'sampling_stability.csv')
    scenarios = ['equal_cells_20', 'equal_cells_10_extended', 'cell_fraction_0.25', 'cell_fraction_0.5',
        'cell_fraction_0.75', 'UMI_fraction_0.25', 'UMI_fraction_0.5', 'UMI_fraction_0.75']
    labels = ['20 cells\n5 types', '10 cells\n7 types', 'Cells\n25%', 'Cells\n50%', 'Cells\n75%', 'UMIs\n25%', 'UMIs\n50%', 'UMIs\n75%']
    heat = sampling.pivot(index='gene', columns='scenario', values='baseline_agreement_fraction').reindex(index=CASES, columns=scenarios)
    fig, ax = plt.subplots(figsize=(10.6, 4.7))
    im = ax.imshow(heat, vmin=0, vmax=1, cmap='YlGnBu', aspect='auto')
    for y in range(len(CASES)):
        for x in range(len(scenarios)):
            v = heat.iloc[y, x]
            ax.text(x, y, f'{v:.0%}', ha='center', va='center', color='white' if v > .65 else '#162a38', fontsize=11)
    ax.set_xticks(range(len(labels)), labels); ax.set_yticks(range(len(CASES)), CASES)
    ax.set_title('Technical agreement with the scenario-specific baseline', fontsize=15, pad=16)
    fig.colorbar(im, ax=ax, fraction=.024, pad=.025).set_label('Fraction of 300 draws')
    fig.text(.02, .01, 'Cell resampling within donor/type; fixed types in every draw. Technical uncertainty, not population confidence intervals.', fontsize=8.8)
    fig.tight_layout(rect=(0,.05,1,1)); figure_save(fig, figures, '01_sampling_stability')

    columns = ['primary_top_type', 'primary_external_top', 'extended_top_type', 'extended_external_top']
    fig, ax = plt.subplots(figsize=(10.6, 4.8)); ax.axis('off')
    contents = [[g] + [str(cases.loc[g,c]) if pd.notna(cases.loc[g,c]) else 'Not mapped' for c in columns] + [str(int(cases.loc[g,'assigned_cell_UMIs']))] for g in CASES]
    table = ax.table(cellText=contents, colLabels=['Candidate', 'V1: 5 types', 'Wang: 5 types', 'V1: 7 types', 'Wang: 7 types', 'V1 UMIs'],
        cellLoc='center', colWidths=[.18,.17,.17,.17,.17,.12], loc='center')
    table.auto_set_font_size(False); table.set_fontsize(10); table.scale(1,2.15)
    for (r,c), cell in table.get_celld().items():
        cell.set_edgecolor('white')
        if r == 0:cell.set_facecolor('#152f47'); cell.set_text_props(color='white', weight='bold')
        elif c in [1,2,3,4]:
            value = contents[r-1][c]; cell.set_facecolor(PALETTE.get(value,'#efefef')); cell.set_text_props(color='white' if value in PALETTE and value else '#333333')
        else:cell.set_facecolor('#eff3f7')
    ax.set_title('Cross-study retinal RNA contexts: six prespecified cases',fontsize=15,pad=12)
    fig.text(.02,.015,'V1: 19,865 whole-cell profiles / 3 donors. Wang: 51,645 nuclear profiles / 4 donors. UMIs are V1 assigned-cell gene counts.',fontsize=8.7)
    figure_save(fig,figures,'02_cross_study_contexts')

    programs = pd.read_csv(stages['triangulation'] / 'cross_study_program_context.csv')
    types = config['primary_types']; matrices=[]
    for study in ['Lukowski2019','Wang2022']:
        f=programs.loc[programs.study.eq(study)&programs.comparison.eq('primary')].set_index('program').reindex(config['programs'])
        matrices.append(f[['score_'+ct for ct in types]].to_numpy())
    bound=max(float(np.max(np.abs(x))) for x in matrices)
    fig,axes=plt.subplots(1,2,figsize=(11.2,5.2),sharey=True)
    for ax,matrix,title in zip(axes,matrices,['V1: whole-cell RNA, 3 donors','Wang: nuclear RNA, 4 donors']):
        im=ax.imshow(matrix,vmin=-bound,vmax=bound,cmap='RdBu_r',aspect='auto')
        ax.set_xticks(range(len(types)),['Rod','Cone','Bipolar','RGC','Müller glia'],rotation=35,ha='right')
        ax.set_yticks(range(6),PROGRAM_LABELS); ax.set_title(title,fontsize=12)
    fig.subplots_adjust(left=.24,right=.9,wspace=.12,bottom=.23,top=.84)
    cb=fig.add_axes([.92,.27,.018,.53]);fig.colorbar(im,cax=cb).set_label('Matched RNA contrast')
    fig.suptitle('Identical shared panel/control genes across studies',fontsize=15)
    fig.text(.02,.03,'Equal-donor pseudobulk log-CPM contrasts. Relative transcriptional programs; no metabolic flux or disease effect estimates.',fontsize=8.6)
    figure_save(fig,figures,'03_shared_program_scores')

    labels=['Published nominations','V1 exact-symbol mapped','RNA context evaluable in both studies','Same top type in both studies','Descriptive context screening flag']
    values=[228,163,replication['primary_evaluable'],replication['primary_context_concordant'],tri['descriptive_primary_context_flag_genes']]
    fig,ax=plt.subplots(figsize=(9.2,4.4)); y=np.arange(len(values))
    ax.barh(y,values,color=['#71879a','#4c819d','#397cbd','#279591','#995eb6'],height=.62)
    ax.set_yticks(y,labels);ax.invert_yaxis();ax.set_xlim(0,255);ax.set_xlabel('Genes')
    for i,v in enumerate(values):ax.text(v+3,i,str(v),va='center',weight='bold')
    ax.spines[['top','right']].set_visible(False);ax.set_title('Candidate screening flow: fixed five-type comparison',pad=14,fontsize=14)
    fig.text(.01,.015,'The screening flag is a descriptive rule; it is not statistical significance, a causal posterior, or a treatment recommendation.',fontsize=8)
    fig.tight_layout(rect=(0,.07,1,1));figure_save(fig,figures,'04_candidate_screening_flow')

    coordinates=[]
    wanted={'NPC2','LTBP2','YLPM1'}
    with gzip.open(root/'data/raw/gencode.v32.annotation.gtf.gz','rt') as handle:
        for line in handle:
            if line.startswith('#'):continue
            fields=line.rstrip().split('\t')
            if fields[2]!='gene':continue
            match=re.search(r'gene_name "([^"]+)"',fields[8])
            if match and match.group(1) in wanted:
                coordinates.append({'gene':match.group(1),'chromosome':fields[0],
                    'start_0based':int(fields[3])-1,'end_exclusive':int(fields[4]),'strand':fields[6],
                    'source':'GENCODE v32 GRCh38; gene body, not enhancer assignment'})
    if len(coordinates)!=3 or any(x['chromosome']!='chr14' for x in coordinates):raise ValueError('Case gene coordinates are ambiguous')
    pd.DataFrame(coordinates).to_csv(directory/'NPC2_case_gene_coordinates.csv',index=False)
    v1=pd.read_csv(stages['expanded']/'donor_candidate_expression.csv')
    ext=pd.read_csv(stages['replication']/'external_donor_candidate_expression.csv')
    fig,axes=plt.subplots(2,1,figsize=(10.6,7.8),gridspec_kw={'height_ratios':[1,1.55]})
    pos=74618126
    for y,record in enumerate(sorted(coordinates,key=lambda x:x['start_0based'])):
        color='#995eb6' if record['gene']=='YLPM1' else '#397cbd'
        axes[0].plot([record['start_0based']/1e6,record['end_exclusive']/1e6],[y,y],lw=7,color=color,solid_capstyle='butt')
        axes[0].text(record['start_0based']/1e6,y+.16,record['gene'],fontsize=10,weight='bold')
    axes[0].axvline(pos/1e6,color='#bf4f46',ls='--',label='rs754458')
    axes[0].set_ylim(-.4,2.55);axes[0].set_yticks([]);axes[0].set_xlabel('Chromosome 14 position (Mb), GRCh38')
    axes[0].set_title('NPC2 locus: RNA context is reproducible; gene-target assignment remains open',fontsize=13,pad=13)
    axes[0].legend(frameon=False,fontsize=9,loc='upper left');axes[0].spines[['top','right','left']].set_visible(False)
    extended=config['extended_types'];x=np.arange(len(extended))
    for frame,offset,color,label in [(v1,-.16,'#397cbd','V1: 3 donors'),(ext,.16,'#c95881','Wang: 4 donors')]:
        selected=frame.loc[frame.gene.eq('NPC2')&frame.cell_type.isin(extended)].copy()
        if 'qc' in selected:selected=selected.loc[selected.qc.eq('primary')]
        mean=selected.groupby('cell_type').log1p_pseudobulk_CPM.mean().reindex(extended)
        axes[1].bar(x+offset,mean,width=.28,color=color,alpha=.55,label=label)
        for t,ct in enumerate(extended):
            vals=selected.loc[selected.cell_type.eq(ct),'log1p_pseudobulk_CPM'].to_numpy()
            jitter=np.linspace(-.06,.06,len(vals));axes[1].scatter(t+offset+jitter,vals,s=24,color=color,zorder=3)
    axes[1].set_xticks(x,['Rod','Cone','Bipolar','Horizontal','RGC','Müller glia','Microglia'],rotation=25,ha='right')
    axes[1].set_ylabel('Donor log1p pseudobulk CPM');axes[1].set_title('Seven-type comparison; each dot is one biological donor',fontsize=11)
    axes[1].legend(frameon=False,fontsize=9);axes[1].spines[['top','right']].set_visible(False)
    fig.text(.02,.016,'rs754458: microglia ATAC overlap; published retinal eQTL targets NPC2/LTBP2; HiChIP target YLPM1; BPNet high-effect call = 0.\nRetinal eQTL resources overlap across publications; these are not independent genetic replications or proven causal targets.',fontsize=8.6)
    fig.tight_layout(rect=(0,.105,1,1));figure_save(fig,figures,'05_NPC2_locus_evidence')

    key={'V1_cells':expanded['cells'],'V1_donors':3,'external_profiles':replication['external_cells'],'external_donors':4,
        'candidates':228,'V1_mapped':163,'primary_evaluable':replication['primary_evaluable'],
        'primary_context_concordant':replication['primary_context_concordant'],
        'technical_draws':2400,'screening_flag_genes':tri['descriptive_primary_context_flag_genes'],
        'predefined_programs':6,'program_top_concordance':tri['primary_program_top_concordance'],
        'regulatory_gene_links':regulatory['candidate_genes_with_any_published_gene_link']}
    save_json(directory/'report_summary.json',key)
    rows=[]
    for gene in CASES:
        r=cases.loc[gene]
        rows.append('| '+gene+' | '+str(r.primary_top_type)+' / '+str(r.primary_external_top)+' | '+str(r.extended_top_type)+' / '+str(r.extended_external_top)+' | '+str(int(r.assigned_cell_UMIs))+' |')
    markdown=f'''# Glaucoma cell-context V2: executed Python pilot

V2 expands the frozen published-genetics pilot to all 228 nominations, technical sampling robustness, independent-study RNA context, predefined metabolic programs and regulatory triangulation. It uses 19,865 V1 whole-cell profiles from three healthy donors and 51,645 nuclear profiles from four different published donors. The two studies are analysed separately.

## Executed results

- 163/228 candidates mapped exactly in V1; 156/228 in the original Wang count distribution.
- 94/163 V1 candidates retained the same unique top type in all three leave-one-donor-out folds; 150 retained their top type under strict QC.
- 2,400 technical draws: 300 per scenario across eight frozen scenarios. Sampling stays within donor/type; comparison types are fixed within each scenario.
- {key['primary_context_concordant']}/{key['primary_evaluable']} jointly evaluable candidates had the same highest-expression type in the fixed five-type comparison. With seven types, 70/140 were concordant.
- {key['screening_flag_genes']} genes passed the descriptive context screening rule. This is not significance or a causal probability.
- Six predefined programs were scored with expression-matched controls. Only {key['program_top_concordance']}/6 had the same relative top type across studies when using identical shared panel/control pairs.
- 25 index variants overlapped reproducible Wang ATAC peaks; 50 were represented in the published glaucoma annotation table; 14 candidates had a source-reported gene link. Computed cell-type peak intersections reproduced all 50 published annotations exactly.
- Four index variants overlapped the recent Yuan 2026 source-listed retina/macula cCRE table. All four overlaps survived plausible one-base endpoint conventions. Microglia were excluded from this supplied atlas table because of low cell numbers.

## Six prespecified cases

| Gene | Five types: V1 / Wang | Seven types: V1 / Wang | V1 assigned-cell UMIs |
|---|---|---|---:|
{chr(10).join(rows)}

NPC2 is technically stable in both comparison sets and concordant across studies. Its five-type maximum is Müller glia; its seven-type maximum is microglia. The same index variant, rs754458, overlaps a microglia ATAC peak. However, the published retinal eQTL table names NPC2 and LTBP2, the HiChIP table links YLPM1, and BPNet does not call this variant high-effect. This is a locus with reproducible RNA context and unresolved target-gene mechanism. It is not a new discovery of NPC2, a microglia-specific eQTL, or a proven causal NPC2 effect.

DGKG changes context across studies. SLC2A12 agrees across RNA studies but has low V1 counts and less than 90% agreement under equal-cell sampling. PSMC3 agrees across studies but is unstable under V1 equal-cell sampling. COL8A2 is especially sparse. RP11-466F5.8 is not exactly symbol-mapped in the external distribution; no alias was inferred.

## Methods and limits

Raw UMIs are summed within donor and published cell type, divided by all-gene UMIs in that group, converted to log1p CPM and averaged with equal donor weights. Primary types are Rod, Cone, Bipolar, RGC and Müller glia (at least 20 cells per donor/type). The extended sensitivity adds Horizontal and Microglia (at least 10 cells per V1 donor/type). Donor omissions reaggregate saved counts; they do not refit clustering.

Technical draws use cell sampling without replacement or binomial thinning of selected-gene and remaining UMIs together. They describe uncertainty conditional on the observed samples, not biological population confidence intervals. The 90% screening cutoff is descriptive.

Six MSigDB Hallmark 2025.1.Hs programs were selected before module execution. Control genes came from 20 expression bins, with up to 20 controls per panel gene and all six panels excluded from controls. Cross-study scoring uses identical intersected panel/control pairs and weights. The scores describe RNA programs, not metabolic flux, mitochondrial function or disease-specific activation. Whole-cell versus nuclear capture, donor characteristics and laboratory effects can contribute to differences.

Regulatory evidence uses exact index rsIDs and GRCh38 coordinates. Nearest-gene marker peaks, target predictions, promoter co-accessibility, HiChIP and retinal eQTL are separate evidence categories. Published BPNet scores and negative calls are imported; no model was trained. Hamel and Wang reuse overlapping retinal eQTL resources, so this is not independent genetic replication. No new GWAS, fine-mapping, colocalization or Mendelian randomization has been fitted.

Yuan's separate peak/gene matrices were not paired by row order without explicit identifiers. A validated enhancer-gene edge list is required for that additional ABC analysis. New colocalization/MR would require complete harmonized GWAS/QTL summary statistics; thresholded published nomination tables are insufficient.

## Reproducibility and disclosure

Python stages are checkpointed with source, software and code hashes. V1 models/checkpoints are read-only inputs. Failed attempts remain diagnostic and are excluded from handover. Code implementation and reporting were AI-assisted; executed outputs and numerical checks are supplied. This is a computational pilot for review, not a manuscript accepted for publication or an experimentally validated mechanism.

## Sources

- Lukowski et al. 2019, DOI https://doi.org/10.15252/embj.2018100811; counts https://zenodo.org/records/5515631
- Hamel et al. 2024, DOI https://doi.org/10.1038/s41467-023-44380-y
- Wang et al. 2022, DOI https://doi.org/10.1016/j.xgen.2022.100164; GSE196235 and CellxGene collection https://cellxgene.cziscience.com/collections/348da6dc-5bf6-435d-adc5-37747b9ae38a
- Yuan et al. 2026, DOI https://doi.org/10.1126/sciadv.adv9162
- MSigDB Hallmark 2025.1.Hs: https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp
- MitoCarta3.0: https://www.broadinstitute.org/mitocarta/mitocarta30-inventory-mammalian-mitochondrial-proteins-and-pathways
- GENCODE v32 GRCh38: https://www.gencodegenes.org/human/release_32.html
'''
    (directory/'REPORT.md').write_text(markdown,encoding='utf-8')
    styles=getSampleStyleSheet();styles['BodyText'].fontSize=9.2;styles['BodyText'].leading=12.5
    styles['Heading1'].textColor=colors.HexColor('#152f47');styles['Heading2'].textColor=colors.HexColor('#152f47')
    story=[]
    def paragraph(text,style='BodyText'):story.append(Paragraph(html.escape(text),styles[style]));story.append(Spacer(1,.085*inch))
    def picture(stem,width=7.0*inch):
        from PIL import Image as PILImage
        image_path=figures/(stem+'.png')
        with PILImage.open(image_path) as im:w,h=im.size
        story.append(Image(str(image_path),width=width,height=width*h/w));story.append(Spacer(1,.08*inch))
    paragraph('Glaucoma cell-context V2','Title');paragraph('Executed Python pilot | 6 October 2026','Heading2')
    paragraph('228 published gene nominations; sampling robustness; cross-study retinal RNA; predefined metabolic programs; regulatory triangulation.')
    metrics=[['Evidence layer','Executed result'],['RNA datasets','19,865 profiles / 3 donors; 51,645 profiles / 4 donors'],['Mapping','163 V1 genes; 156 external genes; 65 V1 symbols missing'],['RNA context agreement',f"{key['primary_context_concordant']}/{key['primary_evaluable']} in five fixed types; 70/140 in seven"],['Technical sampling','2,400 draws across 8 frozen scenarios'],['Descriptive context flag',f"{key['screening_flag_genes']} genes; not a causal probability"],['Matched program agreement',f"{key['program_top_concordance']}/6 programs share the same relative top type"],['Regulatory layers','25 index ATAC overlaps; 14 source gene-linked candidates']]
    t=Table(metrics,colWidths=[1.7*inch,5.3*inch]);t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#152f47')),('TEXTCOLOR',(0,0),(-1,0),colors.white),('FONTNAME',(0,0),(-1,0),'Helvetica-Bold'),('FONTSIZE',(0,0),(-1,-1),9),('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.HexColor('#f0f4f7'),colors.white]),('BOTTOMPADDING',(0,0),(-1,-1),9),('TOPPADDING',(0,0),(-1,-1),9)]));story.append(t)
    paragraph('NPC2 has reproducible RNA context and a microglia-accessible index variant, with unresolved competing gene links. DGKG does not reproduce its top context; sparse genes remain weak even when a donor-omission result appears stable.')
    paragraph('AI-assisted computational reanalysis of published healthy-donor resources. No new disease effect, causal target, fine-mapping, colocalization or MR estimate is claimed.')
    story.append(PageBreak());paragraph('Sampling stability','Heading1');picture('01_sampling_stability')
    paragraph('Each draw retains the same type set as its own baseline. Primary sampling uses five eligible classes; a separate equal-10-cell sensitivity uses seven. Counts are drawn within donor/type, and all donors have equal weight. Depth thinning conserves selected-gene plus other UMIs. Percentages are technical agreement frequencies, not population confidence intervals.')
    paragraph('NPC2 agrees in every draw of these scenarios. COL8A2 has only six assigned-cell UMIs; many equal-20-cell draws contain no detected UMIs. PSMC3 remains unstable under equal-cell sampling despite cross-study top-type agreement.')
    story.append(PageBreak());paragraph('Cross-study RNA context','Heading1');picture('02_cross_study_contexts');picture('04_candidate_screening_flow',width=6.7*inch)
    paragraph('Wang nuclei retain published cell labels and four donor identities. All 51,645 GEO barcodes exactly matched the curated CellxGene metadata; five validation genes matched their UMI counts cell by cell across both distributions. Eight eyes are nested within four donors. These are separate whole-cell and nuclear RNA studies, not a pooled disease comparison.')
    story.append(PageBreak());paragraph('Predefined metabolic RNA programs','Heading1');picture('03_shared_program_scores')
    paragraph('Both studies use the same intersected panel/control pairs and gene weights, with donor means computed separately. Only three of six programs reproduce their relative top type among five classes. The scores are transcriptional contrasts and do not measure metabolic flux. Capture chemistry, nucleus versus whole-cell sampling, donor characteristics and laboratory effects remain potential explanations for differences.')
    paragraph('No pathway enrichment p-values, population tests, or aging/disease models were calculated from these small healthy-donor samples.')
    story.append(PageBreak());paragraph('NPC2: context and competing targets','Heading1');picture('05_NPC2_locus_evidence',width=6.7*inch)
    paragraph('The rs754458 index variant intersects a microglia ATAC peak. Published retinal eQTL annotations include NPC2 and LTBP2; the published HiChIP annotation points to YLPM1. BPNet reports no high-effect call (source FDR approximately 0.882). RNA and regulatory context motivate a hypothesis; they do not resolve the causal gene or establish a microglia-specific eQTL.')
    story.append(PageBreak());paragraph('Scientific scope and remaining inputs','Heading1')
    for text in ['Nearest-gene marker-peak annotations are weaker than explicit gene links. The evidence matrix keeps predicted target, promoter co-accessibility, HiChIP and retinal eQTL categories separate.',
        'Hamel and Wang reuse overlapping published retinal eQTL resources. The independent-study RNA comparison does not make the shared genetic layer an independent replication.',
        'Four index variants intersect source-listed cCREs from Yuan 2026, with boundary-robust overlaps. That supplied retina/macula table excludes microglia because of low cell numbers; absence there is unassessed microglia coverage.',
        'Separate Yuan peak and gene matrices were not joined by row order without explicit pairing identifiers. That optional ABC extension needs a validated enhancer-gene edge list.',
        'New colocalization or MR requires complete GWAS/QTL summary statistics, verified allele/build handling, LD references where needed, and adequate instrument strength. None can be reconstructed faithfully from thresholded nomination tables.',
        'The pilot is suitable as reproducible computational evidence for an application. Publication novelty needs further benchmarking, complete genetic analyses where available, and biological validation; Q1 acceptance and citation counts are not predictable from this pilot.']:
        paragraph(text)
    paragraph('Reproducibility','Heading2');paragraph('Each stage validates source, code, software and output hashes. Model-independent plotting uses saved estimates. V1 checkpoints are preserved; failed and temporary files are excluded from the handover. The supplied code is AI-assisted and the executed results remain subject to scientific review.')
    paragraph('Sources','Heading2')
    for line in markdown.split('## Sources\n\n')[1].splitlines():
        if line.startswith('- '):paragraph(line[2:])
    temporary_pdf = directory / 'Glaucoma_V2_Report.partial.pdf'
    doc=SimpleDocTemplate(str(temporary_pdf),pagesize=(8.27*inch,11.69*inch),rightMargin=.6*inch,leftMargin=.6*inch,topMargin=.52*inch,bottomMargin=.55*inch)
    def footer(canvas,doc):canvas.setFont('Helvetica',8);canvas.setFillColor(colors.grey);canvas.drawString(.6*inch,.29*inch,'Glaucoma cell-context V2 | executed computational pilot');canvas.drawRightString(7.67*inch,.29*inch,str(doc.page))
    doc.build(story,onFirstPage=footer,onLaterPages=footer)
    temporary_pdf.rename(directory / 'Glaucoma_V2_Report.pdf')
