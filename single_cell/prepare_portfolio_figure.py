#!/usr/bin/env python3
"""Render a compact portfolio view from saved UMAP coordinates; no refitting."""
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scsrc.common import publish
from scsrc.reporting import COLORS


def main():
    root = Path(__file__).resolve().parent
    path = root / 'published/figures/02_retina_umap_portfolio.png'
    if path.exists():
        print('REUSE portfolio figure: ' + str(path))
        return
    frame = pd.read_csv(root / 'published/tables/cell_embeddings.csv')
    fig, ax = plt.subplots(figsize=(7.5, 3.71))
    fig.subplots_adjust(left=0.07, right=0.74, top=0.85, bottom=0.15)
    for cell_type in COLORS:
        mask = frame['published_cell_type'].eq(cell_type)
        ax.scatter(frame.loc[mask,'UMAP1'],frame.loc[mask,'UMAP2'],s=2.0,color=COLORS[cell_type],alpha=0.65,linewidths=0,label=cell_type)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel('UMAP 1',fontsize=10); ax.set_ylabel('UMAP 2',fontsize=10)
    ax.spines[['top','right']].set_visible(False)
    legend=ax.legend(loc='center left',bbox_to_anchor=(1.01,0.5),title='Published labels',fontsize=9.5,title_fontsize=9.5,frameon=False,markerscale=4,handletextpad=0.4)
    for handle in legend.legend_handles: handle.set_alpha(1)
    fig.suptitle(f'Count-level retina pilot: {len(frame):,} cells | 3 donors',fontsize=12,weight='semibold')
    fig.text(0.07,0.02,'UMAP computed from filtered UMI counts; published cell labels retained.',fontsize=9,color='#555555')
    publish(path,lambda p:fig.savefig(p,dpi=200,facecolor='white'))
    plt.close(fig)
    print(str(path))


if __name__ == '__main__':
    main()
