"""Reference-style GO/KEGG facets from saved statistics; no enrichment fitting."""
from pathlib import Path
import textwrap
import numpy as np
import pandas as pd

FACETS = ('BP', 'CC', 'MF', 'KEGG')


def render_enrichment_facets(frame, destination, *, title='', top_n=5):
    """Plot recorded ratios, counts and adjusted P; preserve the complete source.

    Required columns: category, term, Count, GeneRatio, p.adjust. Categories must
    be explicitly supplied, never guessed from GO descriptions. Selection is
    display-only: adjusted P <= .05, then top_n per facet, ties by term name.
    """
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.colors import LinearSegmentedColormap, Normalize
    from matplotlib.patches import Rectangle
    from matplotlib.ticker import MaxNLocator, FormatStrFormatter
    required={'category','term','Count','GeneRatio','p.adjust'}
    if not required.issubset(frame):
        raise ValueError('Facet figure requires explicit category, term, Count, GeneRatio and p.adjust')
    data=frame.copy(deep=True)
    for key in ('Count','GeneRatio','p.adjust'):
        data[key]=pd.to_numeric(data[key],errors='raise')
    if not data.category.isin(FACETS).all():
        raise ValueError('GO BP/CC/MF categories must be explicitly recorded; unknown categories cannot be inferred')
    if not np.isfinite(data[['Count','GeneRatio','p.adjust']].to_numpy()).all():
        raise ValueError('Non-finite enrichment plotting values')
    if ((data['p.adjust']<0)|(data['p.adjust']>1)|(data.GeneRatio<0)|(data.GeneRatio>1)|(data.Count<0)).any():
        raise ValueError('Invalid enrichment ratio, count or adjusted P')
    selected=[]
    for category in FACETS:
        part=data.loc[(data.category==category)&(data['p.adjust']<=.05)&(data.Count>0)]
        selected.append(part.sort_values(['p.adjust','term'],kind='stable').head(top_n))
    plotted=pd.concat(selected,ignore_index=True)
    counts=[max(1,len(p)) for p in selected]
    fig=Figure(figsize=(8.6,max(6.5,sum(counts)*.34+1.6)))
    FigureCanvasAgg(fig)
    grid=fig.add_gridspec(4,1,left=.40,right=.77,bottom=.10,top=.91,hspace=.10,height_ratios=counts)
    qmin=float(plotted['p.adjust'].min()) if len(plotted) else 0.
    qmax=float(plotted['p.adjust'].max()) if len(plotted) else .05
    norm=Normalize(qmin,qmax if qmax>qmin else qmin+max(1e-5,qmin*.1))
    cmap=LinearSegmentedColormap.from_list('katana_adjusted_p',['#c63f4a','#2863a5'])
    max_count=max(1,float(plotted.Count.max())) if len(plotted) else 1
    size=lambda values: np.asarray(values)/max_count*180
    xmax=max(.01,float(plotted.GeneRatio.max())*1.12) if len(plotted) else .1
    first=None
    for i,(category,part) in enumerate(zip(FACETS,selected)):
        ax=fig.add_subplot(grid[i],sharex=first)
        if first is None:first=ax
        ax.set_axisbelow(True);ax.grid(color='#e1e1e1',linewidth=.7)
        for spine in ax.spines.values():spine.set_color('#777');spine.set_linewidth(.7)
        ax.set_xlim(0,xmax);ax.xaxis.set_major_locator(MaxNLocator(5));ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        ax.tick_params(axis='both',length=0,labelsize=9,pad=4)
        if len(part):
            y=np.arange(len(part))
            ax.scatter(part.GeneRatio,y,s=size(part.Count),c=part['p.adjust'],norm=norm,cmap=cmap,linewidths=.25,edgecolors='#555',zorder=3)
            labels=['\n'.join(textwrap.wrap(str(t),32,break_long_words=False)) for t in part.term]
            ax.set_yticks(y,labels,fontfamily='DejaVu Serif',fontsize=9)
            ax.set_ylim(len(part)-.5,-.5)
        else:
            ax.set_yticks([]);ax.set_ylim(.5,-.5)
            ax.text(.5,.5,'No significant terms',transform=ax.transAxes,ha='center',va='center',fontsize=9,color='#666')
        ax.add_patch(Rectangle((1,0),.065,1,transform=ax.transAxes,clip_on=False,facecolor='#d6d6d6',edgecolor='#777',linewidth=.7))
        ax.text(1.033,.5,category,rotation=270,ha='center',va='center',transform=ax.transAxes,fontsize=10,fontfamily='DejaVu Serif')
        if i<3:ax.tick_params(labelbottom=False)
        else:ax.set_xlabel('Gene Ratio',fontfamily='DejaVu Serif',fontsize=11)
    fig.suptitle(title,fontsize=12,y=.975)
    if len(plotted):
        from matplotlib.cm import ScalarMappable
        color_ax=fig.add_axes([.865,.55,.022,.19])
        cb=fig.colorbar(ScalarMappable(norm=norm,cmap=cmap),cax=color_ax)
        color_ax.set_title('p.adjust',fontsize=9,loc='left',pad=9)
        cb.ax.tick_params(labelsize=8,length=2);cb.locator=MaxNLocator(4);cb.update_ticks()
        legend_counts=sorted(set([int(plotted.Count.min()),int(round((plotted.Count.min()+max_count)/2)),int(max_count)]))
        handles=[first.scatter([],[],s=size(n),facecolors='none',edgecolors='#333',linewidths=.6) for n in legend_counts]
        fig.legend(handles,[str(n) for n in legend_counts],title='Counts',loc='upper left',bbox_to_anchor=(.825,.46),frameon=False,fontsize=9,title_fontsize=9,labelspacing=1)
    destination=Path(destination);destination.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(destination,dpi=180,facecolor='white')
    fig.savefig(destination.with_suffix('.svg'),facecolor='white')
    fig.clear()
    return plotted
