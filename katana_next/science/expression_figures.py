from pathlib import Path
import numpy as np
import pandas as pd


def expression_figures(directory, labels):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    directory=Path(directory)
    names={
        'A_at_Breference':f"{labels['target_a']} vs {labels['reference_a']} | {labels['reference_b']}",
        'A_at_Btarget':f"{labels['target_a']} vs {labels['reference_a']} | {labels['target_b']}",
        'B_at_Areference':f"{labels['target_b']} vs {labels['reference_b']} | {labels['reference_a']}",
        'B_at_Atarget':f"{labels['target_b']} vs {labels['reference_b']} | {labels['target_a']}",
        'interaction':'Genotype-by-treatment interaction'}
    records=[]
    for name,label in names.items():
        x=pd.read_csv(directory/(name+'.tsv'),sep='\t')
        valid=x.padj.notna()
        records.append(dict(contrast=name,label=label,tested_genes=int(valid.sum()),FDR_pass=int(x.padj.le(.05).sum()),FDR_and_abs_log2FC1=int((x.padj.le(.05)&x.log2FoldChange.abs().ge(1)).sum())))
    summary=pd.DataFrame(records);summary.to_csv(directory/'expression_test_summary.tsv',sep='\t',index=False)
    fig,ax=plt.subplots(figsize=(11,5))
    ax.barh(summary.label,summary.FDR_pass,color='#2166ac',label='FDR <= 0.05')
    ax.barh(summary.label,summary.FDR_and_abs_log2FC1,color='#b2182b',label='Also |log2 effect| >= 1')
    for i,row in summary.iterrows():ax.text(row.FDR_pass,i,f' {row.FDR_pass} / {row.tested_genes} tested',va='center',fontsize=8)
    ax.set_xlabel('Genes; interaction coefficient is a difference of log2 effects')
    ax.set_title('Expression responses: simple effects and formal interaction');ax.legend(loc='upper center',bbox_to_anchor=(.5,-.15),ncol=2)
    fig.tight_layout()
    paths=[]
    for ext in ('png','svg','pdf'):
        p=directory/('expression_effects.'+ext);fig.savefig(p,dpi=300,bbox_inches='tight');paths.append(p)
    plt.close(fig)
    pathways=pd.read_csv(directory/'ranked_pathways.tsv',sep='\t')
    chosen=pathways.sort_values('study_family_FDR').drop_duplicates('pathway').head(16).pathway
    matrix=pathways[pathways.pathway.isin(chosen)].pivot(index='pathway',columns='contrast',values='NES').reindex(index=chosen,columns=names)
    fdr=pathways[pathways.pathway.isin(chosen)].pivot(index='pathway',columns='contrast',values='study_family_FDR').reindex_like(matrix)
    fig,ax=plt.subplots(figsize=(13,max(5,len(matrix)*.32)))
    if len(matrix):
        bound=max(1,np.nanmax(np.abs(matrix.to_numpy())))
        im=ax.imshow(matrix,aspect='auto',cmap='RdBu_r',vmin=-bound,vmax=bound)
        ax.set_yticks(range(len(matrix)),[s.split(' / ',1)[-1] for s in matrix.index],fontsize=8)
        ax.set_xticks(range(5),[names[s] for s in matrix.columns],rotation=45,ha='right',fontsize=8)
        for i in range(len(matrix)):
            for j in range(5):
                if pd.notna(fdr.iloc[i,j]) and fdr.iloc[i,j]<=.05:ax.text(j,i,'*',ha='center',va='center')
        fig.colorbar(im,ax=ax,label='Normalized enrichment score')
    ax.set_title('Ranked expression pathways: top 16 by adjusted P-value\n* Within-design family FDR <= 0.05; unstarred cells are not significant')
    fig.tight_layout()
    for ext in ('png','svg','pdf'):
        p=directory/('ranked_pathways.'+ext);fig.savefig(p,dpi=300,bbox_inches='tight');paths.append(p)
    plt.close(fig)
    return paths
