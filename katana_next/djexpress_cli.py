"""Run source-pinned DJExpress through Katana, retaining both author passes."""
import argparse
from pathlib import Path
import pandas as pd
from katana_next.science.djexpress import run_djexpress
def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['counts','traits','out']:p.add_argument('--'+name,required=True)
    for name in ['rscript','gene-expression','gtf']:p.add_argument('--'+name)
    p.add_argument('--input-kind',choices=['junction_counts','rmats_form_evidence'],default='junction_counts');a=p.parse_args()
    read=lambda path:pd.read_csv(path,sep='\t',index_col=0)
    out=run_djexpress(read(a.counts),read(a.traits),a.out,rscript=a.rscript,gene_expression=read(a.gene_expression) if a.gene_expression else None,gtf=a.gtf,input_kind=a.input_kind)
    print(out)
if __name__=='__main__':main()
