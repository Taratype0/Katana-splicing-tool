"""BAM/index transport checks; native RSeQC is exercised in separate acceptance."""
import json
from pathlib import Path
import sys

import pandas as pd
import pytest

from katana_next.flow import stages


def test_owned_adjacent_index_is_preserved_when_bam_target_lacks_one(tmp_path, monkeypatch):
    from katana_next.sashimi_fixture import write_indexed_bam
    source=tmp_path/'original';source.mkdir()
    bam=source/'input.bam';write_indexed_bam(bam)
    index=tmp_path/'index-task';index.mkdir()
    source_index=Path(str(bam)+'.bai')
    separate=index/'input.bam.bai';source_index.rename(separate)
    owned=tmp_path/'owned';owned.mkdir()
    staged=owned/'input.bam'
    try:
        staged.symlink_to(bam)
        Path(str(staged)+'.bai').symlink_to(separate)
    except OSError as error:
        pytest.skip('Symlink creation is unavailable on this host: '+str(error))
    assert not Path(str(bam)+'.bai').exists()
    bed=tmp_path/'reference.bed';bed.write_text('chr1\t0\t100\tt1\t0\t+\t0\t100\t0\t1\t100\t0\n')
    def native_transport_fixture(argv,out,name):
        assert Path(argv[argv.index('-i')+1])==staged.absolute()
        assert Path(argv[argv.index('-i')+1]+'.bai').is_file()
        (out/'gene_body.geneBodyCoverage.txt').write_text('Percentile\t'+'\t'.join(map(str,range(1,101)))+'\ninput\t'+'\t'.join(['3']*100)+'\n')
        (out/(name+'.command.json')).write_text(json.dumps({'fixture':'transport only','argv':argv}))
    monkeypatch.setattr(stages,'_run',native_transport_fixture)
    result=stages.execute_stage({'adapter':'qc.genebody','inputs':{'bam':str(staged),'reference_bed':str(bed)},
        'parameters':{'sample_id':'sample1','rseqc_executable':sys.executable}},tmp_path/'result')
    assert result['evidence']['bam']==str(staged.absolute())
    assert result['evidence']['bam_resolved_source']==str(bam.resolve())
    assert len(pd.read_csv(tmp_path/'result/gene_body_coverage.tsv',sep='\t'))==100
    assert not Path(str(bam)+'.bai').exists()
