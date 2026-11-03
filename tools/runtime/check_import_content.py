"""Installed-runtime regression for completed-result import; generated data only."""
from pathlib import Path
import argparse
import csv
import hashlib
import io
import json
import os
import subprocess
from unittest.mock import patch


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--outdir', type=Path, required=True)
    root = parser.parse_args().outdir.resolve()
    root.mkdir(parents=True, exist_ok=False)
    from katana_next.flow import file_io
    checks = []
    source = root / 'stable.tsv'
    source.write_bytes(b'original')
    stamp = file_io.identity(source)
    phase = [0]
    reads = []
    native_digest = file_io.hashlib.file_digest
    def identity(path):
        return [*stamp[:4], stamp[4] + phase[0]]
    def metadata_only(handle, algorithm):
        result = native_digest(handle, algorithm)
        reads.append(True)
        phase[0] = 1
        return result
    with patch.dict(os.environ, {'KATANA_VERIFY_CONTENT': '1'}), patch.object(file_io, 'identity', identity), patch.object(file_io.hashlib, 'file_digest', metadata_only):
        assert file_io.digest(source) == hashlib.sha256(b'original').hexdigest()
        assert len(reads) == 2
    checks.append('metadata transition requires identical stable content')
    phase[0] = 0
    def changed(handle, algorithm):
        result = native_digest(handle, algorithm)
        if phase[0] == 0:
            source.write_bytes(b'modified')
            phase[0] = 1
        return result
    with patch.dict(os.environ, {'KATANA_VERIFY_CONTENT': '1'}), patch.object(file_io, 'identity', identity), patch.object(file_io.hashlib, 'file_digest', changed):
        try:
            file_io.digest(source)
        except ValueError as error:
            assert 'File changed during content verification' in str(error)
        else:
            raise AssertionError('Changed scientific input was accepted')
    checks.append('same-size changed bytes rejected')
    generated = root / 'generated.tsv'
    with patch.object(file_io.hashlib, 'file_digest', side_effect=AssertionError('unnecessary projection reread')):
        assert file_io.write_verified(generated, b'G1\tE1\t1e-310\n') == hashlib.sha256(generated.read_bytes()).hexdigest()
    checks.append('generated projection checksum uses write pass')
    bundle = root / 'upstream'; bundle.mkdir()
    native = bundle / 'dexseq_exon' / 'DEXSeqResults.B_vs_A.csv'
    native.parent.mkdir()
    content = 'groupID,featureID,pvalue,padj,log2fold_A_B\nG1,E001,1e-310,2e-310,-0.7\nG1,E002,NA,NA,0.5\n'
    native.write_text(content)
    native_sha = hashlib.sha256(native.read_bytes()).hexdigest()
    manifest = dict(schema='katana-pipeline-bundle/v1', status='complete', reference={'genome_build':'generated-build'},
        samples=[dict(sample_id='A1',condition='A',subject='donor-a'), dict(sample_id='B1',condition='B',subject='donor-b')],
        comparisons=[dict(comparison_id='B_vs_A',treatment='B',control='A',rmats_group1='A',rmats_group2='B',deg_target='B',deg_reference='A')],
        artifacts=[dict(kind='auxiliary',relative_path='dexseq_exon/'+native.name,comparison_id='B_vs_A',status='available',bytes=native.stat().st_size,sha256=native_sha)])
    (bundle/'manifest.json').write_text(json.dumps(manifest))
    template=root/'template.json';template.write_text(json.dumps({'metadata':{'species':'generated'}}))
    result=subprocess.run(['katana-flow','import','--owned-link-root',str(root),'--bundle',str(bundle),'--project-template',str(template),'--outdir',str(root/'workspace')],text=True,capture_output=True)
    (root/'import.stdout.txt').write_text(result.stdout);(root/'import.stderr.txt').write_text(result.stderr)
    assert result.returncode==0,result.stderr
    receipt=json.loads((root/'workspace/import.receipt.json').read_text());assert receipt['status']=='completed'
    projection=receipt['format_projections'][0];path=root/'workspace/bundle'/projection['path']
    assert hashlib.sha256(path.read_bytes()).hexdigest()==projection['sha256']
    assert list(csv.reader(io.StringIO(content)))==list(csv.reader(io.StringIO(path.read_text()),delimiter='\t'))
    assert hashlib.sha256(native.read_bytes()).hexdigest()==native_sha
    checks.append('public installed import preserves native CSV statistics and TSV checksum')
    record=dict(status='passed',checks=checks,science_recomputed=False,file_io_sha256=hashlib.sha256(Path(file_io.__file__).read_bytes()).hexdigest())
    (root/'acceptance.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record))


if __name__=='__main__':
    main()
