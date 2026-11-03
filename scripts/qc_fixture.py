"""Generated native QC inputs for installed pipeline/report handoff acceptance."""
from pathlib import Path
import hashlib
import zipfile


def add_native_qc(project):
    from katana_next.infrastructure.manifest_store import ManifestStore
    from katana_next.domain.models import InputAsset
    path, manifest = ManifestStore().load(project)
    root = Path(project).parent / 'generated-native-qc'
    root.mkdir(exist_ok=False)
    assets = []
    def asset(file, kind, phase, sample=None):
        metadata = {'qc_type':kind, 'qc_phase':phase, 'upstream_sha256':hashlib.sha256(file.read_bytes()).hexdigest()}
        if sample: metadata['sample_id'] = sample
        assets.append(InputAsset(asset_id='generated-qc-'+str(len(assets)), kind='sequencing_qc', path=file, metadata=metadata))
    strand = ['sample\trequested\tstrandedness\tmethod\tinformative_forward_fraction\tinformative_reverse_fraction\tfailed_to_determine_fraction\tusable_reads_sampled\n']
    for i, sample in enumerate(manifest.samples):
        sid = sample.sample_id
        for mate in (1,2):
            for phase, name in [('raw',f'{sid}_{mate}_fastqc.zip'), ('trimmed',f'{sid}_{mate}_val_{mate}_fastqc.zip')]:
                file = root/name
                content = ('##FastQC\t0.12.1\n>>Basic Statistics\tpass\n#Measure\tValue\n'
                    f'Filename\t{sid}_{mate}.fastq.gz\nTotal Sequences\t1000\nSequence length\t150\n%GC\t50\n>>END_MODULE\n'
                    '>>Per base sequence quality\tpass\n#Base\tMean\tMedian\n1\t39\t40\n10-14\t38\t40\n>>END_MODULE\n'
                    '>>Sequence Length Distribution\tpass\n#Length\tCount\n150\t1000\n>>END_MODULE\n')
                with zipfile.ZipFile(file,'w') as z:z.writestr(sid+'/fastqc_data.txt',content)
                asset(file,'fastqc',phase,sid)
                if phase=='raw':
                    html=file.with_suffix('.html');html.write_text('<html><body>Generated native FastQC companion</body></html>');asset(html,'fastqc',phase,sid)
            file=root/f'{sid}_{mate}.fastq.gz_trimming_report.txt'
            file.write_text('Total reads processed: 1,000\nReads with adapters: 250 (25.0%)\nReads written (passing filters): 900 (90.0%)\n')
            asset(file,'trim_galore','trimming',sid)
        file=root/f'{sid}Log.final.out'
        file.write_text('Number of input reads | 1000\nUniquely mapped reads % | 80%\n% of reads mapped to multiple loci | 10%\n% of reads mapped to too many loci | 1%\n% of reads unmapped: too many mismatches | 1%\n% of reads unmapped: too short | 7%\n% of reads unmapped: other | 1%\nAverage input read length | 297\n')
        asset(file,'star','alignment',sid)
        if i%2:strand.append(f'{sid}\treverse\treverse\tuser_declared_not_inferred\t\t\t\t\n')
        else:strand.append(f'{sid}\tauto\treverse\trseqc_infer_experiment\t0.04\t0.96\t0.02\t200000\n')
    file=root/'strandedness.summary.tsv';file.write_text(''.join(strand));asset(file,'strandedness','library')
    manifest=manifest.model_copy(update={'assets':tuple(a for a in manifest.assets if a.kind!='sequencing_qc')+tuple(assets)})
    ManifestStore().save(manifest,project)
    return project
