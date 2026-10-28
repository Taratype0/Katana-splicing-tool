"""Install a new Katana wheel into a verified release runtime on Linux HPC.

Uses an owned, writable Singularity sandbox and produces a new immutable SIF.
The release archive and existing images remain unchanged. Dependency locks must
match the release: this is an application update, not an environment upgrade.
"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--archive', type=Path, required=True)
    ap.add_argument('--release-manifest', type=Path, required=True)
    ap.add_argument('--wheel', type=Path, required=True)
    ap.add_argument('--outdir', type=Path, required=True)
    ap.add_argument('--engine', default='singularity')
    ap.add_argument('--scratch-dir', type=Path, default=Path('/tmp'),
                    help='Node-local build filesystem; only the final SIF is retained in outdir')
    args = ap.parse_args()
    archive, wheel = args.archive.resolve(strict=True), args.wheel.resolve(strict=True)
    manifest = json.loads(args.release_manifest.read_text())
    item = next(v for v in manifest['containers'] if v['filename'] == archive.name)
    if archive.stat().st_size != item['bytes'] or sha(archive) != item['sha256']:
        raise ValueError('Release archive checksum/size differs from manifest')
    for field, name in [('native_dependency_lock_sha256', 'linux-64.explicit.txt'),
                        ('pip_dependency_lock_sha256', 'pip-linux-64.lock')]:
        if sha(ROOT / 'environments/locks' / name) != manifest[field]:
            raise ValueError('Dependency locks changed; use a full environment build instead')
    out = args.outdir.resolve()
    if out == ROOT or ROOT in out.parents:
        raise ValueError('Build outside the checkout')
    out.mkdir(parents=True, exist_ok=False)
    receipt = dict(schema='katana.application_sif_update/v1', status='building',
                   base_archive_sha256=item['sha256'], base_image_id=item['image_id'],
                   wheel_sha256=sha(wheel), wheel=wheel.name, dependency_locks_changed=False,
                   science_data_used=False)
    receipt_path = out / 'build.receipt.json'
    receipt_path.write_text(json.dumps(receipt, indent=2) + '\n')
    def run(argv):
        with (out / 'build.log').open('a') as log:
            log.write('\n' + repr(argv) + '\n'); log.flush()
            subprocess.run(argv, stdout=log, stderr=subprocess.STDOUT, check=True)
    sif = out / 'katana.sif'
    # Large unpacked trees stay on the node, not in a quota-limited project.
    with tempfile.TemporaryDirectory(prefix='katana-sif-build-', dir=args.scratch_dir.resolve(strict=True)) as directory:
        work = Path(directory)
        tar = work / 'base.docker.tar'
        with gzip.open(archive, 'rb') as source, tar.open('xb') as target:
            shutil.copyfileobj(source, target, 8 * 1024 * 1024)
        sandbox = work / 'sandbox'
        run([args.engine, 'build', '--disable-cache', '--sandbox', str(sandbox), 'docker-archive://' + str(tar)])
        (sandbox / 'katana-wheel').mkdir()
        run([args.engine, 'exec', '--writable', '--no-home', '--cleanenv', '--pwd', '/tmp',
             '--bind', str(wheel.parent) + ':/katana-wheel', str(sandbox),
             '/opt/conda/bin/python', '-m', 'pip', 'install', '--no-cache-dir', '--no-deps',
             '--force-reinstall', '/katana-wheel/' + wheel.name])
        run([args.engine, 'exec', '--writable', '--no-home', '--cleanenv', '--pwd', '/tmp', str(sandbox),
             '/opt/conda/bin/python', '-m', 'pip', 'check'])
        local_sif = work / 'katana.sif'
        run([args.engine, 'build', '--disable-cache', str(local_sif), str(sandbox)])
        shutil.copyfile(local_sif, sif)
        if sha(local_sif) != sha(sif):
            raise RuntimeError('Copied SIF checksum differs')
    receipt.update(status='built', sif=str(sif), sif_sha256=sha(sif),
                   installation_acceptance='pending', temporary_build_paths_removed=True)
    receipt_path.write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
