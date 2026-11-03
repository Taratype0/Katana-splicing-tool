#!/usr/bin/env bash
# CI-only release checks. Uses generated fixtures and never study data.
set -euo pipefail
ROOT="$RUNNER_TEMP/katana-release"
IMAGE="katana-downstream:$RUNTIME_VERSION"
case "$1" in
prepare)
  [[ "$RUNTIME_VERSION" =~ ^[0-9][0-9a-zA-Z.+-]*$ ]]
  [[ "$SOFTWARE_RUN" =~ ^[0-9]+$ ]]
  test "$(gh release view "$RELEASE_TAG" --json isDraft --jq .isDraft)" = true
  gh run view "$SOFTWARE_RUN" --json headSha,jobs > /tmp/katana-software-checks.json
  python - <<'PYCHECK'
import json,os
from pathlib import Path
run=json.loads(Path('/tmp/katana-software-checks.json').read_text())
assert run['headSha']==os.environ['GITHUB_SHA'], 'Software checks must match release source'
required={'package-wheel','installed-native-import (linux)','installed-native-import (windows)'}
passed={j['name'] for j in run['jobs'] if j['status']=='completed' and j['conclusion']=='success'}
assert required <= passed, (required-passed)
PYCHECK
  mkdir -p "$ROOT/context" "$ROOT/output"
  cd "$ROOT"
  gh release download "$RELEASE_TAG" -R "$GITHUB_REPOSITORY" -p '*.whl' -D context
  test "$(find context -name '*.whl' | wc -l)" = 1
  printf '%s  %s\n' "$WHEEL_SHA" "$(find context -name '*.whl')" | sha256sum --check --strict
  gh run download "$SOFTWARE_RUN" -R "$GITHUB_REPOSITORY" -n installed-native-linux -D accepted-controller
  python - <<'PY'
import json,os,hashlib
from pathlib import Path
receipt=json.loads(next(Path('accepted-controller').rglob('build.receipt.json')).read_text())
assert receipt['wheel_sha256']==os.environ['WHEEL_SHA'], receipt
archive=next(Path('accepted-controller').rglob('katana-controller-*.tar.gz'))
assert hashlib.file_digest(archive.open('rb'),'sha256').hexdigest()==receipt['archive_sha256']
PY
  python - <<'PYCONTROLLER'
import tarfile,hashlib,json,os
from pathlib import Path
archive=next(Path('accepted-controller').rglob('katana-controller-*.tar.gz'))
checked=[]
with tarfile.open(archive) as tar:
 for source in ('katana_next/flow/file_io.py','katana_next/flow/importer.py','katana_next/flow/_pipeline/nextflow.config'):
  matches=[m for m in tar.getmembers() if m.name.endswith('/'+source)];assert len(matches)==1,source
  local=source.replace('katana_next/flow/_pipeline/','pipeline/')
  content=tar.extractfile(matches[0]).read();assert content==(Path(os.environ['GITHUB_WORKSPACE'])/local).read_bytes(),source
  checked.append(dict(path=source,sha256=hashlib.sha256(content).hexdigest()))
record=dict(status='passed',source_commit=os.environ['GITHUB_SHA'],software_run=os.environ['SOFTWARE_RUN'],controller=dict(filename=archive.name,bytes=archive.stat().st_size,sha256=hashlib.file_digest(archive.open('rb'),'sha256').hexdigest()),packaged_sources=checked,rebuild=False)
Path('output/CONTROLLER_VALIDATION.json').write_text(json.dumps(record,indent=2)+'\n')
PYCONTROLLER
  mkdir controller
  tar xzf "$(find accepted-controller -name 'katana-controller-*.tar.gz')" -C controller
  echo "$ROOT/controller/katana-controller/bin" >> "$GITHUB_PATH"
  echo "CONTROLLER=$ROOT/controller/katana-controller" >> "$GITHUB_ENV"
  gh release download v2.0.0-beta.5.dev20261102 -R "$GITHUB_REPOSITORY" -p katana-downstream-2.0.0b5.dev20261102-linux-amd64.docker.tar.gz
  echo '43161f00122043848cd7b3432d504551c19b2d2f2f3f89d0b6a484d64c5a989d  katana-downstream-2.0.0b5.dev20261102-linux-amd64.docker.tar.gz' | sha256sum --check --strict
  echo 'cc6b9d9c17b98fdcbf59b0dfffe77c20356b786830c7dc503ac2dc5f4b434b99  environments/locks/linux-64.explicit.txt' | (cd "$GITHUB_WORKSPACE"; sha256sum --check --strict)
  echo '4bbfdd6b3510821fce032d52d42a52bef1313238fad82b7c2d831adda5cb87f1  environments/locks/pip-linux-64.lock' | (cd "$GITHUB_WORKSPACE"; sha256sum --check --strict)
  docker load -i katana-downstream-2.0.0b5.dev20261102-linux-amd64.docker.tar.gz
  # The verified archive supplies this unique tag. Docker storage backends may
  # expose a manifest digest or a config digest as their local image identifier.
  docker inspect katana-downstream:2.0.0b5.dev20261102 > output/base-image-inspect.json
  docker tag katana-downstream:2.0.0b5.dev20261102 katana-runtime-base:verified
  ;;
build)
  cd "$ROOT"
  cat > context/Dockerfile <<'DOCKER'
FROM katana-runtime-base:verified
COPY *.whl /tmp/katana-candidate-install/
RUN python -m pip install --no-cache-dir --no-deps --force-reinstall /tmp/katana-candidate-install/*.whl && python -m pip check && rm -rf /tmp/katana-candidate-install /tmp/katana-wheel
DOCKER
  docker build --network=none -t "$IMAGE" context
  test "$(docker run --rm --network=none "$IMAGE" python -c 'from importlib.metadata import version; print(version("katana-splicing-tool"))')" = "$RUNTIME_VERSION"
  docker run --rm --network=none --user "$(id -u):$(id -g)" -e NUMBA_CACHE_DIR=/tmp/numba -v "$ROOT/output:/acceptance" -v "$GITHUB_WORKSPACE/scripts:/checks:ro" "$IMAGE" python /checks/check_downstream_environment.py --out /acceptance/environment
  docker inspect "$IMAGE" > output/image-inspect.json
  cp "$CONTROLLER/CONTROLLER_MANIFEST.json" output/controller-manifest.json
  cp "$(find accepted-controller -name build.receipt.json)" output/controller-build.receipt.json
  ;;
models)
  S="$ROOT/output"
  docker run --rm --network=none --user "$(id -u):$(id -g)" -e NUMBA_CACHE_DIR=/tmp/numba -e MPLCONFIGDIR=/tmp/matplotlib -e OMP_NUM_THREADS=2 -e OPENBLAS_NUM_THREADS=2 -e MKL_NUM_THREADS=2 -v "$S:$S" -v "$GITHUB_WORKSPACE/scripts:/checks:ro" -w "$S" "$IMAGE" python /checks/check_unified_engines.py --out "$S/standalone-models" 2>&1 | tee "$S/standalone-models.log"
  ;;
collect-docker)
  cd "$ROOT/output"
  docker save "$IMAGE" | gzip -1 > "katana-downstream-$RUNTIME_VERSION-linux-amd64.docker.tar.gz"
  sha256sum *.docker.tar.gz > docker-SHA256SUMS
  python - <<'PYDOCKER'
from pathlib import Path
import json,os,tarfile
summaries={name:json.loads(Path(path).read_text()) for name,path in {
  'native_runtime':'native-runtime/acceptance.receipt.json',
  'standalone_models':'standalone-models/acceptance.json',
  'native_models':'native-models/acceptance.json',
  'environment':'environment/acceptance.json',
  'import_content':'import-content/acceptance.json'}.items()}
assert all(s['status']=='passed' for s in summaries.values())
Path('docker.acceptance.json').write_text(json.dumps(dict(status='passed',candidate=os.environ['GITHUB_SHA'],wheel_sha256=os.environ['WHEEL_SHA'],software_run=os.environ['SOFTWARE_RUN'],checks=summaries),indent=2))
with tarfile.open('docker-validation.tar.gz','w:gz') as archive:
  for path in Path('.').rglob('*'):
    if path.is_file() and path.suffix in ('.json','.tsv','.log','.txt') and path.stat().st_size<5000000:
      archive.add(path,arcname=str(path))
PYDOCKER
  test "$(gh release view "$RELEASE_TAG" -R "$GITHUB_REPOSITORY" --json isDraft --jq .isDraft)" = true
  gh release upload "$RELEASE_TAG" -R "$GITHUB_REPOSITORY" *.docker.tar.gz docker-SHA256SUMS docker-validation.tar.gz docker.acceptance.json
  gh release upload "$RELEASE_TAG" -R "$GITHUB_REPOSITORY" "$(find "$ROOT/accepted-controller" -name 'katana-controller-*.tar.gz')" CONTROLLER_VALIDATION.json
  ;;
apptainer)
  mkdir -p "$ROOT/output"
  cd "$ROOT/output"
  gh release download "$RELEASE_TAG" -R "$GITHUB_REPOSITORY" -p "katana-downstream-$RUNTIME_VERSION-linux-amd64.docker.tar.gz" -p docker-SHA256SUMS -p docker.acceptance.json
  sha256sum --check --strict docker-SHA256SUMS
  python - <<'PYACCEPT'
import json,os
from pathlib import Path
r=json.loads(Path('docker.acceptance.json').read_text())
assert r['status']=='passed' and r['wheel_sha256']==os.environ['WHEEL_SHA']
assert r['candidate']==os.environ['GITHUB_SHA']
PYACCEPT
  cd "$ROOT"
  curl --fail --location --retry 3 -o apptainer.deb https://github.com/apptainer/apptainer/releases/download/v1.4.3/apptainer_1.4.3_amd64.deb
  echo 'c3c05120baab7ece2e81d40c1d1be63cc3fd1f5a25f379a4f2b2f99a0f6a1100  apptainer.deb' | sha256sum --check --strict
  sudo apt-get update -qq
  sudo apt-get install -y ./apptainer.deb
  gzip -dc "output/katana-downstream-$RUNTIME_VERSION-linux-amd64.docker.tar.gz" > runtime.tar
  sudo apptainer build "output/katana-downstream-$RUNTIME_VERSION-linux-amd64.sif" "docker-archive://$ROOT/runtime.tar"
  sudo chown "$(id -u):$(id -g)" output/*.sif
  apptainer exec --cleanenv --bind "$ROOT/output:/acceptance,$GITHUB_WORKSPACE/scripts:/checks:ro" "output/katana-downstream-$RUNTIME_VERSION-linux-amd64.sif" python /checks/check_downstream_environment.py --out /acceptance/apptainer-environment
  apptainer exec --cleanenv --bind "$ROOT/output:/acceptance,$GITHUB_WORKSPACE/tools/runtime:/checks:ro" "output/katana-downstream-$RUNTIME_VERSION-linux-amd64.sif" python /checks/check_import_content.py --outdir /acceptance/apptainer-import-content
  ;;
collect)
  cd "$ROOT/output"
  sha256sum *.docker.tar.gz *.sif > runtime-SHA256SUMS
  python - <<'PYACCEPT'
from pathlib import Path
import json,os,tarfile
receipt=json.loads(Path('docker.acceptance.json').read_text())
sif=json.loads(Path('apptainer-environment/acceptance.json').read_text())
assert receipt['status']=='passed' and sif['status']=='passed'
receipt['apptainer']=sif
import_check=json.loads(Path('apptainer-import-content/acceptance.json').read_text()); assert import_check['status']=='passed'
receipt['apptainer_import_content']=import_check
Path('runtime.acceptance.json').write_text(json.dumps(receipt,indent=2))
with tarfile.open('apptainer-validation.tar.gz','w:gz') as archive:
  archive.add('apptainer-environment')
  archive.add('apptainer-import-content')
PYACCEPT
  test "$(gh release view "$RELEASE_TAG" -R "$GITHUB_REPOSITORY" --json isDraft --jq .isDraft)" = true
  gh release upload "$RELEASE_TAG" -R "$GITHUB_REPOSITORY" *.sif runtime-SHA256SUMS apptainer-validation.tar.gz runtime.acceptance.json
  ;;
*) echo 'Unknown release acceptance phase' >&2; exit 2;;
esac
