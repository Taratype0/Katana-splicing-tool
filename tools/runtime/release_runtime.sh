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
  test "$(gh run view "$SOFTWARE_RUN" --json conclusion --jq .conclusion)" = success
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
  mkdir controller
  tar xzf "$(find accepted-controller -name 'katana-controller-*.tar.gz')" -C controller
  echo "$ROOT/controller/katana-controller/bin" >> "$GITHUB_PATH"
  echo "CONTROLLER=$ROOT/controller/katana-controller" >> "$GITHUB_ENV"
  gh release download v2.0.0-beta.5.dev20261102 -R "$GITHUB_REPOSITORY" -p katana-downstream-2.0.0b5.dev20261102-linux-amd64.docker.tar.gz
  echo '43161f00122043848cd7b3432d504551c19b2d2f2f3f89d0b6a484d64c5a989d  katana-downstream-2.0.0b5.dev20261102-linux-amd64.docker.tar.gz' | sha256sum --check --strict
  echo 'cc6b9d9c17b98fdcbf59b0dfffe77c20356b786830c7dc503ac2dc5f4b434b99  environments/locks/linux-64.explicit.txt' | (cd "$GITHUB_WORKSPACE"; sha256sum --check --strict)
  echo '4bbfdd6b3510821fce032d52d42a52bef1313238fad82b7c2d831adda5cb87f1  environments/locks/pip-linux-64.lock' | (cd "$GITHUB_WORKSPACE"; sha256sum --check --strict)
  docker load -i katana-downstream-2.0.0b5.dev20261102-linux-amd64.docker.tar.gz
  docker tag sha256:79c0746d637155fb35aa4137f11b91f5050f85bf715847ae3d680de9a69001d3 katana-runtime-base:verified
  ;;
build)
  cd "$ROOT"
  cat > context/Dockerfile <<'DOCKER'
FROM katana-runtime-base:verified
COPY *.whl /tmp/katana-wheel/
RUN python -m pip install --no-cache-dir --no-deps --force-reinstall /tmp/katana-wheel/*.whl && python -m pip check && rm -rf /tmp/katana-wheel
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
apptainer)
  cd "$ROOT"
  curl --fail --location --retry 3 -o apptainer.deb https://github.com/apptainer/apptainer/releases/download/v1.4.3/apptainer_1.4.3_amd64.deb
  echo 'c3c05120baab7ece2e81d40c1d1be63cc3fd1f5a25f379a4f2b2f99a0f6a1100  apptainer.deb' | sha256sum --check --strict
  sudo apt-get update -qq
  sudo apt-get install -y ./apptainer.deb
  docker save "$IMAGE" -o runtime.tar
  sudo apptainer build "output/katana-downstream-$RUNTIME_VERSION-linux-amd64.sif" "docker-archive://$ROOT/runtime.tar"
  sudo chown "$(id -u):$(id -g)" output/*.sif
  apptainer exec --cleanenv --bind "$ROOT/output:/acceptance,$GITHUB_WORKSPACE/scripts:/checks:ro" "output/katana-downstream-$RUNTIME_VERSION-linux-amd64.sif" python /checks/check_downstream_environment.py --out /acceptance/apptainer-environment
  gzip -1 runtime.tar
  mv runtime.tar.gz "output/katana-downstream-$RUNTIME_VERSION-linux-amd64.docker.tar.gz"
  ;;
collect)
  cd "$ROOT/output"
  sha256sum *.docker.tar.gz *.sif > runtime-SHA256SUMS
  python - <<'PY'
from pathlib import Path
import json,os,tarfile
summaries={name:json.loads(Path(path).read_text()) for name,path in {
  'native_runtime':'native-runtime/acceptance.receipt.json',
  'standalone_models':'standalone-models/acceptance.json',
  'native_models':'native-models/acceptance.json'}.items()}
assert all(s['status']=='passed' for s in summaries.values())
Path('runtime.acceptance.json').write_text(json.dumps(dict(status='passed',candidate=os.environ['GITHUB_SHA'],wheel_sha256=os.environ['WHEEL_SHA'],software_run=os.environ['SOFTWARE_RUN'],checks=summaries),indent=2))
with tarfile.open('runtime-validation.tar.gz','w:gz') as archive:
  for path in Path('.').rglob('*'):
    if path.is_file() and path.suffix in ('.json','.tsv','.log','.txt') and path.stat().st_size<5000000:
      archive.add(path,arcname=str(path))
PY
  test "$(gh release view "$RELEASE_TAG" -R "$GITHUB_REPOSITORY" --json isDraft --jq .isDraft)" = true
  gh release upload "$RELEASE_TAG" -R "$GITHUB_REPOSITORY" *.docker.tar.gz *.sif runtime-SHA256SUMS runtime-validation.tar.gz runtime.acceptance.json
  ;;
*) echo 'Unknown release acceptance phase' >&2; exit 2;;
esac
