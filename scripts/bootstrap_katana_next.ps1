param(
    [switch]$SkipFrontendBuild
)

$ErrorActionPreference = "Stop"
$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$venvRoot = Join-Path $repoRoot ".venv-next"
$python = Join-Path $venvRoot "Scripts\python.exe"

if (-not (Test-Path -LiteralPath $python)) {
    python -m venv $venvRoot
}

& $python -m pip install --disable-pip-version-check -e "$repoRoot[next,next-dev]"

$pnpm = Get-Command pnpm -ErrorAction SilentlyContinue
if ($null -eq $pnpm) {
    $fallback = Join-Path $env:USERPROFILE ".cache\codex-runtimes\codex-primary-runtime\dependencies\bin\fallback\pnpm.cmd"
    if (Test-Path -LiteralPath $fallback) {
        $pnpm = Get-Item $fallback
    } else {
        throw "pnpm is required to build the Katana Next interface."
    }
}

Push-Location (Join-Path $repoRoot "apps\web")
try {
    $env:CI = "true"
    & $pnpm.Source install
    if (-not $SkipFrontendBuild) {
        & $pnpm.Source build
    }
} finally {
    Pop-Location
}

Write-Host "Katana Next development environment is ready at $repoRoot"
