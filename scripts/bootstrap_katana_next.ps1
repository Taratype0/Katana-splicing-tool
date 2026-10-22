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
if ($LASTEXITCODE -ne 0) { throw "Katana Python dependency installation failed." }

if ($SkipFrontendBuild) {
    if (-not (Test-Path -LiteralPath (Join-Path $repoRoot "apps\web\dist\index.html"))) {
        throw "The prebuilt interface is missing. Run this script without -SkipFrontendBuild."
    }
} else {
    $pnpm = Get-Command pnpm -ErrorAction SilentlyContinue
    if ($null -eq $pnpm) {
        throw "pnpm is required to build the Katana Next interface."
    }
    $previousCI = $env:CI
    Push-Location (Join-Path $repoRoot "apps\web")
    try {
        $env:CI = "true"
        & $pnpm.Source install --frozen-lockfile
        if ($LASTEXITCODE -ne 0) { throw "Katana interface dependency installation failed." }
        & $pnpm.Source build
        if ($LASTEXITCODE -ne 0) { throw "Katana interface build failed." }
    } finally {
        $env:CI = $previousCI
        Pop-Location
    }
}

Write-Host "Katana Next development environment is ready at $repoRoot"
