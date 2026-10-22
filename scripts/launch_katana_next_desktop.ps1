param(
    [int]$Port = 8765,
    [switch]$RebuildFrontend
)

$ErrorActionPreference = "Stop"
$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$python = Join-Path $repoRoot ".venv-next\Scripts\python.exe"
$frontend = Join-Path $repoRoot "apps\web\dist\index.html"

if (-not (Test-Path -LiteralPath $python)) {
    throw "Katana Next environment is missing. Run scripts\bootstrap_katana_next.ps1 first."
}

if ($RebuildFrontend -or -not (Test-Path -LiteralPath $frontend)) {
    & (Join-Path $PSScriptRoot "bootstrap_katana_next.ps1")
}

Push-Location $repoRoot
try {
    & $python -m katana_next.launcher --mode desktop --port $Port
    if ($LASTEXITCODE -ne 0) { throw "Katana desktop startup failed." }
} finally { Pop-Location }
