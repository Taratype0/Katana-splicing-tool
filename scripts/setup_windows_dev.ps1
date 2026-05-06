param(
    [string]$PythonExe = "python"
)

$ErrorActionPreference = "Stop"
$root = Split-Path -Parent $PSScriptRoot
Set-Location $root

if (-not (Test-Path ".venv")) {
    & $PythonExe -m venv .venv
}

& ".\.venv\Scripts\python.exe" -m pip install --upgrade pip | Out-Host
& ".\.venv\Scripts\python.exe" -m pip install -e . | Out-Host

Write-Host ""
Write-Host "Development environment is ready."
Write-Host "Run with: .\\launch_katana.bat"
