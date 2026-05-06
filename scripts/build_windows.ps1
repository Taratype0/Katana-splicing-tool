param(
    [string]$PythonExe = "python"
)

$ErrorActionPreference = "Stop"
$root = Split-Path -Parent $PSScriptRoot
Set-Location $root

Write-Host "Building KatanaSplicingTool for Windows..."

& $PythonExe -m pip install pyinstaller | Out-Host
& $PythonExe -m PyInstaller build.spec --noconfirm --clean | Out-Host

Write-Host ""
Write-Host "Build finished."
Write-Host "Output: $root\\dist\\KatanaSplicingTool"
