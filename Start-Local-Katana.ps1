param([switch]$NoOpen, [int]$Port = 18775)
$ErrorActionPreference = 'Stop'
if (-not $env:KATANA_WORKSPACE_HOME) {
    $env:KATANA_WORKSPACE_HOME = Join-Path ([Environment]::GetFolderPath('LocalApplicationData')) 'Katana/workspaces'
}
$katanaPython = if ($env:KATANA_PYTHON) {
    $env:KATANA_PYTHON
} elseif (Test-Path (Join-Path $PSScriptRoot '.venv/Scripts/python.exe')) {
    Join-Path $PSScriptRoot '.venv/Scripts/python.exe'
} elseif (Test-Path (Join-Path $PSScriptRoot '.venv-next/Scripts/python.exe')) {
    Join-Path $PSScriptRoot '.venv-next/Scripts/python.exe'
} else {
    (Get-Command python -ErrorAction SilentlyContinue).Source
}
if (-not (Test-Path -LiteralPath $katanaPython -PathType Leaf)) { throw 'Katana Python is unavailable. Set KATANA_PYTHON to the prepared Python environment.' }
$katanaUrl = "http://127.0.0.1:$Port"
$katanaReady = $false
$katanaIdentity = $null
try { $katanaIdentity = Invoke-RestMethod "$katanaUrl/api/report-workspace/version" -TimeoutSec 3 } catch {}
if ($katanaIdentity) {
    $katanaReady = $katanaIdentity.release -eq 'katana-v2.7.0' -and [IO.Path]::GetFullPath($katanaIdentity.application_root) -eq [IO.Path]::GetFullPath($PSScriptRoot)
    if (-not $katanaReady) { throw "A different Katana service is using port $Port. Close that service or choose another -Port." }
}
if (-not $katanaReady) {
    if (Test-Path -LiteralPath (Join-Path $PSScriptRoot 'release-integrity.json')) {
        & $katanaPython -B (Join-Path $PSScriptRoot 'scripts/verify_local_release.py')
        if ($LASTEXITCODE -ne 0) { throw 'Katana release files failed verification.' }
    }
    $katanaLogs = Join-Path ([Environment]::GetFolderPath('LocalApplicationData')) 'Katana/katana-v2.7.0'
    New-Item -ItemType Directory -Force -Path $katanaLogs | Out-Null
    $katanaProcess = Start-Process -FilePath $katanaPython -ArgumentList @('-B','-X','utf8','-m','katana_next.api.main','--port',"$Port") -WorkingDirectory $PSScriptRoot -WindowStyle Hidden -PassThru -RedirectStandardOutput (Join-Path $katanaLogs "service-$Port.stdout.log") -RedirectStandardError (Join-Path $katanaLogs "service-$Port.stderr.log")
    for ($katanaAttempt = 0; $katanaAttempt -lt 60; $katanaAttempt++) {
        Start-Sleep -Milliseconds 500
        if ($katanaProcess.HasExited) { throw "Katana exited during startup. See $katanaLogs" }
        try {
            $katanaIdentity = Invoke-RestMethod "$katanaUrl/api/report-workspace/version" -TimeoutSec 2
            if ($katanaIdentity.release -eq 'katana-v2.7.0' -and [IO.Path]::GetFullPath($katanaIdentity.application_root) -eq [IO.Path]::GetFullPath($PSScriptRoot)) { $katanaReady=$true; break }
        } catch {}
    }
}
if (-not $katanaReady) { throw 'Katana failed to become ready. Check the local service logs.' }
if (-not $NoOpen) { Start-Process $katanaUrl }
Write-Output $katanaUrl
