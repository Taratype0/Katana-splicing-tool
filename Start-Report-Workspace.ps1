param([switch]$NoOpen)
$ErrorActionPreference = 'Stop'
$katanaPython = 'E:/Project/katana-splice-enhanced-20261017/app/.venv-next/Scripts/python.exe'
$katanaUrl = 'http://127.0.0.1:18774'
$katanaReady = $false
try {
    $katanaIdentity = Invoke-RestMethod "$katanaUrl/api/report-workspace/version" -TimeoutSec 3
    $katanaReady = $katanaIdentity.release -eq 'report-workspace-20261108'
} catch {}
if (-not $katanaReady) {
    $katanaLogs = Join-Path $env:LOCALAPPDATA 'Katana/report-workspace-20261108'
    New-Item -ItemType Directory -Force -Path $katanaLogs | Out-Null
    Start-Process -FilePath $katanaPython -ArgumentList @('-B','-X','utf8','-m','katana_next.api.main','--port','18774') -WorkingDirectory $PSScriptRoot -WindowStyle Hidden -RedirectStandardOutput (Join-Path $katanaLogs 'service.stdout.log') -RedirectStandardError (Join-Path $katanaLogs 'service.stderr.log') | Out-Null
    for ($katanaAttempt = 0; $katanaAttempt -lt 30; $katanaAttempt++) {
        Start-Sleep -Milliseconds 500
        try {
            $katanaIdentity = Invoke-RestMethod "$katanaUrl/api/report-workspace/version" -TimeoutSec 2
            if ($katanaIdentity.release -eq 'report-workspace-20261108') { $katanaReady=$true; break }
        } catch {}
    }
}
if (-not $katanaReady) { throw 'Report workspace failed to start. Check port 18774 and the service logs.' }
if (-not $NoOpen) { Start-Process $katanaUrl }
Write-Output $katanaUrl
