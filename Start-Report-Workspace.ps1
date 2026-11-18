param([switch]$NoOpen)
& (Join-Path $PSScriptRoot "Start-Local-Katana.ps1") -NoOpen:$NoOpen
