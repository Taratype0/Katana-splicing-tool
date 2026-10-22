@echo off
setlocal
cd /d "%~dp0"
if not exist ".venv-next\Scripts\python.exe" (
  echo Run powershell -File scripts\bootstrap_katana_next.ps1 -SkipFrontendBuild first.
  exit /b 1
)
".venv-next\Scripts\python.exe" -m katana_next.launcher --mode browser %*
