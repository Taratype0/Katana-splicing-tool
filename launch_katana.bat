@echo off
setlocal

set ROOT=%~dp0

if exist "%ROOT%dist\KatanaSplicingTool\KatanaSplicingTool.exe" (
  start "" "%ROOT%dist\KatanaSplicingTool\KatanaSplicingTool.exe"
  exit /b 0
)

if exist "%ROOT%.venv\Scripts\python.exe" (
  "%ROOT%.venv\Scripts\python.exe" -m app.main
  exit /b %errorlevel%
)

python -m app.main
