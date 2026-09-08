@echo off
setlocal
cd /d "%~dp0"

if not exist "%~dp0index.html" (
    echo Building standalone LL Analyzer...
    call node build-standalone.mjs
)

start "" "%~dp0index.html"
