@echo off
rem ============================================================================
rem test.bat - Windows end-to-end module test (cmd/PowerShell). ASCII-only + CRLF.
rem
rem Builds and runs test_process_task: process_task() computes delay polynomials for
rem all stations of the example session (.cfx + .scf) and compares P0 of every block
rem with the reference example\*.TXT. Prints per-station error and PASS/FAIL.
rem
rem Run from the repository root:
rem   .\test.bat
rem
rem Prefer test.ps1 in PowerShell/VS Code. Needs MinGW in PATH and the ephemeris file
rem in external\dephem-master\ (see BUILD.md). Exit code 0 = PASS, non-zero = FAIL.
rem ============================================================================
setlocal
cd /d "%~dp0"

call build.bat tests\test_process_task.cpp
set "CODE=%ERRORLEVEL%"
echo.
if "%CODE%"=="0" (echo TEST PASSED [exit 0].) else (echo TEST FAILED [exit %CODE%].)
endlocal & exit /b %CODE%
