# test.ps1 - Windows end-to-end test (PowerShell, incl. VS Code terminal). ASCII-only.
#
# Builds and runs the integration test: process_task() computes delay polynomials for
# ALL stations of the example session (.cfx + .scf) and compares P0 of every block with
# the reference example/*.TXT. Prints per-station error and PASS/FAIL.
#
# Run from the repository root:
#   powershell -ExecutionPolicy Bypass -File .\test.ps1
#
# Needs MinGW in PATH (g++, gcc, ar) and the ephemeris file in external\dephem-master\
# (see BUILD.md). Exit code 0 = PASS, non-zero = FAIL.

$ErrorActionPreference = "Stop"
Set-Location $PSScriptRoot

Write-Host "Building and running the end-to-end module test (test_process_task) ..."
& powershell -ExecutionPolicy Bypass -File .\build.ps1 "tests\test_process_task.cpp"
$code = $LASTEXITCODE
Write-Host ""
if ($code -eq 0) { Write-Host "TEST PASSED (exit 0)." } else { Write-Host "TEST FAILED (exit $code)." }
exit $code
