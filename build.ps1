# build.ps1 - Windows build (PowerShell, incl. VS Code terminal). ASCII-only on purpose:
# Windows PowerShell 5.1 reads BOM-less .ps1 as ANSI, so non-ASCII text breaks parsing.
#
# Run from the repository root:
#   powershell -ExecutionPolicy Bypass -File .\build.ps1
#   powershell -ExecutionPolicy Bypass -File .\build.ps1 tests\test_xxx.cpp
#   (or simply  .\build.ps1  if local scripts are allowed)
#
# Needs MinGW in PATH (g++, gcc, ar). Check:  g++ --version
# Builds libsofa.a (if missing) + the given main (default tests\verify.cpp) and runs it.
# The programs themselves set the console to UTF-8, so Russian output shows correctly.

try { [Console]::OutputEncoding = [System.Text.Encoding]::UTF8 } catch {}
$null = & chcp 65001

$ErrorActionPreference = "Stop"
Set-Location $PSScriptRoot

$Main = if ($args.Count -ge 1) { $args[0] } else { "delay_tool.cpp" }
$Out  = [System.IO.Path]::ChangeExtension($Main, ".exe")

$SofaDir = "external\sofa\20190722\c"
$SofaLib = "$SofaDir\build\libsofa.a"

# --- 1) libsofa.a (once) ---
if (-not (Test-Path $SofaLib)) {
    Write-Host "[1/2] Building libsofa.a ..."
    New-Item -ItemType Directory -Force -Path "$SofaDir\build" | Out-Null
    Get-ChildItem "$SofaDir\src\*.c" | Where-Object { $_.Name -ne "t_sofa_c.c" } | ForEach-Object {
        gcc -O2 -c -I"$SofaDir\src" $_.FullName -o "$SofaDir\build\$($_.BaseName).o"
    }
    if (Test-Path $SofaLib) { Remove-Item $SofaLib }
    $objs = (Get-ChildItem "$SofaDir\build\*.o").FullName
    ar rcs $SofaLib $objs
    Write-Host "  libsofa.a ready."
} else {
    Write-Host "[1/2] libsofa.a exists, skipping."
}

# --- 2) main build ---
Write-Host "[2/2] Building $Out ..."
$inc = @("-I.\external\eigen", "-I.\external", "-I.\external\dephem-master\include", "-I.\external\sofa\20190722\c\src")
$src = @(
  "src\process_ariadna.cpp","src\process_obs.cpp","src\site_pair.cpp","src\site_tide_solid.cpp",
  "src\site_tide_OC.cpp","src\pole_tide.cpp","src\site_atm40.cpp","src\therm_def40.cpp","src\site_inst.cpp",
  "src\baseline.cpp","src\aber_source.cpp","src\trop_delay.cpp","src\nhmf2.cpp","src\nwmf2.cpp",
  "src\sast_dry.cpp","src\sast_wet.cpp","src\mount_tel.cpp","src\sbend.cpp","src\theor_delay.cpp",
  "src\jpleph.cpp","src\ephemeris.cpp","src\fund_arg.cpp","src\GEOD.cpp","src\site_functions.cpp",
  "src\rotation.cpp","src\dmeteo.cpp","src\orbit_interp.cpp","src\t_eph40.cpp","src\tai_time40.cpp",
  "src\nsec.cpp","src\interp_eop.cpp","src\terms_71.cpp","src\terms_lib.cpp","src\UT1R_2010.cpp",
  "src\READ_CAT.cpp","src\catalog_bridge.cpp"
)
$dep = ".\external\dephem-master\include\dephem\EphemerisRelease.cpp"

$gccArgs = @("-std=c++17") + $inc + $src + @($dep, $Main, $SofaLib, "-o", $Out)
& g++ @gccArgs
if ($LASTEXITCODE -ne 0) { Write-Host "BUILD FAILED."; exit 1 }
Write-Host "  OK -> $Out"

Write-Host "Running:"
& ".\$Out"
