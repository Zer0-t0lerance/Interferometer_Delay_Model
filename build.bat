@echo off
rem ============================================================================
rem build.bat - Windows build (cmd/PowerShell). MinGW g++/gcc/ar in PATH.
rem ASCII-only + CRLF on purpose (cmd mis-parses non-ASCII / LF-only batch files).
rem
rem Run from the repository root:
rem   .\build.bat                     build libsofa.a (if missing) + tests\verify.exe and run
rem   .\build.bat tests\test_xxx.cpp  build another main with the same module set
rem
rem Prefer build.ps1 in PowerShell/VS Code. Programs set the console to UTF-8 themselves,
rem so their Russian output shows correctly.
rem ============================================================================
setlocal
chcp 65001 >nul
cd /d "%~dp0"

set "MAIN=%~1"
if "%MAIN%"=="" set "MAIN=delay_tool.cpp"
for %%F in ("%MAIN%") do set "OUT=%%~dpnF.exe"

set "SOFADIR=external\sofa\20190722\c"
set "SOFALIB=%SOFADIR%\build\libsofa.a"

if not exist "%SOFALIB%" (
  echo [1/2] Building libsofa.a ...
  if not exist "%SOFADIR%\build" mkdir "%SOFADIR%\build"
  for %%f in ("%SOFADIR%\src\*.c") do (
    if /I not "%%~nxf"=="t_sofa_c.c" gcc -O2 -c -I"%SOFADIR%\src" "%%f" -o "%SOFADIR%\build\%%~nf.o"
  )
  if exist "%SOFALIB%" del "%SOFALIB%"
  for %%o in ("%SOFADIR%\build\*.o") do ar rs "%SOFALIB%" "%%o" >nul
  echo   libsofa.a ready.
) else (
  echo [1/2] libsofa.a exists, skipping.
)

echo [2/2] Building %OUT% ...
set "INC=-I.\external\eigen -I.\external -I.\external\dephem-master\include -I.\external\sofa\20190722\c\src"
set "SRC=src\process_ariadna.cpp src\process_obs.cpp src\site_pair.cpp src\site_tide_solid.cpp src\site_tide_OC.cpp src\pole_tide.cpp src\site_atm40.cpp src\therm_def40.cpp src\site_inst.cpp src\baseline.cpp src\aber_source.cpp src\trop_delay.cpp src\nhmf2.cpp src\nwmf2.cpp src\sast_dry.cpp src\sast_wet.cpp src\mount_tel.cpp src\sbend.cpp src\theor_delay.cpp src\jpleph.cpp src\ephemeris.cpp src\fund_arg.cpp src\GEOD.cpp src\site_functions.cpp src\rotation.cpp src\dmeteo.cpp src\orbit_interp.cpp src\t_eph40.cpp src\tai_time40.cpp src\nsec.cpp src\interp_eop.cpp src\terms_71.cpp src\terms_lib.cpp src\UT1R_2010.cpp src\READ_CAT.cpp src\catalog_bridge.cpp"
set "DEP=.\external\dephem-master\include\dephem\EphemerisRelease.cpp"

g++ -std=c++17 %INC% %SRC% "%DEP%" "%MAIN%" "%SOFALIB%" -o "%OUT%"
if errorlevel 1 (
  echo BUILD FAILED.
  exit /b 1
)
echo   OK -^> %OUT%

echo Running:
"%OUT%"
endlocal
