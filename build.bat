@echo off
rem ============================================================================
rem build.bat - сборка под Windows (cmd/PowerShell, компилятор MinGW g++/gcc/ar).
rem
rem Запуск из корня репозитория:
rem   .\build.bat                       собрать libsofa.a (если нет) + verify.exe и запустить
rem   .\build.bat tests\test_xxx.cpp    собрать другой main тем же набором модулей
rem
rem Требуется MinGW в PATH (g++, gcc, ar). Проверить: g++ --version
rem ============================================================================
setlocal enabledelayedexpansion
cd /d "%~dp0"

set "MAIN=%~1"
if "%MAIN%"=="" set "MAIN=tests\verify.cpp"
for %%F in ("%MAIN%") do set "OUT=tests\%%~nF.exe"

set "SOFADIR=external\sofa\20190722\c"
set "SOFALIB=%SOFADIR%\build\libsofa.a"

rem --- 1) libsofa.a (один раз) ---
if not exist "%SOFALIB%" (
  echo [1/2] Собираю libsofa.a ...
  if not exist "%SOFADIR%\build" mkdir "%SOFADIR%\build"
  for %%f in ("%SOFADIR%\src\*.c") do (
    if /I not "%%~nxf"=="t_sofa_c.c" gcc -O2 -c -I"%SOFADIR%\src" "%%f" -o "%SOFADIR%\build\%%~nf.o"
  )
  if exist "%SOFALIB%" del "%SOFALIB%"
  for %%o in ("%SOFADIR%\build\*.o") do ar rs "%SOFALIB%" "%%o" >nul
  echo   libsofa.a готова.
) else (
  echo [1/2] libsofa.a уже есть, пропускаю.
)

rem --- 2) итоговый билд ---
echo [2/2] Собираю %OUT% ...
set "INC=-I.\external\eigen -I.\external -I.\external\dephem-master\include -I.\external\sofa\20190722\c\src"
set "SRC=src\process_ariadna.cpp src\process_obs.cpp src\site_pair.cpp src\site_tide_solid.cpp src\site_tide_OC.cpp src\pole_tide.cpp src\site_atm40.cpp src\therm_def40.cpp src\site_inst.cpp src\baseline.cpp src\aber_source.cpp src\trop_delay.cpp src\nhmf2.cpp src\nwmf2.cpp src\sast_dry.cpp src\sast_wet.cpp src\mount_tel.cpp src\sbend.cpp src\theor_delay.cpp src\jpleph.cpp src\ephemeris.cpp src\fund_arg.cpp src\GEOD.cpp src\site_functions.cpp src\rotation.cpp src\dmeteo.cpp src\orbit_interp.cpp src\t_eph40.cpp src\tai_time40.cpp src\nsec.cpp src\interp_eop.cpp src\terms_71.cpp src\terms_lib.cpp src\UT1R_2010.cpp src\READ_CAT.cpp src\catalog_bridge.cpp"
set "DEP=.\external\dephem-master\include\dephem\EphemerisRelease.cpp"

g++ -std=c++17 %INC% %SRC% "%DEP%" "%MAIN%" "%SOFALIB%" -o "%OUT%"
if errorlevel 1 (
  echo СБОРКА НЕ УДАЛАСЬ.
  exit /b 1
)
echo   OK -^> %OUT%

echo Запуск:
"%OUT%"
endlocal
