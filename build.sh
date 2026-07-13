#!/bin/sh
# build.sh — сборка итогового конвейера задержки (process_ariadna) с тестом.
#
# Запуск (из корня репозитория, в Git Bash / sh):
#   sh build.sh                      # собрать delay_tool.exe и запустить (покажет справку)
#   sh build.sh tests/test_xxx.cpp   # собрать другой main/тест с тем же набором модулей
#
# ВНИМАНИЕ: это sh-скрипт. НЕ вставляйте его тело построчно в PowerShell/cmd —
# там перенос строки '\' не работает и линкер ругается "cannot find \".
# В PowerShell/cmd запускайте через:  sh build.sh
#
# Требуется заранее собранная libsofa.a:  sh external/sofa/build_libsofa.sh

set -e

MAIN="${1:-delay_tool.cpp}"
OUT="${MAIN%.cpp}.exe"

INC="-I./external/eigen -I./external -I./external/dephem-master/include -I./external/sofa/20190722/c/src"
DEP="./external/dephem-master/include/dephem/EphemerisRelease.cpp"
SOFA="./external/sofa/20190722/c/build/libsofa.a"

SRC="src/process_ariadna.cpp src/process_obs.cpp \
src/site_pair.cpp src/site_tide_solid.cpp src/site_tide_OC.cpp src/pole_tide.cpp \
src/site_atm40.cpp src/therm_def40.cpp src/site_inst.cpp \
src/baseline.cpp src/aber_source.cpp \
src/trop_delay.cpp src/nhmf2.cpp src/nwmf2.cpp src/sast_dry.cpp src/sast_wet.cpp \
src/mount_tel.cpp src/sbend.cpp src/theor_delay.cpp \
src/jpleph.cpp src/ephemeris.cpp src/fund_arg.cpp src/GEOD.cpp \
src/site_functions.cpp src/rotation.cpp \
src/dmeteo.cpp src/orbit_interp.cpp src/t_eph40.cpp src/tai_time40.cpp src/nsec.cpp \
src/interp_eop.cpp src/terms_71.cpp src/terms_lib.cpp src/UT1R_2010.cpp \
src/READ_CAT.cpp src/catalog_bridge.cpp \
src/delay_poly.cpp src/cfx_parser.cpp src/scf_reader.cpp src/delay_api.cpp"

echo "Собираю $OUT ..."
g++ -std=c++17 $INC $SRC "$DEP" "$MAIN" "$SOFA" -o "$OUT"
echo "OK -> $OUT"

echo "Запуск:"
./"$OUT"
