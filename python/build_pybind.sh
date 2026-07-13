#!/bin/sh
# build_pybind.sh — сборка Python-модуля ariadna (pybind11) -> ariadna.pyd в корне репозитория.
#
# Требуется: g++ (MinGW) в PATH; pybind11 (pip install pybind11); Python с dev-заголовками.
# Пути к Python берутся автоматически из текущего `python`. libsofa.a должна быть собрана
# (её собирает build.sh/ps1 при первом запуске).
#
# Запуск (из любой папки):
#   sh python/build_pybind.sh
# Затем из корня репозитория:  python -c "import ariadna; print(ariadna.__doc__)"

set -e
cd "$(dirname "$0")/.."

PYBIND_INC=$(python -c "import pybind11; print(pybind11.get_include())")
PY_INC=$(python -c "import sysconfig; print(sysconfig.get_path('include'))")
PY_LIBDIR=$(python -c "import sys,os; print(os.path.join(sys.base_prefix,'libs'))")
PY_LIB=$(python -c "import sysconfig; print('python'+sysconfig.get_config_var('py_version_nodot'))")

INC="-I./external/eigen -I./external -I./external/dephem-master/include -I./external/sofa/20190722/c/src -I$PYBIND_INC -I$PY_INC"
DEP="./external/dephem-master/include/dephem/EphemerisRelease.cpp"
SOFA="./external/sofa/20190722/c/build/libsofa.a"

echo "Собираю ariadna.pyd (линкую с $PY_LIB) ..."
g++ -O2 -shared -std=c++17 $INC \
    python/ariadna_pybind.cpp src/*.cpp "$DEP" "$SOFA" \
    -L"$PY_LIBDIR" -l"$PY_LIB" -o ariadna.pyd
echo "OK -> ariadna.pyd"
echo "Проверка:  python python/example.py"
