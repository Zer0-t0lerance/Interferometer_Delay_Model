#!/bin/sh
# Собирает статическую библиотеку SOFA (libsofa.a) из исходников c/src.
# Запуск (один раз):  sh external/sofa/build_libsofa.sh
# Результат: external/sofa/20190722/c/build/libsofa.a
set -e
cd "$(dirname "$0")/20190722/c"
mkdir -p build
for f in src/*.c; do
  case "$f" in
    */t_sofa_c.c) continue ;;   # пропускаем собственный тест SOFA
  esac
  gcc -O2 -c -Isrc "$f" -o "build/$(basename "${f%.c}").o"
done
ar rcs build/libsofa.a build/*.o
echo "built: $(pwd)/build/libsofa.a"
