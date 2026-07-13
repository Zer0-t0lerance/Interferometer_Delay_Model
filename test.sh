#!/bin/sh
# test.sh — сквозной тест готового модуля (Git Bash / Linux).
#
# Собирает и запускает test_process_task: process_task() по example/ (.cfx + .scf)
# считает полиномы задержки ВСЕХ станций сеанса и сверяет P0 каждого блока с эталонами
# example/*.TXT. Печатает ошибку по станциям и PASS/FAIL.
#
# Запуск из корня репозитория:
#   sh test.sh
#
# Требуются g++ в PATH и файл эфемерид в external/dephem-master/ (см. BUILD.md).
# Код возврата 0 = PASS, иначе FAIL.

sh build.sh tests/test_process_task.cpp
code=$?
echo ""
if [ "$code" -eq 0 ]; then echo "TEST PASSED (exit 0)."; else echo "TEST FAILED (exit $code)."; fi
exit $code
