# Сборка и запуск модели задержки (ARIADNA port)

Перенос модели геометрической задержки VLBI ARIADNA (Fortran → C++). Сборка простыми
командами `g++` — **единого CMake нет**. Ниже: что нужно, иерархия, как собрать под
Windows и Linux, как запускать и как добавлять свой код.

---

## 1. Что нужно (и что взять отдельно)

**Приходит с репозиторием (git):**
- исходники `src/`, тесты `tests/`;
- вендоренные зависимости `external/`: Eigen, dephem, SOFA (исходники), spline.h;
- каталоги `external/catalogs/`: EOPC04, океаническая/атмосферная нагрузка, antenna-info.

**НЕ в git — скопировать вручную (иначе тесты делают SKIP):**
- **Файл эфемерид JPL DE440** `linux_p1550p2650.440t` (~112 МБ) — положить в
  `external/dephem-master/`. (В `.gitignore` стоит `*.440t`.) Без него расчёт эфемерид
  недоступен.

**Компилятор:** MinGW-w64 (`g++`, `gcc`, `ar`) в `PATH`. Проверка: `g++ --version`.

---

## 2. Иерархия проекта

```
delay_model/
├─ build.bat            ← сборка под Windows (cmd/PowerShell)
├─ build.sh             ← сборка под Linux / Git Bash
├─ BUILD.md
├─ .vscode/tasks.json   ← сборка из VS Code (Ctrl+Shift+B)
├─ src/                 ← все модули модели (*.cpp, *.h)
├─ tests/               ← тесты и verify.cpp; Вывод.txt (эталон ARIADNA)
└─ external/
   ├─ eigen/                         (header-only)
   ├─ dephem-master/                 (+ сюда linux_p1550p2650.440t вручную)
   ├─ sofa/20190722/c/{src,build}/   (build/libsofa.a собирается локально)
   └─ catalogs/                      (EOPC04, ocload, atmload, antenna-info)
```

---

## 3. Сборка под Windows (PowerShell / VS Code / cmd)

Из корня репозитория. **В терминале VS Code (PowerShell) — предпочтительно `build.ps1`:**

```powershell
powershell -ExecutionPolicy Bypass -File .\build.ps1
powershell -ExecutionPolicy Bypass -File .\build.ps1 tests\test_xxx.cpp
```
(если локальные скрипты разрешены — просто `.\build.ps1`).

**Из cmd — `build.bat`:**
```bat
.\build.bat
.\build.bat tests\test_xxx.cpp
```

Оба скрипта: 1) один раз собирают `libsofa.a` из исходников SOFA; 2) собирают
переданный `main` (по умолчанию `tests\verify.cpp`) со всеми ~35 модулями; 3) запускают.

**VS Code:** `Ctrl+Shift+B` → задача **«build (Windows)»** (запускает `build.ps1`).
Есть «build current file (Windows)» — собрать открытый `tests/*.cpp`.

> **Кодировка (кракозябры):** программы сами ставят консоль в UTF-8 (`SetConsoleOutputCP`),
> поэтому русский вывод читается. Отдельно всё пишется в `tests/verify_report.txt`
> (UTF-8, открыть в VS Code). Скрипты сборки — ASCII (иначе PowerShell 5.1 / cmd ломают
> разбор). Файлы `.bat` — CRLF (через `.gitattributes`), иначе cmd не разбирает.
> **НЕ** вставляйте команду `g++` построчно с `\` в PowerShell/cmd — там `\` не
> продолжение строки → `ld.exe: cannot find \`. Используйте скрипт.

---

## 4. Сборка под Linux / Git Bash

```sh
sh external/sofa/build_libsofa.sh   # один раз: собрать libsofa.a
sh build.sh                         # собрать tests/verify.exe и запустить
sh build.sh tests/test_xxx.cpp      # другой main
```

---

## 5. Запуск

Exe запускать **из корня репозитория** (модель ищет эфемериды и каталоги по
относительным путям). Бинарники (`*.exe`, `*.a`) в git не кладутся.

Главные программы:
- **`tests/verify.exe`** — ПОДРОБНАЯ сверка с дампом ARIADNA `tests/Вывод.txt`:
  таблица «величина | референс | посчитано | ошибка» по всему циклу (элевации,
  зенитные/тропосферные задержки, итог). Отчёт в консоль И в `tests/verify_report.txt`.
- **`tests/test_final_build.exe`** — итоговый прогон конвейера (8 станций из Вывод.txt).

---

## 6. Как добавить свой код

Модель — библиотека функций в `namespace ariadna` (объявления — `src/functions.h`,
структуры — `src/structures.h`, константы — `src/constants.h`). Свой `main`:

1. создайте `tests/my_run.cpp`, включите `#include "../src/functions.h"`;
2. соберите: `.\build.bat tests\my_run.cpp` (Windows) или `sh build.sh tests/my_run.cpp`.

Верхняя функция — `process_ariadna(...)`: принимает станции/источники/наблюдения/EOP,
**возвращает результаты В ПАМЯТИ** (`std::vector<DelayResult>`), опционально —
промежуточные величины (`std::vector<CompDebug>`). Файл-отчёт пишется ТОЛЬКО если
передать непустой `output_path` (по умолчанию — без файлов).

> **Без промежуточных файлов.** Данные между стадиями конвейера передаются в памяти
> (не через временные файлы). Файлы — только для итоговых человекочитаемых отчётов.

---

## 7. Тесты и замечания

Команды сборки каждого модульного теста — в
[`tests/test compile commands.txt`](tests/test%20compile%20commands.txt).
Тесты печатают `РЕЗУЛЬТАТ: PASS/FAIL` (rc 0/1; rc 2 = SKIP, нет эфемерид).

- **Точность obs.utc:** UTC наблюдения нужен с ПОЛНОЙ точностью — округление теряет
  секунды → поворот Земли → ошибка задержки в километры.
- **Вращение и TDB — на SOFA** (стандарт IAU 2006/2000A), не ариадновская реализация.
- **EOP — INTERP_EOP40** (сверен с дампом ARIADNA на MJD 57402 до ~1e-13).
- **Космический VLBI:** путь есть (флаг `is_space`, `orbit_interp` по готовой орбите),
  но нет орбитальной эфемериды `.scf` для сквозной валидации.
