# Interferometer Delay Model

Модель геометрической задержки сигнала для РСДБ-коррелятора (перенос ARIADNA с Fortran на
C++). По заданию коррелятора `.cfx` и орбите космического телескопа `.scf` считает
**полиномы задержки** и **полиномы координат `u,v,w`** всех станций сеанса и пишет их в
файлы формата эталона.

Учитывается: геометрическая задержка в вакууме + релятивизм (эффект Шапиро от Солнца, Луны и Земли, аберрация,
члены для наземно-космических баз), твердотельные приливы и прилив полюса, океаническая
и атмосферная нагрузка, термодеформация антенн, тропосфера (стандартная атмосфера),
а для космического телескопа — релятивистский **ход бортовых часов** относительно наземных.
Расчёт **посканный**: у каждого скана свой источник; станция считается только в тех сканах,
где она участвует. Для космического телескопа пересчитываются **TIMEOFS** (задержка сброса
сигнала на пункт приёма) и пишется новое задание `<cfx>_p.cfx`.

## Что нужно (один раз)

1. **MinGW-w64** в `PATH` (нужны `g++`, `gcc`, `ar`). Проверка: `g++ --version`.
2. **Файл эфемерид** JPL DE440 `linux_p1550p2650.440t` в `external/dephem-master/`
   (не в git из-за размера ~112 МБ). Каталоги (EOP, нагрузки, ITRF) уже лежат в
   `external/catalogs/`.

Подробности сборки, состав модулей и диагностика — в [BUILD.md](BUILD.md).

## Сборка

**Windows, PowerShell** (основной путь, в т.ч. терминал VS Code):
```powershell
powershell -ExecutionPolicy Bypass -File .\build.ps1
```
**Windows, cmd:** `.\build.bat`  **Linux / Git Bash:** `sh build.sh`

Соберётся `delay_tool.exe` (при первом запуске также один раз соберётся `libsofa.a`).

## Использование

```
delay_tool <cfx> [scf] [out_dir] [block_sec=60] [degree=5] [tropo=1] [recv=PUSHCH22]
```

| аргумент | смысл |
|---|---|
| `<cfx>` | файл задания коррелятора |
| `[scf]` | орбита космического телескопа (CCSDS OEM). **Необязателен**: если не указан (или `-`), путь берётся из cfx (`ORB_FILE`); указывайте, чтобы подставить свой файл (напр. при проверке на своём компьютере) |
| `[out_dir]` | каталог для выходных файлов. **Необязателен**: если не указан (или `-`), берётся из cfx (`%W`), а имена файлов — из cfx (`POLY_FILE`); если у станции имя не задано, генерируется `<станция>_<диапазон>.txt` |
| `block_sec` | длина блока полинома, с (по умолчанию 60) |
| `degree` | степень полинома (по умолчанию 5 → коэффициенты P0..P5) |
| `tropo` | `1` — с тропосферой (по умолчанию), `0` — только геометрия в вакууме |
| `recv` | имя пункта приёма в `ITRF2005_2.CAT`. По умолчанию `auto` — определяется по префиксу файлов данных космоса (`PUSH…`→`PUSHCH22`, `GBT…`→`GBT_VLBA`); можно задать явно |

Эфемериды и каталоги находятся автоматически (относительно текущего каталога и
родительских уровней), запускать можно из любой папки.

**Пример** (данные лежат в `example/`) — со своим файлом орбиты:
```powershell
.\delay_tool.exe example\GVLBI_RAKS01A0_L_20140423T130000_ASC_V1.cfx `
                 example\RA140423_1200_v02.scf  out_poly
```

**На корреляторе** можно не указывать ни орбиту, ни каталог — всё возьмётся из cfx
(`ORB_FILE` и `%W`):
```powershell
.\delay_tool.exe task.cfx
```

На выходе для каждой станции — два файла: полиномы задержки (`<станция>.txt`) и полиномы
координат `u,v,w` (`<станция>_uvw.txt`, строка `Pk = u, v, w`, в секундах — коррелятор
кладёт их на uv-плоскость). При наличии космоса дополнительно пишется `<cfx>_p.cfx`
с пересчитанными TIMEOFS.

## Проверка (сквозной тест)

```powershell
powershell -ExecutionPolicy Bypass -File .\test.ps1
```
(cmd: `.\test.bat`; Linux/Git Bash: `sh test.sh`.)

Собирает и запускает `test_process_task`: считает полиномы всех станций примера и сверяет
P0 каждого блока с эталонами `example/*.TXT`. Ожидаемо PASS, максимальная невязка ~1–4 м
по всем станциям (включая космический RASTRON). Код возврата 0 = PASS.

## Встраивание в свой код (библиотека C/C++)

Модель — это C++-библиотека в пространстве имён `ariadna`. Весь публичный API объявлен в
**`src/functions.h`** (структуры — в `src/structures.h`, константы — в `src/constants.h`).
Чтобы встроить, ничего в модели менять не нужно: подключаешь заголовок, собираешь свой код
со всеми модулями `src/*.cpp` (+ SOFA + dephem) и вызываешь функции.

**Точки входа (от высокоуровневой к низкой):**

| функция | что делает |
|---|---|
| `init_ephemeris(path)` | загрузить эфемериды JPL — **один раз** за процесс, до остальных вызовов |
| `process_task(cfx, scf, out_dir, eop, block, degree, sample, tropo, recv)` | по заданию + орбите посчитать полиномы задержки и `u,v,w` **всех** станций и записать файлы (то же, что делает `delay_tool`) |
| `compute_task_polys(cfx, scf, eop, block, degree, sample, tropo, out)` | то же **в памяти, без файлов** → `TaskPolys{delay, uvw}` (векторы `StationPoly`/`StationUvw` по станциям). Отдельный модуль-адаптер [`src/delay_api.h`](src/delay_api.h); на нём же построен Python-модуль |
| `compute_station_poly(st, src, …, &uvw)` | полиномы **одной** станции в памяти → `StationPoly` + `StationUvw`; коэффициенты в `poly.blocks[b].coef[k]`, `uvw.blocks[b].u/v/w[k]` |
| `process_ariadna(…, results)` / `compute_delay_obs(…)` | сами **значения задержки** (не полиномы) на заданные моменты → `std::vector<DelayResult>` / одно наблюдение |

То есть можно получать **и полиномы, и сами задержки — по потребности**. Рабочий пример —
[`examples/use_as_library.cpp`](examples/use_as_library.cpp).

**Сборка своего кода.** Внутри репозитория build-скрипт собирает любой main с полным набором
модулей:
```powershell
powershell -ExecutionPolicy Bypass -File .\build.ps1 examples\use_as_library.cpp
```
Для **внешнего** проекта соберите статическую библиотеку из всех модулей `src/*.cpp` и
линкуйте её вместе с `libsofa.a`:
```sh
g++ -std=c++17 -c -I./external/eigen -I./external -I./external/dephem-master/include \
    -I./external/sofa/20190722/c/src \
    src/*.cpp external/dephem-master/include/dephem/EphemerisRelease.cpp
ar rcs libariadna.a *.o
# ваша программа:
g++ -std=c++17 -I<пути выше> your_program.cpp libariadna.a \
    external/sofa/20190722/c/build/libsofa.a -o your_program
```

## Использование из Python

Есть **готовый нативный модуль** `ariadna` (pybind11) поверх адаптера
[`compute_task_polys`](src/delay_api.h) — считает полиномы **в памяти, без файлов**, и
отдаёт их прямо в Python. Отдельный модуль-адаптер (`src/delay_api.cpp`) не меняет остальной
код модели.

**Сборка модуля** (нужны `g++`, `pip install pybind11`, Python с dev-заголовками):
```powershell
powershell -ExecutionPolicy Bypass -File .\python\build_pybind.ps1   # Windows
sh python/build_pybind.sh                                            # Linux / Git Bash
```
Соберётся `ariadna.pyd` (на Linux — `ariadna.so`) в корне репозитория.

**Использование:**
```python
import ariadna
ariadna.init_ephemeris("external/dephem-master/linux_p1550p2650.440t")   # один раз
res = ariadna.compute_task_polys("task.cfx", "orbit.scf",               # orbit="" -> из cfx
        eop="external/catalogs/EOPC04_14_IAU2000_62-now.cat",
        block=60, degree=5, sample=6.0, tropo=True)
for sp, uv in zip(res.delay, res.uvw):       # по станциям
    for b in sp.blocks:                      # блоки полинома задержки
        b.mjd, b.utc_start, b.source, b.coef # coef = [P0..P5] (список float -> np.array)
    for b in uv.blocks:                      # блоки координат u,v,w
        b.u, b.v, b.w                        # по осям, списки [P0..P5]
```
Полный пример — [`python/example.py`](python/example.py) (`python python/example.py`).
Так можно получать в Python **и полиномы, и сами задержки** (значение полинома в момент `t`:
`sum(P_k * t**k)`).

**Альтернативы без сборки модуля:**
- **CLI + `subprocess`** — вызвать `delay_tool` и прочитать `<станция>.txt` / `*_uvw.txt`
  (парсинг: `telescope/order`, затем блоки `source/start/stop/P0..P5`). Работает всегда.
- **ctypes/cffi** — потребует `extern "C"` слоя; менее удобно, чем готовый `ariadna`.

Где полезно: пакетная обработка заданий, построение uv-покрытия (u,v из uvw), сверка
с эталоном, встраивание в питоновый пайплайн планировщика/коррелятора.
