# Сборка модели задержки (ARIADNA port)

Проект — перенос модели геометрической задержки VLBI ARIADNA (Fortran → C++).
Сборка **модульная**: единого CMake нет, каждый модуль/тест компилируется командой
`g++` (см. `tests/test compile commands.txt`). Ниже — как собрать зависимости,
«итоговый билд» (полный конвейер задержки) и как гонять тесты.

## 1. Требования

- **Компилятор** C++17 (проверено на MinGW-w64 g++ под Windows; Git Bash для shell).
- Вендоренные зависимости уже лежат в `external/` (в git):
  - `external/eigen` — линейная алгебра (header-only).
  - `external/dephem-master` — чтение эфемерид JPL DE.
  - `external/sofa/20190722/c` — IAU SOFA (прецессия/нутация/вращение, TDB).
  - `external/spline.h` — кубический сплайн (интерполяция EOP).
  - `external/catalogs/` — каталоги: EOPC04, океаническая/атмосферная нагрузка, antenna-info.
- **Файл эфемерид** DE440 (бинарный dephem): `external/dephem-master/linux_p1550p2650.440t`.
  Тесты его ищут; без него делают SKIP (rc=2).

## 2. Одноразово: статическая библиотека SOFA

Цепочка вращения (ITRF↔GCRS) и TDB строятся на SOFA. Собрать `libsofa.a`:

```bash
sh external/sofa/build_libsofa.sh
# -> external/sofa/20190722/c/build/libsofa.a
```

**ВАЖНО:** при линковке `libsofa.a` указывать ПОСЛЕДНИМ (после всех `.cpp`),
иначе undefined reference на `iauC2t06a`/`iauGst06a`/`iauDtdb`.

## 3. Итоговый билд — полный конвейер задержки (process_ariadna)

`process_ariadna` — верхний слой: читает каталоги, гоняет цикл по наблюдениям и на
каждое считает теоретическую задержку через `compute_delay_obs`
(site_pair → baseline → aber_source → trop_delay → mount_tel → theor_delay).
Запускать из **корня репозитория** (нужны эфемериды и каталоги по относительным путям).

```bash
g++ -std=c++17 \
  -I./external/eigen -I./external -I./external/dephem-master/include \
  -I./external/sofa/20190722/c/src \
  src/process_ariadna.cpp src/process_obs.cpp \
  src/site_pair.cpp src/site_tide_solid.cpp src/site_tide_OC.cpp src/pole_tide.cpp \
  src/site_atm40.cpp src/therm_def40.cpp src/site_inst.cpp \
  src/baseline.cpp src/aber_source.cpp \
  src/trop_delay.cpp src/nhmf2.cpp src/nwmf2.cpp src/sast_dry.cpp src/sast_wet.cpp \
  src/mount_tel.cpp src/sbend.cpp src/theor_delay.cpp \
  src/jpleph.cpp src/ephemeris.cpp src/fund_arg.cpp src/GEOD.cpp \
  src/site_functions.cpp src/rotation.cpp \
  src/dmeteo.cpp src/t_eph40.cpp src/tai_time40.cpp src/nsec.cpp \
  src/interp_eop.cpp src/terms_71.cpp src/terms_lib.cpp src/UT1R_2010.cpp \
  src/READ_CAT.cpp src/catalog_bridge.cpp \
  ./external/dephem-master/include/dephem/EphemerisRelease.cpp \
  tests/test_final_build.cpp \
  ./external/sofa/20190722/c/build/libsofa.a \
  -o tests/test_final_build.exe

./tests/test_final_build.exe   # запуск из корня репозитория
```

Замена `tests/test_final_build.cpp` на свой `main` даёт целевой исполняемый файл.

## 4. Модульные тесты

Полный список команд сборки каждого модуля/теста — в
[`tests/test compile commands.txt`](tests/test%20compile%20commands.txt).
Все тесты самодостаточны (сами ищут файлы) и печатают `РЕЗУЛЬТАТ: PASS/FAIL` (rc 0/1;
rc 2 = SKIP, нет эфемерид). Ключевые:

| Тест | Что проверяет |
|------|---------------|
| `test_final_build`     | Итоговый билд: полный конвейер vs эталон ARIADNA (Вывод.txt) |
| `test_process_ariadna` | Оркестрация: site+EOP(INTERP_EOP40)+нагрузки из каталогов → задержка |
| `test_compute_delay`   | Полный конвейер одного наблюдения vs дамп (2e-12 с) |
| `test_interp_eop`      | EOP: INTERP_EOP40 vs дамп ARIADNA (MJD 57402), ~1e-13 |
| `test_therm_def40`     | Термодеформация (Nothnagel 2009) |
| `test_dmeteo`          | Производные метео (МНК-полином) |
| `test_theor_delay`     | Реляти­вистская задержка vs дамп (покомпонентно) |

## 5. Замечания

- **Точность obs.utc:** UTC наблюдения нужен с ПОЛНОЙ точностью — округление теряет
  секунды → поворот Земли → ошибка задержки в километры.
- **Цепочка вращения и TDB — на SOFA** (стандарт IAU 2006/2000A), не ариадновская
  реализация Жарова (её исходников в архиве нет).
- **Космический VLBI (наземно-космические базы) НЕ поддержан:** нужны
  `THEOR_DELAYcorr_ORB` (задержка орбитальной станции) и `INTEGR8asc` (интегрирование
  орбиты) — не перенесены.
- Каталоги и эфемериды не переупакованы в один архив; пути — относительно корня репо.
