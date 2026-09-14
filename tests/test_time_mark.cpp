// test_time_mark.cpp
//
// Проверка точки расширения для метки времени TIMEOFS (см. functions.h, раздел
// "Метка времени для TIMEOFS"). Проверяются три состояния переключателя:
//
//   1) по умолчанию   — источник .cfx, метка разбирается из строки FILExx (YYYYDDDHHMMSS);
//   2) выключен, без обработчика — работает встроенная заглушка, метки нет;
//   3) выключен, с обработчиком  — метку выдаёт посторонний код, и она доходит до вызова.
//
// Файлы данных и эфемериды не нужны: проверяется только маршрутизация запроса.

#include "../src/functions.h"
#include <cstdio>
#include <cmath>

using namespace ariadna;

static int fails = 0;

static void check(bool ok, const char* what) {
    std::printf("  [%s] %s\n", ok ? " ok " : "ОШИБ", what);
    if (!ok) ++fails;
}

// Посторонний код: выдаёт заранее известную метку и запоминает, что его позвали.
static bool g_hook_called = false;
static TimeMarkRequest g_seen;
static const double HOOK_MJD = 56770.5;

static bool test_hook(const TimeMarkRequest& req, double& mjd_utc) {
    g_hook_called = true;
    g_seen = req;
    mjd_utc = HOOK_MJD;
    return true;
}

int main() {
    std::printf("test_time_mark: точка расширения для метки времени TIMEOFS\n");

    TimeMarkRequest req;
    req.cfx_path   = "example/task.cfx";
    req.station    = "RASTRON";
    req.file_index = "00";
    // 2014, 113-й день года, 13:00:00 UTC -> MJD 56770.541666...
    req.file_value = " %P:RA_2014113130000.rdf";

    // MJD для 2014-04-23 = 56770 (23 апреля — 113-й день невисокосного 2014 года).
    const double expect_cfx = 56770.0 + 13.0 / 24.0;

    // --- 1) состояние по умолчанию: источник .cfx
    check(time_mark_from_cfx(), "по умолчанию переключатель стоит на .cfx");
    double mjd = 0.0;
    bool got = get_time_mark(req, mjd);
    check(got, "метка из .cfx получена");
    check(got && std::fabs(mjd - expect_cfx) < 1e-9, "метка из .cfx разобрана верно");
    if (got) std::printf("        MJD = %.9f (ожидалось %.9f)\n", mjd, expect_cfx);

    // --- 2) переключатель выключен, обработчик не установлен: заглушка
    set_time_mark_from_cfx(false);
    check(!time_mark_from_cfx(), "переключатель выключен");
    mjd = -1.0;
    got = get_time_mark(req, mjd);
    check(!got, "заглушка метку не выдаёт (TIMEOFS не пишется)");
    check(mjd == -1.0, "заглушка не трогает выходное значение");

    // --- 3) переключатель выключен, обработчик установлен: метку даёт посторонний код
    set_time_mark_hook(&test_hook);
    mjd = 0.0;
    got = get_time_mark(req, mjd);
    check(got, "обработчик вызван и метку выдал");
    check(g_hook_called, "обработчик действительно был вызван");
    check(got && std::fabs(mjd - HOOK_MJD) < 1e-12, "метка от обработчика дошла без изменений");
    check(g_seen.station == req.station && g_seen.file_index == req.file_index &&
          g_seen.file_value == req.file_value && g_seen.cfx_path == req.cfx_path,
          "запрос дошёл до обработчика полностью");

    // --- 4) снятие обработчика возвращает заглушку
    set_time_mark_hook(nullptr);
    mjd = -1.0;
    got = get_time_mark(req, mjd);
    check(!got, "set_time_mark_hook(nullptr) возвращает заглушку");

    // --- 5) возврат к .cfx восстанавливает исходное поведение
    set_time_mark_from_cfx(true);
    mjd = 0.0;
    got = get_time_mark(req, mjd);
    check(got && std::fabs(mjd - expect_cfx) < 1e-9, "возврат на .cfx восстанавливает поведение");

    std::printf(fails ? "\nПРОВАЛЕНО проверок: %d\n" : "\nвсе проверки пройдены\n", fails);
    return fails ? 1 : 0;
}
