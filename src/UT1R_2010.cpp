#include "functions.h"

namespace ariadna {

void ut1r_2010(const Eigen::VectorXd& f, double& dut, double& dlod, double& domega) {
    // Проверка входного массива
    if (f.size() != 5) {
        throw std::invalid_argument("Input vector f must have 5 elements (l, l', F, D, Omega)");
    }

    // Инициализация выходных параметров
    dut = 0.0;
    dlod = 0.0;
    domega = 0.0;

    // Цикл по приливным членам
    for (int j = 0; j < cnst::N_TIDAL_TERMS; ++j) {
        // Вычисление аргумента (в радианах)
        double arg = 0.0;
        for (int i = 0; i < 5; ++i) {
            arg += cnst::TIDAL_COEFFS[j][i] * f(i);
        }
        arg = std::fmod(arg, cnst::SEC360) * cnst::CARCRAD;

        // Суммирование поправок
        dut += cnst::TIDAL_COEFFS[j][5] * std::sin(arg) + cnst::TIDAL_COEFFS[j][6] * std::cos(arg);
        dlod += cnst::TIDAL_COEFFS[j][7] * std::cos(arg) + cnst::TIDAL_COEFFS[j][8] * std::sin(arg);
        domega += cnst::TIDAL_COEFFS[j][9] * std::cos(arg) + cnst::TIDAL_COEFFS[j][10] * std::sin(arg);
    }

    // Масштабирование результатов
    dut *= cnst::UT1_SCALE;
    dlod *= cnst::LOD_SCALE;
    domega *= cnst::OMEGA_SCALE;
}

} // namespace ariadna