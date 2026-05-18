#pragma once
#include "structures.h"

namespace ariadna {

// ============================================================================
// Конвертация времени и параметров ориентации Земли (EOP)
// ============================================================================

/**
 * @brief Конвертация времени из шкал MJD и UTC в атомное (TAI) и земное время (TT).
 */
void tai_time(double mjd, double UTC, double &TAI, double &TT);

/**
 * @brief Определение количества високосных секунд (leap seconds) на заданную дату.
 */
void nsec(double mjd, double& idelt);

/**
 * @brief Интерполяция данных параметров ориентации Земли (IERS) для заданного наблюдения.
 */
void interp_iers(const Observation& obs, Eigen::VectorXd& eop_int);

/**
 * @brief Вычисляет фундаментальные аргументы Луны и Солнца (IERS 2003).
 *
 * @param[in]  jd    Юлианская дата (дни).
 * @param[in]  ct    Координатное время (доля дня).
 * @param[out] cent  Юлианские столетия с J2000.0 (T).
 * @param[out] f     5 аргументов (l, l', F, D, Omega) [арксекунды].
 * @param[out] fd    Производные аргументов [арксекунды/столетие].
 */
void fund_arg(double jd, double ct, double& cent, Eigen::VectorXd& f, Eigen::VectorXd& fd);

/**
 * @brief Интерполяция параметров ориентации Земли (EOP) на момент наблюдения.
 * * Вычисляет UT1-UTC, координаты полюса (x, y) и углы нутации (dpsi, deps), 
 * а также их производные. Включает поправки за зональные приливы, суточные 
 * приливы и либрацию.
 *
 * @param[in]  k_int       Метод: 0 - кубический сплайн, 1 - полином Лагранжа.
 * @param[in]  obs         Структура наблюдения (содержит MJD и UTC).
 * @param[in]  tt          Земное время (TT) в долях суток.
 * @param[out] ut1         Вычисленное время UT1 (сутки).
 * @param[out] eop_int     Вектор интерполированных EOP: [UT1-UTC, x, y, dpsi, deps].
 * @param[out] deop_int    Производные EOP по времени.
 * @param[out] arg_oc_tide Аргументы океанических приливов (8 значений).
 * @param[out] deop_diu    Суточные приливные поправки (матрица 3x2).
 * @param[out] deop_lib    Поправки за либрацию (матрица 3x2).
 * @param[in]  eop_data    Вектор опорных данных EOP (7 точек вокруг момента).
 */
void interp_eop(int k_int, const Observation& obs, double tt, double& ut1,
                Eigen::VectorXd& eop_int, Eigen::VectorXd& deop_int,
                Eigen::VectorXd& arg_oc_tide, Eigen::MatrixXd& deop_diu,
                Eigen::MatrixXd& deop_lib, const std::vector<EOPData>& eop_data);

/**
 * @brief Вычисляет поправки за зональные приливы в ротации Земли (IERS 2010).
 * * Рассчитывает вариации UT1-UT1R, изменения длительности суток (LOD) и 
 * угловой скорости (Omega) для 62 гармоник с периодами от 5 дней до 18.6 лет.
 *
 * @param[in]  f       Вектор 5 фундаментальных аргументов (l, l', F, D, Omega) [арксекунды].
 * @param[out] dut     Поправка UT1-UT1R [секунды].
 * @param[out] dlod    Поправка длительности суток (LOD) [секунды].
 * @param[out] domega  Поправка угловой скорости [рад/сек].
 */
void ut1r_2010(const Eigen::VectorXd& f, double& dut, double& dlod, double& domega);

/**
 * @brief Расчет инструментальной задержки, вызванной смещением осей монтировки телескопа.
 * Учитывает влияние атмосферной рефракции на видимое положение источника.
 * * @param[in]  obs       Структура наблюдения (содержит метеоданные и индексы станций).
 * @param[in]  r2000     Матрица перехода от земной системы (Crust-fixed) к J2000.0 и ее производные (9x3).
 * @param[in]  stations  Вектор данных станций (тип монтировки, смещение, координаты).
 * @param[in]  k_star    Единичный вектор направления на источник в J2000.0.
 * @param[in]  vw        Матрицы перехода от топоцентрической системы (VEN) к геоцентрической (Crust-fixed).
 * @param[in]  e         Матрица углов возвышения (элевации) источника и их производных [рад, рад/с].
 * @param[in]  az        Матрица азимутов источника и их производных [рад, рад/с].
 * @param[out] doff_dl   Частные производные вектора смещения монтировки по величине самого смещения (блок 2x2).
 * @param[out] d_dax     Частные производные задержки и ее скорости по смещению оси [с/м, с/с/м] (блок 2x2).
 * @param[out] dtau_off  Вклад смещения оси монтировки в общую задержку и скорость изменения задержки [с, с/с] (блок 2x2).
 */
void mount_tel(const Observation& obs, const Eigen::MatrixXd& r2000, const std::vector<Station>& stations, const std::vector<Eigen::Vector3d>& k_star, const std::vector<Eigen::Matrix3d>& vw, const Eigen::MatrixXd& e, const Eigen::MatrixXd& az, Eigen::MatrixXd& doff_dl, Eigen::MatrixXd& d_dax, Eigen::MatrixXd& dtau_off);

/**
 * @brief Вычисление угла изгиба луча (атмосферной рефракции).
 * * @param[in] el_rad    Угол возвышения (элевация) наблюдаемого объекта [радианы].
 * @param[in] temp_k    Температура окружающего воздуха [Кельвины].
 * @param[in] humid_f   Относительная влажность воздуха [доли единицы, 0.0 - 1.0].
 * @param[in] press_hg  Атмосферное давление [мм ртутного столба].
 * @return              Угол изгиба луча (рефракции) [радианы].
 */
double sbend(double el_rad, double temp_k, double humid_f, double press_hg);

/**
 * @brief Вычисление гидростатической (сухой) картирующей функции Нилла (NMF).
 * * @param[in]  epoch    Юлианская дата наблюдения (JD).
 * @param[in]  latitude Геодезическая широта станции [радианы].
 * @param[in]  height   Геодезическая высота станции над эллипсоидом [метры].
 * @param[in]  elev     Угол возвышения (элевация) источника [радианы].
 * @param[out] hmf      Вектор (размер 2): [0] - значение картирующей функции, 
 * [1] - её производная по углу возвышения.
 */
void nhmf2(double epoch, double latitude, double height, double elev, Eigen::Vector2d& hmf);

/**
 * @brief Вычисление влажной картирующей функции Нилла (NMF).
 * * @param[in]  latitude Геодезическая широта станции [радианы].
 * @param[in]  elev     Угол возвышения (элевация) источника [радианы].
 * @param[out] wmf      Вектор (размер 2): [0] - значение влажной картирующей функции, 
 * [1] - её производная по углу возвышения.
 */
void nwmf2(double latitude, double elev, Eigen::Vector2d& wmf);

/**
 * @brief Расчет сухой компоненты зенитной тропосферной задержки по модели Саастамойнена.
 * * @param[in]  pres     Атмосферное давление на станции [мбар].
 * @param[in]  dot_pres Скорость изменения давления [мбар/с].
 * @param[in]  lat_geod Геодезическая широта станции [радианы].
 * @param[in]  height   Высота станции [метры].
 * @param[in]  dpdh     Производная давления по высоте.
 * @param[out] z_d      Сухая зенитная задержка [метры].
 * @param[out] dot_z_d  Производная сухой задержки по времени [м/с].
 * @param[out] dz_ddh   Производная сухой задержки по высоте.
 */
void sast_dry(double pres, double dot_pres, double lat_geod, double height, double dpdh, double& z_d, double& dot_z_d, double& dz_ddh);

/**
 * @brief Расчет влажной компоненты зенитной тропосферной задержки по модели Саастамойнена.
 * * @param[in]  rel_hum     Относительная влажность воздуха [%].
 * @param[in]  tc          Температура окружающего воздуха [градусы Цельсия].
 * @param[in]  dot_rel_hum Скорость изменения влажности [%/с].
 * @param[in]  dot_tc      Скорость изменения температуры [Цельсий/с].
 * @param[out] z_w         Влажная зенитная задержка [метры].
 * @param[out] dot_z_w     Производная влажной задержки по времени [м/с].
 */
void sast_wet(double rel_hum, double tc, double dot_rel_hum, double dot_tc, double& z_w, double& dot_z_w);

/**
 * @brief Комплексный расчет тропосферной задержки для двух станций интерферометра.
 * Вычисляет зенитные задержки (сухую и влажную), картирующие функции Нилла (NMF),
 * а также горизонтальные градиенты задержки (Север и Восток). Значения для первой 
 * станции идут со знаком минус (для формирования разности задержек на базе).
 *
 * @param[in]  obs       Структура наблюдения (содержит давление, температуру, влажность).
 * @param[in]  jd        Юлианская дата наблюдения.
 * @param[in]  ct        Координатное время (доля дня).
 * @param[in]  sta1      Данные первой станции (широта, высота).
 * @param[in]  sta2      Данные второй станции (широта, высота).
 * @param[in]  e         Матрица 2x2: элевации [рад] и их производные [рад/с].
 * @param[in]  az        Матрица 2x2: азимуты [рад] и их производные [рад/с].
 * @param[out] datmc_d   Матрица 2x2: сухая тропосферная задержка и её производная [с, с/с].
 * @param[out] datmc_w   Матрица 2x2: влажная тропосферная задержка и её производная [с, с/с].
 * @param[out] datmp_hmf Матрица 2x2: сухая (гидростатическая) картирующая функция и её производная.
 * @param[out] datmp_wmf Матрица 2x2: влажная картирующая функция и её производная.
 * @param[out] dgrad_n   Матрица 2x2: северный градиент картирующей функции и его производная.
 * @param[out] dgrad_e   Матрица 2x2: восточный градиент картирующей функции и его производная.
 * @param[out] zen_dry   Матрица 2x2: сухая зенитная задержка [с] и её производная [с/с].
 * @param[out] zen_wet   Матрица 2x2: влажная зенитная задержка [с] и её производная [с/с].
 */
void trop_delay(const Observation& obs, double jd, double ct, const Station& sta1, const Station& sta2, const Eigen::MatrixXd& e, const Eigen::MatrixXd& az, Eigen::MatrixXd& datmc_d, Eigen::MatrixXd& datmc_w, Eigen::MatrixXd& datmp_hmf, Eigen::MatrixXd& datmp_wmf, Eigen::MatrixXd& dgrad_n, Eigen::MatrixXd& dgrad_e, Eigen::MatrixXd& zen_dry, Eigen::MatrixXd& zen_wet);

/**
 * @brief Вычисление суточных и полусуточных приливных поправок (71 гармоника) к EOP.
 */
void terms_71(double cent, const Eigen::VectorXd& f, const Eigen::VectorXd& fd, Eigen::MatrixXd& dEOP_diu, Eigen::VectorXd& arg_oc_tide);

/**
 * @brief Вычисление поправок к EOP, обусловленных либрацией нежесткой Земли.
 */
void terms_lib(double cent, const Eigen::VectorXd& f, const Eigen::VectorXd& fd, Eigen::MatrixXd& dEOP_lib);


// ============================================================================
// Вычисление координат, систем отсчета и геофизических эффектов
// ============================================================================

/**
 * @brief Инициализация и расчет базовых геоцентрических и геодезических координат для списка станций.
 */
void site(const std::vector<Station>& stations, const std::vector<SpaceStation>& space_stations, int n_stations, double dyear,
          std::vector<Eigen::Vector3d>& site_xyz, std::vector<Eigen::Vector3d>& site_vel, std::vector<double>& lat_geod,
          std::vector<double>& h_geod, std::vector<double>& lat_gcen,
          std::vector<double>& lon_gcen, std::vector<double>& sph_rad,
          std::vector<double>& u_site, std::vector<double>& v_site,
          std::vector<Eigen::Matrix3d>& vw);

/**
 * @brief Вычисляет геоцентрические координаты, скорости и ускорения станций в системе J2000 с учетом всех поправок.
 */
void site(const std::vector<Station>& stations, int j1, int j2, const Observation& observation,
          double cent, const Eigen::VectorXd& f, const Eigen::VectorXd& fd,
          double gast, const Eigen::MatrixXd& r2000_full,
          const Eigen::VectorXd& eop_int, const Eigen::MatrixXd& deop_diu,
          const Eigen::MatrixXd& deop_lib,
          std::vector<Eigen::Vector3d>& xsta_j2000t, std::vector<Eigen::Vector3d>& vsta_j2000t,
          std::vector<Eigen::Vector3d>& asta_j2000);

/**
 * @brief Вычисление теоретической задержки и ее производной (модель IERS 2000/2010).
 * Учитывает релятивистскую гравитационную задержку от Земли, Солнца и Луны (эффект Шапиро),
 * геометрическую часть задержки с поправками на скорости движения станций, 
 * а также прибавляет тропосферные и инструментальные смещения.
 *
 * @param[in]  base_line   Матрица 3x2: вектор базы (от ст.1 к ст.2) и его производная [м, м/с].
 * @param[in]  xsta        Векторы (2): геоцентрические координаты станций 1 и 2 в J2000 [м].
 * @param[in]  vsta        Векторы (2): геоцентрические скорости станций 1 и 2 в J2000 [м/с].
 * @param[in]  asta        Векторы (2): геоцентрические ускорения станций 1 и 2 в J2000 [м/с^2].
 * @param[in]  K_s         Единичный вектор направления на источник в J2000.
 * @param[in]  Earth       Матрица 3x3: координаты, скорость и ускорение барицентра Земли в SSB [м, м/с, м/с^2].
 * @param[in]  Sun         Матрица 3x3: координаты, скорость и ускорение Солнца в SSB [м, м/с, м/с^2].
 * @param[in]  Moon        Матрица 3x3: координаты, скорость и ускорение Луны в SSB [м, м/с, м/с^2].
 * @param[in]  Datmc_d     Матрица 2x2: сухая тропосфера (строка 0 - задержка, строка 1 - производная; столбцы - ст1, ст2) [с, с/с].
 * @param[in]  Datmc_w     Матрица 2x2: влажная тропосфера [с, с/с].
 * @param[in]  dtau_off    Матрица 2x2: инструментальная задержка (оси монтировки) [с, с/с].
 * @param[in]  dt_temp     Матрица 2x2: температурная поправка антенны [с, с/с].
 * @param[out] t2_t1       Полная теоретическая задержка [секунды].
 * @param[out] dt2_t1      Производная теоретической задержки по времени [секунды/секунду].
 */
void theor_delay(const Eigen::Matrix<double, 3, 2>& base_line,
                 const std::vector<Eigen::Vector3d>& xsta,
                 const std::vector<Eigen::Vector3d>& vsta,
                 const std::vector<Eigen::Vector3d>& asta,
                 const Eigen::Vector3d& K_s,
                 const Eigen::Matrix3d& Earth,
                 const Eigen::Matrix3d& Sun,
                 const Eigen::Matrix3d& Moon,
                 const Eigen::Matrix2d& Datmc_d,
                 const Eigen::Matrix2d& Datmc_w,
                 const Eigen::Matrix2d& dtau_off,
                 const Eigen::Matrix2d& dt_temp,
                 double& t2_t1, double& dt2_t1);

/**
 * @brief Вычисляет поправки к координатам станции из-за температурных деформаций монтировки.
 */
void therm_def(const Station& station, const Observation& obs, double dtdt, const Eigen::Matrix3d& vw, const Eigen::MatrixXd& r2000, Eigen::Vector3d& dx_temp, Eigen::Vector3d& dv_temp);

/**
 * @brief Вычисляет смещение и скорость станции, вызванные приливом полюса (Pole Tide) в ITRF->J2000.0.
 */
void POLE_TIDE(double cent, double lat_geod, double lon_geod,
               double xp, double yp, double xp_rate, double yp_rate,
               const Eigen::Matrix3d& vw_i, const Eigen::MatrixXd& r2000,
               Eigen::Vector3d& dx_poltide, Eigen::Vector3d& dv_poltide);

/**
 * @brief Расчет смещения коры и скорости из-за твердых земных приливов (Solid Earth Tides).
 * Вычисляет гравитационное влияние Луны и Солнца по модели IERS 2000.
 * Включает базовый прилив 2 и 3 степени, поправки на неупругость мантии,
 * широтную зависимость, а также частотно-зависимые поправки (суточные и долгопериодические).
 *
 * @param[in]  xsta_itrf Геоцентрические координаты станции в земной системе ITRF [метры].
 * @param[in]  lat_gcen  Геоцентрическая широта станции [радианы].
 * @param[in]  lon_gcen  Геоцентрическая долгота станции [радианы].
 * @param[in]  sun       Матрица 3x2: вектор положения и скорости Солнца в J2000 [м, м/с].
 * @param[in]  moon      Матрица 3x2: вектор положения и скорости Луны в J2000 [м, м/с].
 * @param[in]  f         Вектор (5) фундаментальных аргументов (IERS) [арксекунды].
 * @param[in]  fd        Вектор (5) производных фундаментальных аргументов [арксек/столетие].
 * @param[in]  vw_i      Матрица перехода 3x3 из топоцентрической системы (VEN) в ITRF.
 * @param[in]  gast      Вектор (2): Истинное звездное время по Гринвичу (GAST) и его производная [рад, рад/с].
 * @param[in]  r2000     Матрица 9x3 (или блок матриц): переходы r0, r1, r2 из ITRF в J2000.
 * @param[out] dxtide    Вектор смещения станции в небесной системе J2000 [метры].
 * @param[out] dvtide    Вектор скорости смещения станции в небесной системе J2000 [м/с].
 */
void SITE_TIDE_SOLID(const Eigen::Vector3d& xsta_itrf, double lat_gcen, double lon_gcen,
                     const Eigen::Matrix<double, 3, 2>& sun, const Eigen::Matrix<double, 3, 2>& moon,
                     const Eigen::VectorXd& f, const Eigen::VectorXd& fd,
                     const Eigen::Matrix3d& vw_i, const Eigen::Vector2d& gast,
                     const Eigen::MatrixXd& r2000,
                     Eigen::Vector3d& dxtide, Eigen::Vector3d& dvtide);

/**
 * @brief Вычисляет смещение и скорость станции, вызванные океаническими приливами (Ocean Tide Loading).
 */
void SITE_TIDE_OC(int mjd_start, double ut1_sec, const ariadna::OceanTideData& tide_data,
                  const Eigen::Matrix3d& vw_i, const Eigen::MatrixXd& r2000,
                  Eigen::Vector3d& dx_octide, Eigen::Vector3d& dv_octide);

/**
 * @brief Расчет координат станции на эпоху наблюдения с учетом движения тектонических плит.
 * Вычисляет геодезические координаты (GEOD) и формирует матрицу перехода 
 * из топоцентрической системы (VEN) в геоцентрическую (ITRF).
 * Для космических телескопов тектоника и геодезия обнуляются.
 *
 * @param[in]  station  Объект станции (содержит базовые xyz и скорости).
 * @param[in]  dYear    Разница в годах между базовой эпохой каталога и эпохой наблюдения.
 * @param[out] out      Структура SiteCorrData с обновленными параметрами и матрицей vw_i.
 */
void SITE_CORR(const Station& station, double dYear, SiteCorrData& out);

/**
 * @brief Преобразование декартовых геоцентрических координат в геодезические (Borkowski, 1989).
 * Функция реализует точное (не итерационное) решение для пересчета координат точки 
 * относительно земного эллипсоида.
 * * @param[in]  equatorial_radius_r  Расстояние от оси вращения Земли (sqrt(X^2 + Y^2)) [метры].
 * @param[in]  z_polar              Координата по оси Z (расстояние от экваториальной плоскости) [метры].
 * @param[out] geodetic_latitude_fi Рассчитанная геодезическая широта точки [радианы].
 * @param[out] geodetic_height_h    Рассчитанная геодезическая высота над поверхностью эллипсоида [метры].
 */
void GEOD(double equatorial_radius_r, double z_polar, double& geodetic_latitude_fi, double& geodetic_height_h);

/**
 * @brief Вспомогательная функция для расчета аргументов океанических приливов.
 */
void calc_tide_angles(int mjd_start, double ut1_sec, Eigen::VectorXd& angle, Eigen::VectorXd& speed_angle);

/**
 * @brief Вычисляет смещение и скорость станции, вызванные атмосферной нагрузкой (Atmospheric Loading).
 */
void SITE_ATM40(double dPdt, const Eigen::Matrix3d& vw_i, const Eigen::MatrixXd& r2000,
                Eigen::Vector3d& dx_atm, Eigen::Vector3d& dv_atm);

/**
 * @brief Вычисляет влияние атмосферной нагрузки на координаты двух станций интерферометра.
 */
void site_atm40(int j1, int j2, const std::vector<Station>& stations, const Observation& observation,
                const Eigen::Vector2d& dPdt, const std::vector<Eigen::Matrix3d>& vw,
                const Eigen::MatrixXd& r2000, Eigen::Matrix<double, 3, 2>& dx_atm, Eigen::Matrix<double, 3, 2>& dv_atm);

/**
 * @brief Вычисляет итоговую позицию и скорость станции в системе J2000.0 
 * путем суммирования всех поправок (приливы, полюс, атмосфера, температура).
 */
void SITE_INST(const Eigen::Vector3d& xsta_itrf, const Eigen::Vector3d& vsta_itrf,
               const Eigen::MatrixXd& r2000,
               const Eigen::Vector3d& dx_tide, const Eigen::Vector3d& dv_tide,
               const Eigen::Vector3d& dx_octide, const Eigen::Vector3d& dv_octide,
               const Eigen::Vector3d& dx_poltide, const Eigen::Vector3d& dv_poltide,
               const Eigen::Vector3d& dx_atm, const Eigen::Vector3d& dv_atm,
               const Eigen::Vector3d& dx_temp, const Eigen::Vector3d& dv_temp,
               Eigen::Vector3d& xsta_j2000t, Eigen::Vector3d& vsta_j2000t,
               Eigen::Vector3d& asta_j2000);

/**
 * @brief (Перегрузка) Вычисляет итоговую позицию и скорость станции в системе J2000.0 (внутренний вызов).
 */
void site_inst(const Eigen::Vector3d& xsta_itrf, const Eigen::Vector3d& vsta_itrf,
               const Eigen::MatrixXd& r2000,
               const Eigen::Vector3d& dx_tide, const Eigen::Vector3d& dv_tide,
               const Eigen::Vector3d& dx_octide, const Eigen::Vector3d& dv_octide,
               const Eigen::Vector3d& dx_poltide, const Eigen::Vector3d& dv_poltide,
               const Eigen::Vector3d& dx_atm, const Eigen::Vector3d& dv_atm,
               const Eigen::Vector3d& dx_temp, const Eigen::Vector3d& dv_temp,
               Eigen::Vector3d& xsta_j2000t, Eigen::Vector3d& vsta_j2000t,
               Eigen::Vector3d& asta_j2000);

// ============================================================================
// Астрометрия и матрицы переходов (Эфемериды, прецессия, нутация)
// ============================================================================

/**
 * @brief Вычисляет единичные векторы направлений на радиоисточники.
 */
void source_vec(const std::vector<Source>& sources, double t_mean, std::vector<Eigen::Vector3d>& k_star);

/**
 * @brief Получение координат Земли, Солнца и Луны на основе эфемерид JPL (матричный формат).
 */
void jpl_eph(double jd, double ct, Eigen::Matrix3d& earth, Eigen::MatrixXd& sun, Eigen::MatrixXd& moon);

/**
 * @brief Получение координат Земли, Солнца и Луны на основе эфемерид JPL (векторный формат).
 */
void jpleph(double jd, double ct, Eigen::Vector3d& earth, Eigen::Vector3d& sun, Eigen::Vector3d& moon);

/**
 * @brief Вычисление наклона эклиптики к экватору по модели IAU 2006.
 */
void eps_a06(double jd, double ct, double& eps2000, double& eps_p03_2000, Eigen::VectorXd& e_mn);

/**
 * @brief Формирование матрицы прецессии и ее производных.
 */
void prec_matrix(double jd, double ct, double eps_p03_2000, Eigen::MatrixXd& pr, Eigen::MatrixXd& dpdp_ls, Eigen::MatrixXd& dpdp_pl, Eigen::MatrixXd& dpdp_om);

/**
 * @brief Формирование матрицы смещения (bias) для перехода от GCRS к J2000.
 */
void bias(double eps2000, Eigen::Matrix3d& bias_matr);

/**
 * @brief Вычисление углов нутации и связанных матриц по моделям IAU 2000/2006.
 */
void iau_2000_2006(double jd, double ct, const Eigen::VectorXd& f, const Eigen::VectorXd& fd, double eop_int_x, double deop_int_x, double eop_int_y, double deop_int_y, double eps_p03_2000, const Eigen::VectorXd& e_mn, Eigen::VectorXd& dpsir, Eigen::VectorXd& depsr, Eigen::VectorXd& eps, Eigen::MatrixXd& rn, Eigen::MatrixXd& dndpsi, Eigen::MatrixXd& dndeps);

/**
 * @brief Расчет Гринвичского истинного звездного времени (GAST) и матрицы вращения Земли (Earth Rotation).
 */
void gas_time(double jd, double ut1, double ct, const Eigen::VectorXd& f, const Eigen::VectorXd& fd, const Eigen::VectorXd& dpsir, const Eigen::VectorXd& depsr, const Eigen::VectorXd& e_mn, const Eigen::VectorXd& eps, double dtaidct, double deop_int_ut1, double& diurnv, Eigen::VectorXd& gast, Eigen::MatrixXd& rs, Eigen::MatrixXd& drsdp_ls, Eigen::MatrixXd& drsdp_pl);

/**
 * @brief Формирование матрицы движения полюса (Wobble) и ее производных.
 */
void wobble(double cent, double eop_x, double eop_y, double deop_x, double deop_y, Eigen::MatrixXd& ryx, Eigen::Matrix3d& ydxdx, Eigen::Matrix3d& dydyx, Eigen::Matrix3d& ddxdyx, Eigen::Matrix3d& ddydyx);

/**
 * @brief Построение полной матрицы перехода от ITRF к J2000 (учет вращения, движения полюса, прецессии и нутации).
 */
void r2000_matrix(double mjd, double ut1, const Eigen::VectorXd& eop_int, const Eigen::VectorXd& deop_int, int i_choice, Eigen::Matrix3d r2000[3], double& gast);


// ============================================================================
// Расчет задержек, производных и тропосферы
// ============================================================================

/**
 * @brief Вычисление метеорологических градиентов/производных во времени (для всех станций).
 */
void dmeteo1_dt(const std::vector<Observation>& observations, const std::vector<Station>& stations, double t_mean, std::vector<Eigen::Vector3d>& site_meteo, int& ndeg, std::vector<Eigen::VectorXd>& t_coef, std::vector<Eigen::VectorXd>& p_coef, std::vector<Eigen::VectorXd>& hum_coef);

/**
 * @brief Вычисление метеорологических градиентов/производных во времени для конкретной станции.
 */
void dmeteo2_dt(const Station& station, int n_stations, int ndeg, const Observation& obs, double t_mean, const std::vector<Eigen::VectorXd>& t_coef, const std::vector<Eigen::VectorXd>& p_coef, const std::vector<Eigen::VectorXd>& hum_coef, double& dtdt, double& dpdt, double& dhumdt);

/**
 * @brief Вычисляет положение источника, скорректированное за годовую и суточную аберрацию.
 */
void aber_source(const Observation& obs, const std::vector<Eigen::Matrix3d>& r2000, const Eigen::Vector3d& k_s, const Eigen::Matrix<double, 3, 3>& earth, const std::vector<Eigen::Vector3d>& vsta_j2000t, const std::vector<Eigen::Matrix3d>& vw, Eigen::Matrix2d& e, Eigen::Matrix2d& az);

/**
 * @brief Вычисляет вектор базы в инерциальной (J2000) и земной (Crust-fixed) системах координат.
 */
void baseline(const Eigen::Matrix3d& r2000, const Eigen::MatrixXd& xsta_j2000t, const Eigen::MatrixXd& vsta_j2000t, const Eigen::MatrixXd& asta_j2000, Eigen::MatrixXd& base_line, Eigen::Vector3d& b_cfs);

/**
 * @brief Проекция вектора базы интерферометра на UV-плоскость.
 */
void uv_plane(const Source& source, const std::vector<Eigen::Vector3d>& base_line, const std::vector<Eigen::Vector3d>& xsta_j2000t, double scale, Eigen::Vector3d& uv_coor);


/**
 * @brief Расчет теоретической задержки для баз с участием космических аппаратов (Space VLBI).
 */
void theor_delay_orb(const std::vector<Eigen::Vector3d>& base_line, const std::vector<Eigen::Vector3d>& xsta_j2000t, const std::vector<Eigen::Vector3d>& vsta_j2000t, const std::vector<Eigen::Vector3d>& asta_j2000, const Eigen::Vector3d& k_s, const Eigen::Vector3d& earth, const Eigen::MatrixXd& sun, const Eigen::MatrixXd& moon, const Eigen::MatrixXd& datmc_d, const Eigen::MatrixXd& datmc_w, const Eigen::MatrixXd& dtau_off, double& t2_t1, double& dt2_t1);


// ============================================================================
// Частные производные для МНК (Least Squares)
// ============================================================================

/**
 * @brief Вычисление частных производных задержки по координатам радиоисточника.
 */
void der_star(double jd, double ct, double dyear, const Eigen::Vector3d& k_s, const Eigen::Matrix3d& earth, const std::vector<Eigen::Vector3d>& base_line, Eigen::MatrixXd& dstar, Eigen::MatrixXd& dstar_rate);

/**
 * @brief Вычисление частных производных задержки по координатам станций.
 */
void der_site(double dyear, const Eigen::MatrixXd& r2000, const Eigen::Vector3d& k_s, const Eigen::Matrix3d& earth, const std::vector<Eigen::Vector3d>& base_line, std::vector<Eigen::MatrixXd>& dsite, std::vector<Eigen::MatrixXd>& dsite_v);

/**
 * @brief Вычисление частных производных задержки по координатам полюса (x_p, y_p).
 */
void der_polar(const Eigen::MatrixXd& r2000, const Eigen::Vector3d& k_s, const Eigen::Matrix3d& earth, const std::vector<Eigen::Vector3d>& base_line, const Eigen::Vector3d& b_cfs, std::vector<Eigen::MatrixXd>& pr, std::vector<Eigen::MatrixXd>& rn, std::vector<Eigen::MatrixXd>& rs, std::vector<Eigen::MatrixXd>& ryx, std::vector<Eigen::MatrixXd>& ydxdx, std::vector<Eigen::MatrixXd>& dydyx, std::vector<Eigen::MatrixXd>& ddxdyx, std::vector<Eigen::MatrixXd>& ddydyx, std::vector<Eigen::VectorXd>& dx_pol_dx, std::vector<Eigen::VectorXd>& dx_pol_dy, const Eigen::VectorXd& arg_oc_tide, Eigen::MatrixXd& dwob, std::vector<Eigen::MatrixXd>& dx_aj, std::vector<Eigen::MatrixXd>& dx_bj, std::vector<Eigen::MatrixXd>& dy_aj, std::vector<Eigen::MatrixXd>& dy_bj);

/**
 * @brief Вычисление частных производных задержки по параметру вращения Земли (UT1-TAI).
 */
void der_ut1(const Eigen::MatrixXd& r2000, const Eigen::Vector3d& k_s, const Eigen::Matrix3d& earth, const std::vector<Eigen::Vector3d>& base_line, const Eigen::Vector3d& b_cfs, const Eigen::VectorXd& gast, std::vector<Eigen::MatrixXd>& pr, std::vector<Eigen::MatrixXd>& rn, std::vector<Eigen::MatrixXd>& rs, std::vector<Eigen::MatrixXd>& ryx, const Eigen::VectorXd& arg_oc_tide, double diurnv, Eigen::VectorXd& dut1_tai, std::vector<Eigen::MatrixXd>& dut1_aj, std::vector<Eigen::MatrixXd>& dut1_bj);

/**
 * @brief Вычисление частных производных задержки по углам нутации.
 */
void der_nut(const Eigen::MatrixXd& r2000, const Eigen::Vector3d& k_s, const Eigen::Matrix3d& earth, const std::vector<Eigen::Vector3d>& base_line, const Eigen::Vector3d& b_cfs, const Eigen::VectorXd& gast, std::vector<Eigen::MatrixXd>& pr, std::vector<Eigen::MatrixXd>& rn, std::vector<Eigen::MatrixXd>& rs, std::vector<Eigen::MatrixXd>& ryx, const Eigen::VectorXd& e_mn, const Eigen::MatrixXd& dndpsi, const Eigen::MatrixXd& dndeps, Eigen::MatrixXd& dnut);

/**
 * @brief Вычисление частных производных задержки по параметрам прецессии.
 */
void der_prec(const Eigen::MatrixXd& r2000, const Eigen::Vector3d& k_s, const Eigen::Matrix3d& earth, const std::vector<Eigen::Vector3d>& base_line, const Eigen::Vector3d& b_cfs, const Eigen::MatrixXd& pr, const Eigen::MatrixXd& rn, const Eigen::MatrixXd& rs, const Eigen::MatrixXd& ryx, const Eigen::MatrixXd& dpdp_ls, const Eigen::MatrixXd& dpdp_pl, const Eigen::MatrixXd& drsdp_ls, const Eigen::MatrixXd& drsdp_pl, Eigen::MatrixXd& dpr_lspl);

/**
 * @brief Вычисление частных производных задержки по числам Лява (h, l) твердых земных приливов.
 */
void der_love_number(const Station& sta1, const Station& sta2, const Eigen::MatrixXd& r2000, const Eigen::Vector3d& k_s, const Eigen::Matrix3d& earth, const std::vector<Eigen::Vector3d>& base_line, const std::vector<Eigen::Matrix3d>& drdh0_3, const std::vector<Eigen::Matrix3d>& drdh02_3, const std::vector<Eigen::Matrix3d>& drdrl0_3, const std::vector<Eigen::Matrix3d>& drdl02_3, const std::vector<Eigen::Matrix3d>& drdh3_3, const std::vector<Eigen::Matrix3d>& drdl3_3, const std::vector<Eigen::Matrix3d>& drdl1_1_2000_3, const std::vector<Eigen::Matrix3d>& drdl1_2_2000_3, const std::vector<Eigen::Matrix3d>& drdhi_1_2000_3, const std::vector<Eigen::Matrix3d>& drdhi_2_2000_3, const std::vector<Eigen::Matrix3d>& drdli_1_2000_3, const std::vector<Eigen::Matrix3d>& drdli_2_2000_3, Eigen::VectorXd& d_dh0, Eigen::VectorXd& d_dh02, Eigen::VectorXd& d_dl0, Eigen::VectorXd& d_dl02, Eigen::VectorXd& d_dh3, Eigen::VectorXd& d_dl3, Eigen::VectorXd& d_dl1_1, Eigen::VectorXd& d_dl1_2, Eigen::VectorXd& d_dhi_1, Eigen::VectorXd& d_dhi_2, Eigen::VectorXd& d_dli_1, Eigen::VectorXd& d_dli_2);

/**
 * @brief Вычисление частных производных задержки по параметрам атмосферной нагрузки.
 */
void der_atm_load(const Station& sta1, const Station& sta2, const Eigen::MatrixXd& r2000, const Eigen::Vector3d& k_s, const Eigen::Matrix3d& earth, const std::vector<Eigen::Vector3d>& base_line, const std::vector<Eigen::Matrix3d>& dr1_da2000, const std::vector<Eigen::Matrix3d>& dr1_db2000, const std::vector<Eigen::Matrix3d>& dv1_da2000, const std::vector<Eigen::Matrix3d>& dv1_db2000, const std::vector<Eigen::Matrix3d>& dr2_da2000, const std::vector<Eigen::Matrix3d>& dr2_db2000, const std::vector<Eigen::Matrix3d>& dv2_da2000, const std::vector<Eigen::Matrix3d>& dv2_db2000, const std::vector<Eigen::Vector3d>& dr_dreg2000, const std::vector<Eigen::Vector3d>& dv_dreg2000, std::vector<Eigen::MatrixXd>& dt_da, std::vector<Eigen::MatrixXd>& dt_db, std::vector<Eigen::MatrixXd>& df_da, std::vector<Eigen::MatrixXd>& df_db, Eigen::VectorXd& dt_dreg, Eigen::VectorXd& df_dreg);

/**
 * @brief Формирование матрицы условных уравнений (design matrix) и вектора невязок для метода наименьших квадратов (МНК).
 */
void create_matr(const Observation& obs, const Station& sta1, const Station& sta2, int i_good, int nobs, int n_sources, int n_stations, double t_mean, double t2_t1, double dt2_t1, const Eigen::MatrixXd& dwob, const Eigen::VectorXd& dut1_tai, const Eigen::MatrixXd& dnut, const std::vector<Eigen::MatrixXd>& dsite, const std::vector<Eigen::MatrixXd>& dsite_v, const Eigen::MatrixXd& datmp_hmf, const Eigen::MatrixXd& datmp_wmf, const Eigen::MatrixXd& dgrad_n, const Eigen::MatrixXd& dgrad_e, const Eigen::MatrixXd& dstar, const Eigen::MatrixXd& dstar_rate, const Eigen::MatrixXd& dpr_lspl, const Eigen::VectorXd& d_dh0, const Eigen::VectorXd& d_dh02, const Eigen::VectorXd& d_dl0, const Eigen::VectorXd& d_dl02, const Eigen::VectorXd& d_dh3, const Eigen::VectorXd& d_dl3, const Eigen::VectorXd& d_dl1_1, const Eigen::VectorXd& d_dl1_2, const Eigen::VectorXd& d_dhi_1, const Eigen::VectorXd& d_dhi_2, const Eigen::VectorXd& d_dli_1, const Eigen::VectorXd& d_dli_2, const Eigen::MatrixXd& d_dax, const std::vector<Eigen::MatrixXd>& dt_da, const std::vector<Eigen::MatrixXd>& dt_db, const std::vector<Eigen::MatrixXd>& df_da, const std::vector<Eigen::MatrixXd>& df_db, const Eigen::VectorXd& dt_dreg, const Eigen::VectorXd& df_dreg, const std::vector<Eigen::MatrixXd>& dx_aj, const std::vector<Eigen::MatrixXd>& dx_bj, const std::vector<Eigen::MatrixXd>& dy_aj, const std::vector<Eigen::MatrixXd>& dy_bj, const std::vector<Eigen::MatrixXd>& dut1_aj, const std::vector<Eigen::MatrixXd>& dut1_bj, const Eigen::MatrixXd& zen_dry, const Eigen::MatrixXd& zen_wet, Eigen::MatrixXd& m_pd, Eigen::VectorXd& y, Eigen::VectorXd& w);

/**
 * @brief Численное интегрирование орбиты космического аппарата и вычисление фазы для Space VLBI.
 */
void integr8_asc(const Station& station, const Observation& obs, int nobs, int l_segm, int n_wr_tot, double delta_sec, const std::string& track_site, int mjd_beg, double utc_beg, double ct_beg, double dyear, const Eigen::MatrixXd& r2000, std::vector<Eigen::Vector3d>& xsta_j2000t, std::vector<Eigen::Vector3d>& vsta_j2000t, std::vector<Eigen::Vector3d>& asta_j2000, Eigen::MatrixXd& r_sat_pr, const Eigen::Vector3d& k_s, double& phase1);

/**
 * @brief Главная функция обработки наблюдений: расчет геометрических, тропосферных и ионосферных задержек, производных и формирование уравнений для МНК.
 */
void process_ariadna(const std::vector<Station>& stations, const std::vector<Source>& sources, const std::vector<Observation>& observations, const std::vector<SpaceStation>& space_stations, const std::vector<OrbitData>& orbit_data, int n_segm, int k_ch_c, int k_ch_z, double delta_sec, const std::string& output_path, const std::vector<EOPData>& eop_data, double mjd_mean, double utc_mean, double t_mean);

} // namespace ariadna