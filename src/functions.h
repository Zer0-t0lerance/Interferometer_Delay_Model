#pragma once
#include "structures.h"
#include "READ_CAT.h"

namespace ariadna {

// ============================================================================
// Конвертация времени и параметров ориентации Земли (EOP)
// ============================================================================

/**
 * @brief Конвертация времени из шкал MJD/UTC в атомное (TAI) и земное (TT).
 *
 * TAI = UTC + високосные секунды (nsec); TT = TAI + 32.184 с. Всё — доли суток.
 *
 * @param[in]  mjd  Модифицированная юлианская дата на 0h (целые сутки).
 * @param[in]  UTC  Доля суток UTC.
 * @param[out] TAI  Международное атомное время (доля суток).
 * @param[out] TT   Земное время (доля суток).
 */
void tai_time(double mjd, double UTC, double &TAI, double &TT);

/**
 * @brief Число високосных секунд (TAI-UTC) на заданную дату.
 * @param[in]  mjd    Модифицированная юлианская дата.
 * @param[out] idelt  Число високосных секунд TAI-UTC [с].
 */
void nsec(double mjd, double& idelt);

/**
 * @brief Интерполяция параметров ориентации Земли по стандарту IERS (модель Bizouard).
 *
 * 4-точечный Лагранж по 7 узлам + суточные/полусуточные поправки (океан + либрация).
 * ВНИМАНИЕ по единицам: выход НЕ в радианах — eop_int(0)=UT1-UTC [сек],
 * eop_int(1..2)=x,y [угл.сек], eop_int(3..4)=dpsi,deps [угл.сек].
 *
 * @param[in]  obs      Наблюдение (MJD и UTC).
 * @param[in]  eop_data 7 опорных точек EOP вокруг момента наблюдения.
 * @param[out] eop_int  Вектор из 5 значений [UT1-UTC, x, y, dpsi, deps].
 */
void interp_iers(const Observation& obs, const std::vector<EOPData>& eop_data,
                 Eigen::VectorXd& eop_int);

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
 * @brief Подгонка метеопараметров станций полиномом по времени (порт DMETEO1_DT).
 *
 * Для каждой станции собирает точки (температура, давление, влажность) из всех её
 * наблюдений и строит МНК-полином степени ndeg в аргументе (t - t_mean) [сутки].
 * Коэффициенты идут в dmeteo2_dt для вычисления производных.
 *
 * @param[in]  observations Наблюдения сессии (метео t1/p1/e1, t2/p2/e2 по станциям).
 * @param[in]  n_stations   Число станций.
 * @param[in]  t_mean       Средняя эпоха сессии [MJD] (центр аргумента полинома).
 * @param[in]  ndeg         Степень полинома (обычно 2); при нехватке точек снижается.
 * @param[out] t_coef       Коэффициенты полинома температуры по станциям (ndeg+1, по возрастанию).
 * @param[out] p_coef       Коэффициенты полинома давления по станциям.
 * @param[out] hum_coef     Коэффициенты полинома влажности по станциям.
 */
void dmeteo1_dt(const std::vector<Observation>& observations, int n_stations,
                double t_mean, int ndeg,
                std::vector<Eigen::VectorXd>& t_coef,
                std::vector<Eigen::VectorXd>& p_coef,
                std::vector<Eigen::VectorXd>& hum_coef);

/**
 * @brief Производные метеопараметров станции на момент наблюдения (порт DMETEO2_DT).
 *
 * Возвращает dT/dt, dP/dt, dHum/dt как производные полиномов (dmeteo1_dt) в точке
 * (mjd+utc). Единицы — на СУТКИ: dPdt [мбар/сут] -> SITE_ATM40, dTdt [°C/сут] -> THERM_DEF40.
 *
 * @param[in]  ista            Индекс станции.
 * @param[in]  ndeg            Степень полинома (не используется явно; для совместимости).
 * @param[in]  mjd, utc        Эпоха наблюдения.
 * @param[in]  t_mean          Средняя эпоха (центр полинома).
 * @param[in]  t_coef,p_coef,hum_coef  Коэффициенты из dmeteo1_dt.
 * @param[out] dTdt,dPdt,dHumdt Производные [ед./сут].
 */
void dmeteo2_dt(int ista, int ndeg, int mjd, double utc, double t_mean,
                const std::vector<Eigen::VectorXd>& t_coef,
                const std::vector<Eigen::VectorXd>& p_coef,
                const std::vector<Eigen::VectorXd>& hum_coef,
                double& dTdt, double& dPdt, double& dHumdt);

/**
 * @brief Суточные/полусуточные приливные поправки к EOP (71 океаническая гармоника).
 * @param[in]  cent         Юлианские столетия от J2000.0 (T).
 * @param[in]  f            5 фундаментальных аргументов (l, l', F, D, Omega) [рад].
 * @param[in]  fd           Производные аргументов [рад/сут].
 * @param[out] dEOP_diu     Поправки полюса/UT1: 3x2 (x, y, UT1 × значение/скорость).
 * @param[out] arg_oc_tide  Аргументы океанических приливов (для дальнейших расчётов).
 */
void terms_71(double cent, const Eigen::VectorXd& f, const Eigen::VectorXd& fd, Eigen::MatrixXd& dEOP_diu, Eigen::VectorXd& arg_oc_tide);

/**
 * @brief Поправки к EOP от либрации нежёсткой Земли (суточные члены полюса/UT1).
 * @param[in]  cent      Юлианские столетия от J2000.0 (T).
 * @param[in]  f         5 фундаментальных аргументов [рад].
 * @param[in]  fd        Производные аргументов [рад/сут].
 * @param[out] dEOP_lib  Поправки полюса/UT1 от либрации: 3x2 (x, y, UT1 × значение/скорость).
 */
void terms_lib(double cent, const Eigen::VectorXd& f, const Eigen::VectorXd& fd, Eigen::MatrixXd& dEOP_lib);


// ============================================================================
// Вычисление координат, систем отсчета и геофизических эффектов
// ============================================================================

/**
 * @brief Инициализация базовых координат станций (аналог SUBROUTINE SITE в Фортране).
 *
 * Для каждой станции переносит ITRF-координаты на эпоху наблюдения (с учётом
 * скорости плит), вычисляет геоцентрические (широта/долгота/радиус) и
 * геодезические (широта/высота через GEOD) координаты и матрицу перехода
 * VEN -> ITRF. Космический телескоп ("RASTRON") берётся из space_stations,
 * геодезия для него обнуляется.
 *
 * @param[in]  stations       Список станций (базовые ITRF-координаты, скорости, имена).
 * @param[in]  space_stations Данные космического телескопа (позиция/скорость в GCRS).
 * @param[in]  n_stations     Число станций.
 * @param[in]  dyear          Годы от эпохи каталога до эпохи наблюдения.
 * @param[out] site_xyz       Геоцентрические ITRF-координаты станций на эпоху [м].
 * @param[out] site_vel       Скорости станций [м/год].
 * @param[out] lat_geod       Геодезическая широта [рад].
 * @param[out] h_geod         Геодезическая высота над эллипсоидом [м].
 * @param[out] lat_gcen       Геоцентрическая широта [рад].
 * @param[out] lon_gcen       Геоцентрическая долгота [рад].
 * @param[out] sph_rad        Сферический радиус (|xyz|) [м].
 * @param[out] u_site         Расстояние от оси вращения (экваториальный радиус) [км].
 * @param[out] v_site         Расстояние от экваториальной плоскости (Z) [км].
 * @param[out] vw             Матрицы перехода VEN (Vertical, East, North) -> ITRF (3x3).
 */
void site(const std::vector<Station>& stations, const std::vector<SpaceStation>& space_stations, int n_stations, double dyear,
          std::vector<Eigen::Vector3d>& site_xyz, std::vector<Eigen::Vector3d>& site_vel, std::vector<double>& lat_geod,
          std::vector<double>& h_geod, std::vector<double>& lat_gcen,
          std::vector<double>& lon_gcen, std::vector<double>& sph_rad,
          std::vector<double>& u_site, std::vector<double>& v_site,
          std::vector<Eigen::Matrix3d>& vw);

// ПРИМЕЧАНИЕ: перегрузка site(...) для пары станций (посуточная сборка координат
// в J2000 с приливными/нагрузочными поправками) будет восстановлена вместе со
// слоем оркестрации (см. project-ariadna-port-status). Её черновик был удалён
// как несобираемый.

/**
 * @brief Вычисление полной теоретической задержки и её производной (порт THEOR_DELAY).
 *
 * Соответствует операционному прогону ARIADNA (подтверждён отладочным дампом):
 * гравитационная задержка Шапиро (Солнце, Луна, Земля) + геометрия с поправками
 * на движение станций (включая члены term2e/term2f для Радиоастрона) + тропосфера
 * (Datmc) + инструментальная задержка оси монтировки (dtau_off) + термодеформация (dt_temp).
 *
 * Соглашение о матрицах 2x2 (Datmc_d/Datmc_w/dtau_off/dt_temp): СТРОКА = задержка(0)/скорость(1),
 * СТОЛБЕЦ = станция 1(0)/2(1). Т.е. элемент (0,j) — задержка станции j, (1,j) — её производная.
 *
 * @param[in]  base_line   Матрица 3x2: вектор базы (от ст.1 к ст.2) и его производная [м, м/с].
 * @param[in]  xsta        Векторы (2): геоцентрические координаты станций 1 и 2 в J2000 [м].
 * @param[in]  vsta        Векторы (2): геоцентрические скорости станций 1 и 2 в J2000 [м/с].
 * @param[in]  asta        Векторы (2): геоцентрические ускорения станций 1 и 2 в J2000 [м/с^2].
 * @param[in]  K_s         Единичный вектор направления на источник в J2000.
 * @param[in]  Earth       Матрица 3x3: координаты, скорость и ускорение барицентра Земли в SSB [м, м/с, м/с^2].
 * @param[in]  Sun         Матрица 3x3: координаты, скорость и ускорение Солнца в SSB [м, м/с, м/с^2].
 * @param[in]  Moon        Матрица 3x3: координаты, скорость и ускорение Луны в SSB [м, м/с, м/с^2].
 * @param[in]  Datmc_d     Матрица 2x2: сухая тропосфера (строка - задержка/скорость, столбец - станция) [с, с/с].
 * @param[in]  Datmc_w     Матрица 2x2: влажная тропосфера [с, с/с].
 * @param[in]  dtau_off    Матрица 2x2: инструментальная задержка оси монтировки [с, с/с].
 * @param[in]  dt_temp     Матрица 2x2: температурная поправка антенны [с, с/с].
 * @param[out] t2_t1       Полная теоретическая задержка [секунды].
 * @param[out] dt2_t1      Производная задержки по времени [секунды/секунду].
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
 * @brief Вычисляет смещение и скорость станции, вызванные приливом полюса (Pole Tide).
 * Расчет производится в соответствии со стандартами IERS 2000 с учетом векового дрейфа 
 * среднего полюса. Возвращает векторы смещения и скорости в небесной системе J2000.0.
 *
 * @param[in]  cent        Количество юлианских столетий от эпохи J2000.0 (для расчета векового дрейфа).
 * @param[in]  lat_geod    Геодезическая широта станции (радианы).
 * @param[in]  lon_geod    Геодезическая долгота станции (радианы).
 * @param[in]  xp          Координата X полюса Земли из EOP (секунды дуги).
 * @param[in]  yp          Координата Y полюса Земли из EOP (секунды дуги).
 * @param[in]  xp_rate     Скорость изменения координаты X полюса (секунды дуги в секунду).
 * @param[in]  yp_rate     Скорость изменения координаты Y полюса (секунды дуги в секунду).
 * @param[in]  vw_i        Матрица перехода от локальной топоцентрической системы (VEN) к земной ITRF (3x3).
 * @param[in]  r2000       Матрица перехода от земной системы ITRF к небесной J2000.0 (3x3).
 * @param[in]  dr2000_dt   Первая производная матрицы поворота r2000 по времени (3x3).
 * @param[out] dx_poltide  Выходной вектор смещения коры в системе J2000.0 (метры).
 * @param[out] dv_poltide  Выходной вектор скорости коры в системе J2000.0 (метры в секунду).
 */
void POLE_TIDE(double cent, double lat_geod, double lon_geod,
               double xp, double yp, double xp_rate, double yp_rate,
               const Eigen::Matrix3d& vw_i, const Eigen::Matrix3d& r2000,
               const Eigen::Matrix3d& dr2000_dt,
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
 * @brief Переносит коэффициенты океанической нагрузки из сырых данных каталога в структуры станций.
 * * Функция проходит по списку станций, ищет совпадение имени в сыром массиве данных каталога 
 * и перекладывает плоские массивы коэффициентов (33 элемента) в двумерные матрицы Eigen (3x11) 
 * (оси: Up, West, South; 11 приливных волн). Если данных для станции нет, заполняет матрицы нулями.
 *
 * @param[in]     raw_oc_data  Вектор "сырых" записей каталога океанических приливов (прочитанных из ocload40.dat).
 * @param[in,out] stations     Вектор структур станций. При совпадении имен поле `tide_data` будет заполнено.
 */
void map_ocean_tides_to_stations(const std::vector<oc_record>& raw_oc_data, std::vector<Station>& stations);

/**
 * @brief Переносит коэффициенты атмосферной нагрузки из сырого каталога в структуры станций.
 *
 * Ищет совпадение имени станции в массиве atm_record (каталог VLBI_atmload4_12.cat)
 * и заполняет Station.atm_load.coef (3 компоненты VEN × 8 коэффициентов A1 B1 A2 B2 A3 B3 b0 b1,
 * амплитуды в мм). Опорное давление p_0 заполняется отдельно из каталога ANTENNA_INFO.
 * Если станция не найдена — coef обнуляется, has_data = false.
 *
 * @param[in]     raw_atm_data Вектор записей каталога атмосферной нагрузки (ReadATM).
 * @param[in,out] stations     Станции; при совпадении имён заполняется поле atm_load.
 */
void map_atm_loading_to_stations(const std::vector<atm_record>& raw_atm_data, std::vector<Station>& stations);

/**
 * @brief Вычисляет смещение и скорость станции, вызванные океаническими приливами (Ocean Tide Loading).
 * * Вычисляет эффекты от 11 основных океанических приливных волн по модели Швидерского 
 * (Внимание: фазы Швидерского сдвинуты на +/- 90 градусов и не совпадают со стандартами Дудсона/Картрайта).
 * Переводит топоцентрические смещения в геоцентрическую систему (ITRF), а затем в J2000.
 *
 * @param[in]  jd          Юлианская дата на 0 часов UTC текущих суток (дни).
 * @param[in]  ut1_sec     Интерполированное значение времени UT1 в течение суток (секунды).
 * @param[in]  tide_data   Структура амплитуд (метры) и фаз (градусы) 11 волн по 3 осям: Up, West, South.
 * @param[in]  vw_i        Матрица перехода от локальной системы (VEN) к земной (ITRF).
 * @param[in]  r2000       Матрица перехода от ITRF к J2000 (в Фортране r2000(1,1,1)).
 * @param[in]  dr2000_dt   Первая производная матрицы перехода r2000 по времени (в Фортране r2000(1,1,2)).
 * @param[out] dx_octide   Выходной вектор смещения коры в системе J2000 (метры).
 * @param[out] dv_octide   Выходной вектор скорости коры в системе J2000 (метры в секунду).
 */
void SITE_TIDE_OC(double jd, double ut1_sec, const OceanTideData& tide_data,
                  const Eigen::Matrix3d& vw_i, 
                  const Eigen::Matrix3d& r2000, 
                  const Eigen::Matrix3d& dr2000_dt,
                  Eigen::Vector3d& dx_octide, Eigen::Vector3d& dv_octide);

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
 * @brief Смещение и скорость станции из-за атмосферной нагрузки в системе J2000.0 (порт SITE_ATM40).
 *
 * Топоцентрическое смещение по модели каталога VLBI_atmload4_12.cat:
 *   dr(c) = Sum_{i=1..3}[A_i*cos(omega_i*t) + B_i*sin(omega_i*t)]  для всех компонент c=V,E,N,
 *   а для вертикали (V) дополнительно + b0 + b1*(P - p_0),
 * где t = MJD - ATM_EPOCH_MJD [сут], omega_i = TWOPI/ATM_PERIODS[i].
 * Затем VEN -> ITRF (через vw_i) -> J2000 (через r2000), для скорости учитывается
 * вращение системы: dv = r2000*dv_itrf + dr2000_dt*dx_itrf.
 *
 * @param[in]  mjd        Модифицированная юлианская дата на 0h UTC [сут].
 * @param[in]  utc        Доля суток UTC.
 * @param[in]  pres       Атмосферное давление на станции P(MJD) [мбар].
 * @param[in]  dPdt       Производная давления dP/dt (в тех же ед./сут, что ряд; из DMETEO).
 * @param[in]  atm        Коэффициенты нагрузки станции (Station.atm_load): coef(3x8) [мм] и p_0 [мбар].
 * @param[in]  vw_i       Матрица перехода VEN (Vertical, East, North) -> ITRF (3x3).
 * @param[in]  r2000      Матрица перехода ITRF -> J2000.0 (3x3).
 * @param[in]  dr2000_dt  Первая производная r2000 по времени (3x3) [1/с].
 * @param[out] dx_atm     Вектор смещения станции в J2000.0 [м].
 * @param[out] dv_atm     Вектор скорости смещения станции в J2000.0 [м/с].
 */
void SITE_ATM40(int mjd, double utc, double pres, double dPdt,
                const AtmLoadData& atm,
                const Eigen::Matrix3d& vw_i,
                const Eigen::Matrix3d& r2000, const Eigen::Matrix3d& dr2000_dt,
                Eigen::Vector3d& dx_atm, Eigen::Vector3d& dv_atm);

/**
 * @brief Термодеформация антенны (порт THERM_DEF40, модель Nothnagel 2009).
 *
 * Тепловое расширение бетонного фундамента (высота hf) и колонны/опоры (высота hp)
 * поднимает опорную точку станции по вертикали на d_up = (gamma_hf*hf + gamma_hp*hp)*(T - T0).
 * Разность (T - T0) в °C и K одинакова. Только эти члены используются в оригинале
 * (hv/hs/AO закомментированы). Смещение задаётся в VEN (только вертикаль) и переводится
 * VEN -> ITRF (vw_i) -> J2000 (r2000); скорость учитывает вращение и dTdt (из DMETEO2_DT).
 *
 * @param[in]  tC         Локальная температура станции T [°C].
 * @param[in]  dTdt       Производная температуры dT/dt [°C/сут] (из DMETEO2_DT; 0, если не задана).
 * @param[in]  dp         Параметры антенны (ANTENNA_INFO): t_0, hf, gamma_hf, hp, gamma_hp.
 * @param[in]  vw_i       Матрица перехода VEN (Vertical, East, North) -> ITRF (3x3).
 * @param[in]  r2000      Матрица перехода ITRF -> J2000.0 (3x3).
 * @param[in]  dr2000_dt  Первая производная r2000 по времени (3x3) [1/с].
 * @param[out] dx_temp    Вектор смещения станции в J2000.0 [м].
 * @param[out] dv_temp    Вектор скорости смещения станции в J2000.0 [м/с].
 */
void THERM_DEF40(double tC, double dTdt, const DefPar& dp,
                 const Eigen::Matrix3d& vw_i,
                 const Eigen::Matrix3d& r2000, const Eigen::Matrix3d& dr2000_dt,
                 Eigen::Vector3d& dx_temp, Eigen::Vector3d& dv_temp);

// ПРИМЕЧАНИЕ: двухстанционная обёртка site_atm40(j1, j2, ...) (сборка dx_atm/dv_atm
// для пары станций наблюдения) будет добавлена вместе со слоем оркестрации, когда
// появятся посуточные P и dP/dt из DMETEO. Её черновик был удалён как несобираемый.

/**
 * @brief Вычисляет итоговую позицию, скорость и ускорение станции в небесной системе координат J2000.0.
 * Функция суммирует геоцентрические координаты станции с векторами смещений от различных 
 * геофизических эффектов (твердые приливы, океан, полюс, атмосфера, температурные деформации).
 * Результат преобразуется в J2000.0 с использованием матриц поворота и их производных.
 * Внимание: тектоническая скорость станции здесь не нужна, она учтена заранее в xsta_itrf.
 *
 * @param[in]  xsta_itrf   Вектор базовых координат станции в земной системе ITRF [м].
 * @param[in]  r2000       Матрица 3x3 перехода от ITRF к J2000.0.
 * @param[in]  dr2000      Первая производная матрицы r2000 по времени [1/с].
 * @param[in]  d2r2000     Вторая производная матрицы r2000 по времени [1/с^2].
 * @param[in]  dx_tide     Смещение коры от твердых земных приливов [м].
 * @param[in]  dv_tide     Скорость смещения коры от твердых приливов [м/с].
 * @param[in]  dx_octide   Смещение коры от океанической нагрузки [м].
 * @param[in]  dv_octide   Скорость смещения коры от океанической нагрузки [м/с].
 * @param[in]  dx_poltide  Смещение коры из-за прилива полюса [м].
 * @param[in]  dv_poltide  Скорость смещения из-за прилива полюса [м/с].
 * @param[in]  dx_atm      Смещение коры от атмосферной нагрузки [м].
 * @param[in]  dv_atm      Скорость смещения от атмосферной нагрузки [м/с].
 * @param[in]  dx_temp     Смещение из-за температурной деформации монтировки [м].
 * @param[in]  dv_temp     Скорость смещения из-за температурной деформации [м/с].
 * @param[out] xsta_j2000t Итоговый вектор координат станции в J2000.0 [м].
 * @param[out] vsta_j2000t Итоговый вектор скорости станции в J2000.0 [м/с].
 * @param[out] asta_j2000  Итоговый вектор ускорения станции в J2000.0 [м/с^2].
 */
void SITE_INST(const Eigen::Vector3d& xsta_itrf,
               const Eigen::Matrix3d& r2000, const Eigen::Matrix3d& dr2000, const Eigen::Matrix3d& d2r2000,
               const Eigen::Vector3d& dx_tide, const Eigen::Vector3d& dv_tide,
               const Eigen::Vector3d& dx_octide, const Eigen::Vector3d& dv_octide,
               const Eigen::Vector3d& dx_poltide, const Eigen::Vector3d& dv_poltide,
               const Eigen::Vector3d& dx_atm, const Eigen::Vector3d& dv_atm,
               const Eigen::Vector3d& dx_temp, const Eigen::Vector3d& dv_temp,
               Eigen::Vector3d& xsta_j2000t, Eigen::Vector3d& vsta_j2000t,
               Eigen::Vector3d& asta_j2000);

/**
 * @brief Посуточная сборка координат ПАРЫ станций в J2000.0 со всеми геофизическими поправками.
 *
 * Для каждой станции вычисляет твёрдые приливы (SITE_TIDE_SOLID), океаническую (SITE_TIDE_OC)
 * и атмосферную (SITE_ATM40) нагрузки, прилив полюса (POLE_TIDE) и через SITE_INST переводит
 * в J2000 с первой и второй производными. Термодеформация пока не учитывается (нет THERM_DEF40).
 * ВНИМАНИЕ: sun_geo/moon_geo — ГЕОЦЕНТРИЧЕСКИЕ (jpleph), не SSB.
 *
 * @param[in]  s1, s2       Подготовленные данные двух станций (SitePrep).
 * @param[in]  mjd, utc     Дата/время наблюдения.
 * @param[in]  jd           Юлианская дата на 0h [сут] (для океанических приливов).
 * @param[in]  ut1_sec      UT1 в течение суток [сек] (для океанических приливов).
 * @param[in]  cent         Юлианские столетия от J2000 (для прилива полюса).
 * @param[in]  f, fd        Фундаментальные аргументы и их производные.
 * @param[in]  gast         Гринвичское истинное звёздное время и его производная [рад, рад/с].
 * @param[in]  sun_geo      Геоцентрические поз/скор Солнца 3x2 [м, м/с].
 * @param[in]  moon_geo     Геоцентрические поз/скор Луны 3x2 [м, м/с].
 * @param[in]  xp, yp       Координаты полюса [рад] и их скорости xp_rate, yp_rate [рад/с].
 * @param[in]  r2000, dr2000, d2r2000  Матрица ITRF->J2000 и её первая/вторая производные.
 * @param[out] xsta_j2000, vsta_j2000, asta_j2000  Координаты/скорости/ускорения станций в J2000 (2 шт.).
 */
void site_pair(const SitePrep& s1, const SitePrep& s2,
               int mjd, double utc, double jd, double ut1_sec, double cent,
               const Eigen::VectorXd& f, const Eigen::VectorXd& fd,
               const Eigen::Vector2d& gast,
               const Eigen::Matrix<double, 3, 2>& sun_geo,
               const Eigen::Matrix<double, 3, 2>& moon_geo,
               double xp, double yp, double xp_rate, double yp_rate,
               const Eigen::Matrix3d& r2000, const Eigen::Matrix3d& dr2000,
               const Eigen::Matrix3d& d2r2000,
               std::vector<Eigen::Vector3d>& xsta_j2000,
               std::vector<Eigen::Vector3d>& vsta_j2000,
               std::vector<Eigen::Vector3d>& asta_j2000);

// ============================================================================
// ЭФЕМЕРИДЫ (JPL DE)
// ============================================================================

/**
 * @brief Координатное время T_eph (≈ TDB) как доля суток (роль T_EPH40, на SOFA iauDtdb).
 *
 * Даёт аргумент времени для входа в эфемериды JPL DE: ct = tt_frac + (TDB−TT)/86400.
 * Топоцентрические аргументы задают суточный член TDB−TT (нс); для геоцентра — нули.
 *
 * @param[in] jd          Юлианская дата на 0h (целая часть) [сут].
 * @param[in] tt_frac     Доля суток в шкале TT.
 * @param[in] ut1_frac    Доля суток в шкале UT1.
 * @param[in] lon_gcen    Восточная (геоцентрическая) долгота станции [рад].
 * @param[in] u_site_km   Расстояние станции от оси вращения Земли [км].
 * @param[in] v_site_km   Расстояние станции от экваториальной плоскости (Z) [км].
 * @return Доля суток в шкале TDB (координатное время T_eph).
 */
double t_eph40(double jd, double tt_frac, double ut1_frac,
               double lon_gcen, double u_site_km, double v_site_km);

/**
 * @brief Инициализация читалки эфемерид JPL DE. Вызывается один раз при старте программы.
 * @param[in] eph_file Путь к бинарному файлу эфемерид (например, "linux_p1550p2650.440t").
 */
void init_ephemeris(const std::string& eph_file);

/**
 * @brief Получение координат, скоростей и ускорений Земли, Солнца и Луны в барицентрической системе координат Солнечной Системы (SSB).
 * @param[in]  jd     Юлианская дата (целая часть).
 * @param[in]  ct     Доля дня (TDB/TT).
 * @param[out] Earth  Матрица 3x3: pos [м], vel [м/с], acc [м/с^2].
 * @param[out] Sun    Матрица 3x3: pos [м], vel [м/с], acc [м/с^2].
 * @param[out] Moon   Матрица 3x3: pos [м], vel [м/с], acc [м/с^2].
 */
void get_celestial_bodies(double jd, double ct, Eigen::Matrix3d& Earth, Eigen::Matrix3d& Sun, Eigen::Matrix3d& Moon);

/**
 * @brief Матрица перехода ITRF -> J2000(GCRS) по IAU 2006/2000A (SOFA) и её производные по времени.
 *
 * R2000 = transpose(iauC2t06a) (ITRF -> GCRS). Производные считаются центральными
 * разностями по времени (шаг 1 с); движение полюса (xp, yp) на этом интервале постоянно.
 * TDB используется как TT (расхождение <1.7 мс, ~1e-11 на матрице).
 *
 * @param[in]  JD_TDB       Юлианская дата (TDB, используется как TT).
 * @param[in]  JD_UT1       Юлианская дата в шкале UT1.
 * @param[in]  xp           Координата X полюса Земли (рад).
 * @param[in]  yp           Координата Y полюса Земли (рад).
 * @param[out] R2000        Матрица поворота 3x3 (X_J2000 = R2000 * X_ITRF).
 * @param[out] dR2000_dt    Первая производная поворота 3x3 [1/с] (для скорости).
 * @param[out] d2R2000_dt2  Вторая производная поворота 3x3 [1/с^2] (для ускорения, нужна SITE_INST).
 */
void get_r2000_matrices(double JD_TDB, double JD_UT1, double xp, double yp,
                        Eigen::Matrix3d& R2000, Eigen::Matrix3d& dR2000_dt,
                        Eigen::Matrix3d& d2R2000_dt2);

/**
 * @brief Матрица ITRF->J2000 и её первая производная (без второй производной).
 *        Совместимая версия основной перегрузки (см. выше).
 * @param[in]  JD_TDB     Юлианская дата TDB.
 * @param[in]  JD_UT1     Юлианская дата UT1.
 * @param[in]  xp, yp     Координаты полюса Земли [рад].
 * @param[out] R2000      Матрица поворота 3x3 (X_J2000 = R2000 * X_ITRF).
 * @param[out] dR2000_dt  Первая производная поворота 3x3 [1/с].
 */
void get_r2000_matrices(double JD_TDB, double JD_UT1, double xp, double yp,
                        Eigen::Matrix3d& R2000, Eigen::Matrix3d& dR2000_dt);

/**
 * @brief Гринвичское истинное звёздное время (GAST) по IAU 2006/2000A (SOFA) и его производная.
 * @param[in] jd_tt  Юлианская дата TT.
 * @param[in] jd_ut1 Юлианская дата UT1.
 * @return Vector2d: [GAST рад, d(GAST)/dt рад/с (~угловая скорость Земли)].
 */
Eigen::Vector2d gast_iau2006(double jd_tt, double jd_ut1);

/**
 * @brief Полная теоретическая задержка ОДНОГО наблюдения (связка всех блоков конвейера).
 *
 * Собирает site_pair -> baseline -> aber_source -> trop_delay -> mount_tel -> theor_delay.
 * "Средовые" величины (время, эфемериды, EOP, r2000, gast) готовит вызывающий код.
 * Солнце/Луна: для theor_delay — SSB (Earth/Sun/Moon), для приливов — геоцентрические (sun_geo/moon_geo).
 *
 * @param[in]  s1, s2      Подготовленные станции (SitePrep: ITRF, геодезия, vw, нагрузки, метео, монтировка).
 * @param[in]  K_s         Единичный вектор на источник в J2000 (source_vec).
 * @param[in]  obs         Наблюдение (метео, индексы).
 * @param[in]  mjd, utc, jd0, ct, cent, ut1_sec  Время (MJD/UTC, JD@0h, доля TDB, столетия, UT1 в сек суток).
 * @param[in]  f, fd, gast Фундаментальные аргументы и GAST (рад, рад/с).
 * @param[in]  Earth, Sun, Moon      Матрицы 3x3 в SSB (поз/скор/ускор) — для theor_delay.
 * @param[in]  sun_geo, moon_geo     Геоцентрические поз/скор 3x2 — для приливов.
 * @param[in]  xp, yp, xp_rate, yp_rate  Полюс [рад] и скорости [рад/с].
 * @param[in]  r2000, dr2000, d2r2000    ITRF->J2000 и производные.
 * @param[out] tau, dtau   Теоретическая задержка [с] и её производная [с/с].
 */
void compute_delay_obs(const SitePrep& s1, const SitePrep& s2,
                       const Eigen::Vector3d& K_s, const Observation& obs,
                       int mjd, double utc, double jd0, double ct, double cent, double ut1_sec,
                       const Eigen::VectorXd& f, const Eigen::VectorXd& fd,
                       const Eigen::Vector2d& gast,
                       const Eigen::Matrix3d& Earth, const Eigen::Matrix3d& Sun, const Eigen::Matrix3d& Moon,
                       const Eigen::Matrix<double, 3, 2>& sun_geo, const Eigen::Matrix<double, 3, 2>& moon_geo,
                       double xp, double yp, double xp_rate, double yp_rate,
                       const Eigen::Matrix3d& r2000, const Eigen::Matrix3d& dr2000, const Eigen::Matrix3d& d2r2000,
                       double& tau, double& dtau, CompDebug* dbg = nullptr, bool with_tropo = true);

// ============================================================================
// Астрометрия и матрицы переходов (Эфемериды, прецессия, нутация)
// ============================================================================

/**
 * @brief Единичные векторы направлений на радиоисточники в J2000.0.
 * @param[in]  sources  Источники (RA/Dec [рад] + собственное движение).
 * @param[in]  t_mean   Эпоха для учёта собственного движения источников [MJD].
 * @param[out] k_star   Единичные векторы на источники в J2000 (по одному на источник).
 */
void source_vec(const std::vector<Source>& sources, double t_mean, std::vector<Eigen::Vector3d>& k_star);

/**
 * @brief Интерполяция орбиты космического телескопа (RASTRON) на момент наблюдения.
 *
 * Орбита задаётся готовой эфемеридой (таблица .scf: время, положение, скорость в GCRS/J2000);
 * между точками положение/скорость/ускорение получаются кубическим сплайном. Интегратор
 * орбиты (INTEGR8asc) НЕ используется. Кадр GCRS≈J2000 берётся напрямую (без r2000).
 *
 * @param[in]  orbit    Точки орбиты (SpaceStation: mjd, utc, xyz [км], vel [км/с]), по возрастанию времени.
 * @param[in]  mjd_utc  Момент интерполяции (MJD + доля суток UTC).
 * @param[out] x_j2000  Положение в J2000 [м].
 * @param[out] v_j2000  Скорость в J2000 [м/с].
 * @param[out] a_j2000  Ускорение в J2000 [м/с^2].
 */
void orbit_interp(const std::vector<SpaceStation>& orbit, double mjd_utc,
                  Eigen::Vector3d& x_j2000, Eigen::Vector3d& v_j2000, Eigen::Vector3d& a_j2000);

/**
 * @brief Интерфейс к эфемеридам JPL DE (порт JPLEPH_421), матричный формат.
 *
 * Соглашение оригинала: Земля — барицентрическая (SSB); Солнце и Луна —
 * ГЕОЦЕНТРИЧЕСКИЕ (тело_SSB - Земля_SSB). Обёртка над get_celestial_bodies().
 *
 * @param[in]  jd    Юлианская дата на 0h UTC (целая часть) [сут].
 * @param[in]  ct    Доля координатного времени (TDB/TT) [сут].
 * @param[out] earth Матрица 3x3: SSB позиция, скорость, ускорение Земли [м, м/с, м/с^2].
 * @param[out] sun   Матрица 3x2: геоцентрические позиция и скорость Солнца [м, м/с].
 * @param[out] moon  Матрица 3x2: геоцентрические позиция и скорость Луны [м, м/с].
 */
void jpl_eph(double jd, double ct, Eigen::Matrix3d& earth, Eigen::MatrixXd& sun, Eigen::MatrixXd& moon);

// ============================================================================
// Аберрация, вектор базы, теоретическая задержка
// ============================================================================

/**
 * @brief Элевации/азимуты станций с учётом годовой и суточной аберрации источника.
 * @param[in]  obs          Наблюдение (индексы станций).
 * @param[in]  r2000        {R, dR/dt}: матрица ITRF->J2000 и её производная.
 * @param[in]  k_s          Единичный вектор на источник в J2000.
 * @param[in]  earth        Скорость Земли (SSB) в столбце 1 — для годовой аберрации [м/с].
 * @param[in]  vsta_j2000t  Скорости станций в J2000 [м/с] (суточная аберрация).
 * @param[in]  vw           Матрицы VEN->ITRF по станциям (локальный горизонт).
 * @param[out] e            Элевации: 2x2 (станция × значение/скорость) [рад, рад/с].
 * @param[out] az           Азимуты:  2x2 (станция × значение/скорость) [рад, рад/с].
 */
void aber_source(const Observation& obs, const std::vector<Eigen::Matrix3d>& r2000, const Eigen::Vector3d& k_s, const Eigen::Matrix<double, 3, 3>& earth, const std::vector<Eigen::Vector3d>& vsta_j2000t, const std::vector<Eigen::Matrix3d>& vw, Eigen::Matrix2d& e, Eigen::Matrix2d& az);

/**
 * @brief Вектор базы в инерциальной (J2000) и земной (crust-fixed) системах координат.
 * @param[in]  r2000        Матрица ITRF->J2000.
 * @param[in]  xsta_j2000t  Координаты двух станций в J2000, 3x2 [м].
 * @param[in]  vsta_j2000t  Скорости двух станций в J2000, 3x2 [м/с].
 * @param[in]  asta_j2000   Ускорения двух станций в J2000, 3x2 [м/с^2].
 * @param[out] base_line    Вектор базы в J2000, 3x2 (значение/скорость) [м, м/с].
 * @param[out] b_cfs        Вектор базы в земной (crust-fixed) системе [м].
 */
void baseline(const Eigen::Matrix3d& r2000, const Eigen::MatrixXd& xsta_j2000t, const Eigen::MatrixXd& vsta_j2000t, const Eigen::MatrixXd& asta_j2000, Eigen::MatrixXd& base_line, Eigen::Vector3d& b_cfs);

/**
 * @brief Верхний слой оркестрации: пакетный расчёт теоретической задержки всех наблюдений сессии.
 *
 * Порт основного цикла ARIADNA4_5corr (без слоя оценивания der_* и create_matr — вне области переноса).
 * ОДИН РАЗ на сессию: site() (геодезия + VEN->ITRF на среднюю эпоху) и source_vec() (направления).
 * НА КАЖДОЕ НАБЛЮДЕНИЕ: время (tai_time/t_eph40), EOP (окно 7 узлов + interp_iers), эфемериды
 * (get_celestial_bodies/jpl_eph), ориентация Земли (get_r2000_matrices/fund_arg/gast_iau2006),
 * затем весь конвейер задержки compute_delay_obs().
 *
 * РЕЗУЛЬТАТ ВОЗВРАЩАЕТСЯ В ПАМЯТИ (results) — БЕЗ промежуточных файлов. Запись в файл
 * output_path выполняется ТОЛЬКО если путь непустой (необязательный итоговый отчёт).
 * Опционально заполняет debug (промежуточные величины конвейера на каждое наблюдение).
 *
 * @param[in]  stations       Станции с координатами ITRF@эпоха каталога, скоростями, нагрузками, монтировкой.
 * @param[in]  sources        Источники (RA/Dec в рад).
 * @param[in]  observations   Наблюдения (индексы станций/источника, метео, эпоха).
 * @param[in]  space_stations Данные бортовых станций (Space VLBI); передаются в site().
 * @param[in]  orbit_data     Орбитальные данные (Space VLBI); не задействованы.
 * @param[in]  n_segm,k_ch_c,k_ch_z,delta_sec  Параметры сегментации вывода; не задействованы.
 * @param[in]  output_path    Путь к необязательному файлу-отчёту ("" = не писать файл).
 * @param[in]  eop_data       Полная таблица EOP (EOPData, отсортирована по MJD); окно выбирается на наблюдение.
 * @param[in]  mjd_mean,utc_mean  Средняя эпоха сессии (для дрейфа координат в site()).
 * @param[in]  t_mean         Средняя эпоха для собственного движения источников (source_vec).
 * @param[out] results        Задержки всех наблюдений (В ПАМЯТИ).
 * @param[out] debug          (опц.) промежуточные величины на наблюдение; nullptr = не считать.
 */
void process_ariadna(const std::vector<Station>& stations, const std::vector<Source>& sources, const std::vector<Observation>& observations, const std::vector<SpaceStation>& space_stations, const std::vector<OrbitData>& orbit_data, int n_segm, int k_ch_c, int k_ch_z, double delta_sec, const std::string& output_path, const std::vector<EOPData>& eop_data, double mjd_mean, double utc_mean, double t_mean, std::vector<DelayResult>& results, std::vector<CompDebug>* debug = nullptr);

// ============================================================================
// Ввод задания коррелятора (.cfx) и орбиты (.scf)
// ============================================================================

/**
 * @brief Разбор файла задания коррелятора (.cfx) в структуру CfxTask.
 *
 * Читает блоки [$TLSC] (станции), [$SOURCE] (источники), [$skan] (сканы) и [$OUTPAR].
 * Наземные станции берут координаты из TLSC_PAR; космическая (RASTRON) — из ORB_FILE.
 *
 * @param[in]  path  Путь к файлу .cfx.
 * @param[out] task  Заполненное задание (станции, источники, сканы, выходная папка).
 * @return true при успехе (есть станции и сканы).
 */
bool parse_cfx(const std::string& path, CfxTask& task);

/**
 * @brief Чтение орбиты космического телескопа из .scf (CCSDS OEM, EME2000≈J2000, UTC).
 *
 * @param[in]  path   Путь к файлу .scf.
 * @param[out] orbit  Точки орбиты (SpaceStation: mjd, utc, xyz [км], vel [км/с]).
 * @return true при успехе (орбита непуста).
 */
bool read_scf_orbit(const std::string& path, std::vector<SpaceStation>& orbit);

// ============================================================================
// Полиномы задержки для коррелятора
// ============================================================================

/**
 * @brief Полиномы задержки одной НАЗЕМНОЙ станции по сеансу (сшитые кубическим сплайном).
 *
 * Считает геоцентрическую ВАКУУМНУЮ задержку станции (относительно центра Земли,
 * без тропосферы и clock) на сетке с шагом sample_sec, строит единый кубический сплайн
 * (C2-непрерывность = сшивка) и раскладывает по блокам block_sec в полиномы степени degree.
 *
 * @param[in]  st         Станция задания (координаты ITRF, эпоха, монтировка).
 * @param[in]  src        Источник (RA/Dec [рад]).
 * @param[in]  mjd0,utc0  Начало сеанса (MJD, доля суток UTC).
 * @param[in]  dur_sec    Длительность сеанса [с].
 * @param[in]  eop        Узлы EOP (7 шт.) вокруг сеанса.
 * @param[in]  block_sec  Длительность блока полинома [с] (эталон = 60).
 * @param[in]  degree     Степень полинома (эталон = 5, т.е. order=6).
 * @param[in]  sample_sec Шаг сетки расчёта задержки [с] (плотнее блока).
 * @param[in]  orbit      Орбита (для космической станции, is_space); пусто для наземной.
 * @return Полиномы станции (StationPoly).
 */
StationPoly compute_station_poly(const CfxStation& st, const CfxSource& src,
                                 int mjd0, double utc0, double dur_sec,
                                 const std::vector<EOPData>& eop,
                                 double block_sec, int degree, double sample_sec,
                                 const std::vector<SpaceStation>& orbit);

/**
 * @brief Запись полиномов задержки станции в текстовый файл .TXT (формат эталона).
 * @param[in] path  Путь к выходному файлу.
 * @param[in] poly  Полиномы станции.
 */
void write_station_poly(const std::string& path, const StationPoly& poly);

/**
 * @brief ГОТОВЫЙ МОДУЛЬ: по заданию .cfx и орбите .scf считает полиномы задержки всех
 *        станций сеанса и пишет файлы .TXT (формат коррелятора).
 *
 * Читает задание (станции/источник/сканы), орбиту космической станции, узлы EOP из
 * каталога; определяет границы сеанса; для каждой станции считает сшитые полиномы
 * геоцентрической вакуумной задержки и пишет файл (имя из POLY_FILE) в out_dir.
 * Эфемериды должны быть инициализированы заранее (init_ephemeris).
 *
 * @param[in] cfx_path   Путь к файлу задания .cfx.
 * @param[in] orbit_path Путь к файлу орбиты .scf (для космической станции; "" если нет).
 * @param[in] out_dir    Папка для выходных файлов полиномов.
 * @param[in] eop_path   Путь к каталогу EOP (EOPC04).
 * @param[in] block_sec  Длительность блока [с] (по умолчанию 60).
 * @param[in] degree     Степень полинома (по умолчанию 5).
 * @param[in] sample_sec Шаг сетки задержки [с] (по умолчанию 6).
 */
void process_task(const std::string& cfx_path, const std::string& orbit_path,
                  const std::string& out_dir, const std::string& eop_path,
                  double block_sec = 60.0, int degree = 5, double sample_sec = 6.0);

} // namespace ariadna