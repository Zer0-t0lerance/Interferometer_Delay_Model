// Catalog Parsing Functions
// Version: 19.08.2021
//
// EOP, Atmosphere, Ocean, VTRF, ITRF, ICRF, ICRF_NON_VCS, ICRF_VCS, ALTNAMES, ANTENNA_INFO
#pragma once

// Catalog Parsing Functions
// Version: 19.08.2021
#include <vector>
#include <fstream>

using namespace std;

//-----------------------------------------------------------------------------
/* Variables */
//-----------------------------------------------------------------------------

class eop_record
{
public:
	int year, month, day;
	double MJD, x, y, ut1_utc, lod, dpsi, deps, x_err, y_err, ut1_utc_err, lod_err, dpsi_err, depsilon_err;

	eop_record()
	{
		year = 0; month = 0; day = 0;
		MJD = 0; 
		x = 0;
		y = 0;
		ut1_utc = 0;
		lod = 0;
		dpsi = 0;
		deps = 0;
		x_err = 0;
		y_err = 0;
		ut1_utc_err = 0;
		lod_err = 0;
		dpsi_err = 0;
		depsilon_err = 0;
	}
};

class oc_record
{
public:
	char telescope[256];
	int telescope_id;

	double coeff1[33];
	double coeff2[33];

	oc_record()
	{
		telescope_id = 0;
	}
};

class atm_record
{
public:
	char telescope[256];
	
	double coeff[30];

};

class ant_record
{
public:
	/* string parameters */
	char telescope[8]; /* station name */
	char focus_type[7]; /* focus type */
	char mount_type[7]; /* mount type MO_AZEL(azimuthal), MO_EQUA(equatorial),
							MO_XYNO(XY north), MO_XYEA(XY east),
							MO_RICH(misplaced equatorial RICHMOND) */
	char radome[6];		/* Flag whether the station has radome : RA_NO or RA_YES */
	char m_type[6];		/* Measurement type ME_COMP(complete) ME_INCM(incomplete) ME_ROUG(rough) */
	/* numerical parameters */
	double ref_temp;	/* Reference temperature (degree C) */
	
	double sin_amp;	/* sin amplitude of temp variations (degrees C), J2000.0 */
	double cos_amp;	/* cos amplitude of temp variations (degrees C), J2000.0 */
		
	double ref_press;	/* Reference pressure (hPa) */

	double ant_diam;	/* Antenna diameter (m) */
	double ant_h;		/* Height of foundation(m) */
	double ant_d;		/* Depth of foundation(m) */

	double th_exp_f;	/* Foundation thermal expansion coefficient(1 / K) */
	double f_a_len;	/* Length of the fixed axis(m) */
	double th_exp_f_a;	/* Fixed axis thermal expansion coefficient(1 / K) */
	double a_of_len;	/* Length of the axis offset(m) */

	double a_of_thermal;	/* 137:143 -- Axis offset thermal expansion coefficient(1 / K) */
	double dist_mov_a;		/* 146:152 -- Distance from the movable axis to the antenna vertex(m) */
	double therm_exp;		/* 154:160 -- Thermal expansion coefficient of the structure from 
										  the movable axis to the antenna vertex(1 / K) */
	double sub_h;			/* 163:169 --Height of the sub - reflector / primary focus above the vertex(m) */
	double sub_exp;        /* 171:177 --Sub - reflector / primary focus mounting legs thermal expansion coefficient */
};

class src_icrf_record
{
public:
	char src_name_j2000[16];	/* J2000 source name */
	char src_name_1950[8];		/* B1950 source name */

	char inf;

	/* coordinates */
	double ra_h, ra_m, ra_s;
	double dec_d, dec_m, dec_s;

	double ra_err;
	double dec_err;

	double ra_dec_corr;

	double MJD_mean;
	double MJD_first;
	double MJD_last;

	double Nb_sess;
	double Nb_del;
	double Nb_rat;
};

class src_non_vcs_record
{
public:
	char src_name_j2000[16];	/* J2000 source name */
	char src_name_1950[8];		/* B1950 source name */

	char inf;

	/* coordinates */
	double ra_h, ra_m, ra_s;
	double dec_d, dec_m, dec_s;

	double ra_err;
	double dec_err;

	double ra_dec_corr;

	double MJD_mean;
	double MJD_first;
	double MJD_last;

	double Nb_sess;
	double Nb_del;
};

class src_vcs_record
{
public:
	char src_name_j2000[16];	/* J2000 source name */
	char src_name_1950[8];		/* B1950 source name */

	/* coordinates */
	double ra_h, ra_m, ra_s;
	double dec_d, dec_m, dec_s;

	double ra_err;
	double dec_err;

	double ra_dec_corr;

	double MJD_mean;
	double MJD_first;
	double MJD_last;

	double Nb_sess;
	double Nb_del;
};

class alt_record
{
public:
	char alt_name[256];
	char src_name[256];
};

class ant_vtrf_record
{
public:
	char domes_nb[9];
	char ant_name[9];
	
	char tech_id[4];

	double x, y, z;
	double sx, sy, sz;
	double vx, vy, vz;
	double svx, svy, svz;
};

class ant_itrf_record
{
public:
	char domes_nb[9];
	char ant_name[9];

	char tech_id[4];

	char mount_type[4];
	double axis_offset;
};
//-----------------------------------------------------------------------------
//									 Functions								*/
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//                                   Functions                               */
//-----------------------------------------------------------------------------

/**
 * @brief Считывает каталог параметров вращения Земли (EOP).
 * @param[out] EOP_Data      Вектор структур для сохранения данных EOP.
 * @param[in]  filename_EOP  Путь к файлу каталога (например, "external/catalogs/EOP.dat").
 * @return 1 при успешном чтении, -1 если файл не найден.
 */
int ReadEOP(vector<eop_record> &EOP_Data, char filename_EOP[256]);

/**
 * @brief Считывает каталог атмосферной нагрузки.
 * @param[out] ATM_Data      Вектор структур для сохранения данных атмосферной нагрузки.
 * @param[in]  filename_ATM  Путь к файлу каталога.
 * @return 1 при успешном чтении, -1 если файл не найден.
 */
int ReadATM(vector<atm_record> &ATM_Data, char filename_ATM[256]);

/**
 * @brief Считывает каталог океанической нагрузки (амплитуды и фазы 11 приливных волн).
 * @param[out] OC_Data       Вектор структур oc_record (каждая хранит массивы coeff1 и coeff2 по 33 элемента).
 * @param[in]  filename_OC   Путь к файлу каталога (например, "external/catalogs/ocload40.dat").
 * @return 1 при успешном чтении, -1 если файл не найден.
 */
int ReadOC(vector<oc_record> &OC_Data, char filename_OC[256]);

/**
 * @brief Считывает каталог информации об антеннах (фокус, монтировка, размеры).
 * @param[out] ANT_INFO_Data Вектор структур с технической информацией о телескопах.
 * @param[in]  filename_ANT  Путь к файлу каталога.
 * @return 1 при успешном чтении, -1 если файл не найден.
 */
int ReadANT_INFO(vector<ant_record> &ANT_INFO_Data, char filename_ANT[256]);

/**
 * @brief Считывает каталог опорных координат и скоростей станций (VTRF).
 * @param[out] ANT_VTRF_Data Вектор структур с координатами (X, Y, Z) и скоростями.
 * @param[in]  filename_VTRF Путь к файлу каталога.
 * @return 1 при успешном чтении, -1 если файл не найден.
 */
int ReadVTRF(vector<ant_vtrf_record> &ANT_VTRF_Data, char filename_VTRF[256]);

/**
 * @brief Считывает каталог ITRF (типы монтировок, смещения осей).
 * @param[out] ANT_ITRF_Data Вектор структур с параметрами ITRF.
 * @param[in]  filename_ITRF Путь к файлу каталога.
 * @return 1 при успешном чтении, -1 если файл не найден.
 */
int ReadITRF(vector<ant_itrf_record> &ANT_ITRF_Data, char filename_ITRF[256]);

int ReadICRF(vector<src_icrf_record> &SRC_ICRF_Data, char filename_ICRF[256]);
int ReadICRF_NON_VCS(vector<src_non_vcs_record> &SRC_NON_VCS_Data, char filename_ICRF_NON_VCS[256]);
int ReadICRF_VCS(vector<src_vcs_record> &SRC_Data, char filename_ICRF_VCS[256]);

int ReadALT(vector<alt_record> &ALT_Data, char filename_ALT[256]);



