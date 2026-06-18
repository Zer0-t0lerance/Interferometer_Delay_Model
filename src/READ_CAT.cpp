// Catalog Parsing Functions
// Version: 19.08.2021
//
// EOP, Atmosphere, Ocean, VTRF, ITRF, ICRF, ICRF_NON_VCS, ICRF_VCS, ALTNAMES, ANTENNA_INFO
#include "READ_CAT.h"
#include <cstring>   // <--- ДОБАВИТЬ ЭТУ СТРОКУ
#include <iostream>  // Пригодится для вывода ошибок

using namespace std;

//-----------------------------------------------------------------------------
// EOP catalog parsing function 
//-----------------------------------------------------------------------------
/* Format -- ASCII text

First 14 lines are obsolete

Column values:
Date (yyyy mm dd), MJD, x, y, UT1-UTC (s), LOD (s), dPsi, dEps, x Err, y Err, UT1-UTC Err (s), LOD Err (s), dPsi Err, dEpsilon Err
Initial time moment for each calculation - 0h UTC
*/
//-----------------------------------------------------------------------------
int ReadEOP(vector<eop_record> &EOP_Data, char filename_EOP[256])
{
	fstream fs;
	char data[512];

	eop_record temp_record;
	EOP_Data.clear();

	fs.open(filename_EOP, fstream::in);

	if (fs.is_open())
	{
		// Read EOP files
		for (int i = 0; i < 14; i++) //skip 14 lines
		{
			fs.getline(data, 512);
		}

		while (!fs.eof())
		{
			fs.getline(data, 512);
			sscanf(data,"%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf ", &temp_record.year, &temp_record.month, &temp_record.day, &temp_record.MJD, &temp_record.x, &temp_record.y, &temp_record.ut1_utc, &temp_record.lod, &temp_record.dpsi, &temp_record.deps, &temp_record.x_err, &temp_record.y_err, &temp_record.ut1_utc_err, &temp_record.lod_err, &temp_record.dpsi_err, &temp_record.depsilon_err);
			EOP_Data.push_back(temp_record);
		}
		fs.close();
		return 1; // Parsing successful
	}
	else
	{
		return -1;
	}
	return 0; 
}

//-----------------------------------------------------------------------------
// Atmosphere load catalogue
//-----------------------------------------------------------------------------
/* Format -- ASCII text
Column values:
*/
//-----------------------------------------------------------------------------
int ReadATM(vector<atm_record> &ATM_Data, char filename_ATM[256])
{
	fstream fs;
	char data[512];
	
	atm_record temp_record;
	ATM_Data.clear();

	fs.open(filename_ATM, fstream::in);

	if (fs.is_open())
	{
		// Read atmosphere files
		for (int i = 0; i < 26; i++) //skip 26 lines
		{
			fs.getline(data, 512);
		}

		while (!fs.eof())
		{
			fs.getline(data, 512); // telescope name
			sscanf(data, "%s", temp_record.telescope);
			fs.getline(data, 512); // first line of coefficients
			sscanf(data, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &temp_record.coeff[0], &temp_record.coeff[1], &temp_record.coeff[2], &temp_record.coeff[3], &temp_record.coeff[4], &temp_record.coeff[5], &temp_record.coeff[6], &temp_record.coeff[7], &temp_record.coeff[8], &temp_record.coeff[9]);
			fs.getline(data, 512); // second line of coefficients
			sscanf(data, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &temp_record.coeff[10], &temp_record.coeff[11], &temp_record.coeff[12], &temp_record.coeff[13], &temp_record.coeff[14], &temp_record.coeff[15], &temp_record.coeff[16], &temp_record.coeff[17], &temp_record.coeff[18], &temp_record.coeff[19]);
			fs.getline(data, 521); // last line of coefficients
			sscanf(data, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &temp_record.coeff[20], &temp_record.coeff[21], &temp_record.coeff[22], &temp_record.coeff[23], &temp_record.coeff[24], &temp_record.coeff[25], &temp_record.coeff[26], &temp_record.coeff[27], &temp_record.coeff[28], &temp_record.coeff[29]);
			fs.getline(data, 512); // empty line
			ATM_Data.push_back(temp_record);
		}
		fs.close();
		return 1; // Parsing successful
	}
	else
	{
		return -1;
	}
	return 0;
}


/* --------------------------------------------------------- */
/*  Read catalogue with ocean load information ocload40.dat  */
/* --------------------------------------------------------- */
/*															 */
/* Format: */
int ReadOC(vector<oc_record> &OC_Data, char filename_OC[256])
{
	fstream fs;
	char data[512];

	oc_record temp_record;
	OC_Data.clear();

	fs.open(filename_OC, fstream::in);

	if (fs.is_open())
	{
		// Read atmosphere files
		for (int i = 0; i < 30; i++) //skip 30 lines
		{
			fs.getline(data, 512);
		}

		while (!fs.eof())
		{
			fs.getline(data, 512); // telescope name and id
			sscanf(data, "%s %d", temp_record.telescope, &temp_record.telescope_id);
			
			fs.getline(data, 512); // commented line
			fs.getline(data, 512); // commented line
			fs.getline(data, 512); // commented line

			fs.getline(data, 512); // first line of coefficients 1
			sscanf(data, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &temp_record.coeff1[0], &temp_record.coeff1[1], &temp_record.coeff1[2], &temp_record.coeff1[3], &temp_record.coeff1[4], &temp_record.coeff1[5], &temp_record.coeff1[6], &temp_record.coeff1[7], &temp_record.coeff1[8], &temp_record.coeff1[9], &temp_record.coeff1[10]);
			fs.getline(data, 512); // second line of coefficients 1
			sscanf(data, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &temp_record.coeff1[11], &temp_record.coeff1[12], &temp_record.coeff1[13], &temp_record.coeff1[14], &temp_record.coeff1[15], &temp_record.coeff1[16], &temp_record.coeff1[17], &temp_record.coeff1[18], &temp_record.coeff1[19], &temp_record.coeff1[20], &temp_record.coeff1[21]);
			fs.getline(data, 512); // third line of coefficients 1
			sscanf(data, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &temp_record.coeff1[22], &temp_record.coeff1[23], &temp_record.coeff1[24], &temp_record.coeff1[25], &temp_record.coeff1[26], &temp_record.coeff1[27], &temp_record.coeff1[28], &temp_record.coeff1[29], &temp_record.coeff1[30], &temp_record.coeff1[31], &temp_record.coeff1[32]);

			fs.getline(data, 512); // first line of coefficients 2
			sscanf(data, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &temp_record.coeff2[0], &temp_record.coeff2[1], &temp_record.coeff2[2], &temp_record.coeff2[3], &temp_record.coeff2[4], &temp_record.coeff2[5], &temp_record.coeff2[6], &temp_record.coeff2[7], &temp_record.coeff2[8], &temp_record.coeff2[9], &temp_record.coeff2[10]);
			fs.getline(data, 512); // second line of coefficients 2
			sscanf(data, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &temp_record.coeff2[11], &temp_record.coeff2[12], &temp_record.coeff2[13], &temp_record.coeff2[14], &temp_record.coeff2[15], &temp_record.coeff2[16], &temp_record.coeff2[17], &temp_record.coeff2[18], &temp_record.coeff2[19], &temp_record.coeff2[20], &temp_record.coeff2[21]);
			fs.getline(data, 512); // third line of coefficients 2
			sscanf(data, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &temp_record.coeff2[22], &temp_record.coeff2[23], &temp_record.coeff2[24], &temp_record.coeff2[25], &temp_record.coeff2[26], &temp_record.coeff2[27], &temp_record.coeff2[28], &temp_record.coeff2[29], &temp_record.coeff2[30], &temp_record.coeff2[31], &temp_record.coeff2[32]);
			fs.getline(data, 512); // commented line
			OC_Data.push_back(temp_record);
		}
		fs.close();
		return 1; // Parsing successful
	}
	else
	{
		return -1;
	}
	return 0;
}

/* --------------------------------------------------------- */
/*    Read catalogue with antenna information ANTENNA_INFO   */
/* --------------------------------------------------------- */
/*															 */
/* Format:

   1)   1:12  -- Record label : word ANTENNA_INFO
   2)  15:22  -- IVS station name
   3)  25:31  -- Focus type : FO_PRIM -- primary, FO_SECN -- secondary
                 In the case of separate foci for S and X band we refer the
                 focus type to the primary frequency, i.e.X band at this time
   4)  33:39  -- Mounting type : MO_AZEL(azimuthal), MO_EQUA(equatorial),
                 MO_XYNO(XY north), MO_XYEA(XY east),
                 MO_RICH(misplaced equatorial RICHMOND)
   5)  41:47  -- Flag, whether the station has radome : RA_NO or RA_YES
   6)  49:55  -- Measurement type ME_COMP(complete) ME_INCM(incomplete)
                 ME_ROUG(rough)
   7)  58:61  -- Reference temperature(degree C)
   8)  63:66  -- sin amplitude of annual temperature variations with respect
                 to the J2000.0 epoch(degree C)
   9)  68:71  -- cos amplitude of annual temperature variations with respect
                 to the J2000.0 epoch(degree C)
  10)  73:78  -- Reference pressure(hPa)
  11)  81:85  -- Antenna diameter(m)
  12)  87:93  -- Height of foundation(m)
  13)  95:100 -- Depth of foundation(m)
  14) 103:109 -- Foundation thermal expansion coefficient(1 / K)
  15) 112:118 -- Length of the fixed axis(m)
  16) 120:126 -- Fixed axis thermal expansion coefficient(1 / K)
  17) 129:135 -- Length of the axis offset(m)
  18) 137:143 -- Axis offset thermal expansion coefficient(1 / K)
  19) 146:152 -- Distance from the movable axis to the antenna vertex(m)
  20) 154:160 -- Thermal expansion coefficient of the structure from
                 the movable axis to the antenna vertex(1 / K)
  21) 163:169 -- Height of the sub - reflector / primary focus above the vertex(m)
  22) 171:177 -- Sub-reflector / primary focus mounting legs thermal expansion coefficient(1 / K)
*/
int ReadANT_INFO(vector<ant_record> &ANT_INFO_Data, char filename_ANT[256])
{
	fstream fs;
	char data[512];

	ant_record temp_record;
	ANT_INFO_Data.clear();

	fs.open(filename_ANT, fstream::in);

	if (fs.is_open())
	{
		while (!fs.eof())
		{
			fs.getline(data, 512); // telescope name and id
			if (strncmp(data, "ANTENNA_INFO",12)==0)
			{
				sscanf(data, "%*s %s %s %s %s %s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", temp_record.telescope, temp_record.focus_type, temp_record.mount_type, temp_record.radome, temp_record.m_type, &temp_record.ref_temp, &temp_record.sin_amp, &temp_record.cos_amp, &temp_record.ref_press, &temp_record.ant_diam, &temp_record.ant_h, &temp_record.ant_d, &temp_record.th_exp_f, &temp_record.f_a_len, &temp_record.th_exp_f_a, &temp_record.a_of_len, &temp_record.a_of_thermal, &temp_record.dist_mov_a, &temp_record.therm_exp, &temp_record.sub_h, &temp_record.sub_exp);
				ANT_INFO_Data.push_back(temp_record);
			}
		}
		return 1;
	}
	else return -1;

	return 0;
}
/* --------------------------------------------------------- */
/* Read catalogue of alternative source names ALT_SOURCE.cat */
/* Format: [alt_name] [source_name_b1950]					 */
/* --------------------------------------------------------- */
int ReadALT(vector<alt_record> &ALT_Data, char filename_ALT[256])
{
	fstream fs;
	char data[512];

	alt_record temp_record;
	ALT_Data.clear();

	fs.open(filename_ALT, fstream::in);

	if (fs.is_open())
	{
		while (!fs.eof())
		{
			fs.getline(data, 512); // telescope name and id
			sscanf(data, "%s %s", temp_record.alt_name, temp_record.src_name);
			ALT_Data.push_back(temp_record);
		}
		return 1;
	}
	else return -1;

	return 0;
}
/* --------------------------------------------------------- */
/* Read source catalogue ICRF VCS		ICRF-vcs-only.cat    */
/* Format:													 */
/* --------------------------------------------------------- */
int ReadICRF_VCS(vector<src_vcs_record> & SRC_VCS_Data, char filename_ICRF_VCS[256])
{
	fstream fs;
	char data[512];
	char ra_h_tmp[2];
	char ra_m_tmp[2];
	char ra_s_tmp[11];

	char dec_d_tmp[3];
	char dec_m_tmp[2];
	char dec_s_tmp[10];

	char ra_err_tmp[10];
	char dec_err_tmp[10];

	char ra_dec_corr_tmp[7];

	char mjd_mean_tmp[7];
	char mjd_first_tmp[7];
	char mjd_last_tmp[7];

	char nb_sess_tmp[6];
	char nb_del_tmp[6];

	src_vcs_record temp_record;
	SRC_VCS_Data.clear();

	fs.open(filename_ICRF_VCS, fstream::in);

	if (fs.is_open())
	{
		while (!fs.eof())
		{
			fs.getline(data, 512);								/* telescope name and id */

			strncpy(temp_record.src_name_j2000, data + 5, 16);	/* source name J2000 */
			strncpy(temp_record.src_name_1950, data + 23, 8);	/* source name B1950 */
				
			strncpy(ra_h_tmp, data + 33, 2);					/* right ascension hh */
			sscanf(ra_h_tmp, "%lf", &temp_record.ra_h);
			strncpy(ra_m_tmp, data + 36, 2);					/* right ascension mm */
			sscanf(ra_m_tmp, "%lf", &temp_record.ra_m);
			strncpy(ra_s_tmp, data + 39, 11);					/* right ascension ss */
			sscanf(ra_s_tmp, "%lf", &temp_record.ra_s);

			strncpy(dec_d_tmp, data + 52, 3);					/* declination degrees */
			sscanf(dec_d_tmp, "%lf", &temp_record.dec_d);
			strncpy(dec_m_tmp, data + 56, 2);					/* declination mm */
			sscanf(dec_m_tmp, "%lf", &temp_record.dec_m);
			strncpy(dec_s_tmp, data + 59, 10);					/* declination ss */
			sscanf(dec_s_tmp, "%lf", &temp_record.dec_s);

			strncpy(ra_err_tmp, data + 71, 10);					/* right ascension err ss */
			sscanf(ra_err_tmp, "%lf", &temp_record.ra_err);

			strncpy(dec_err_tmp, data + 82, 9);					/* declination err ss */
			sscanf(dec_err_tmp, "%lf", &temp_record.dec_err);

			strncpy(ra_dec_corr_tmp, data + 93, 7);			/* RA-DEC correction */
			sscanf(ra_dec_corr_tmp, "%lf", &temp_record.ra_dec_corr);
			
			strncpy(mjd_mean_tmp, data + 101, 7);			/* MJD mean */
			sscanf(mjd_mean_tmp, "%lf", &temp_record.MJD_mean);
			strncpy(mjd_first_tmp, data + 109, 7);			/* MJD first */
			sscanf(mjd_first_tmp, "%lf", &temp_record.MJD_first);
			strncpy(mjd_last_tmp, data + 117, 7);			/* MJD last */
			sscanf(mjd_last_tmp, "%lf", &temp_record.MJD_last);

			strncpy(nb_sess_tmp, data + 125, 6);			/* Nb sessions */
			sscanf(nb_sess_tmp, "%lf", &temp_record.Nb_sess);
			strncpy(nb_del_tmp, data + 132, 6);			/* Nb del */
			sscanf(nb_del_tmp, "%lf", &temp_record.Nb_del);

			SRC_VCS_Data.push_back(temp_record);
		}
		return 1;
	}
	else return -1;

	return 0;
}
/* --------------------------------------------------------- */
/* Read source catalogue ICRF NON VCS		ICRF-non-vcs.cat */
/* Format:													 */
/* --------------------------------------------------------- */
int ReadICRF_NON_VCS(vector<src_non_vcs_record> & SRC_NON_VCS_Data, char filename_ICRF_NON_VCS[256])
{
	fstream fs;
	char data[512];
	char ra_h_tmp[2];
	char ra_m_tmp[2];
	char ra_s_tmp[11];

	char dec_d_tmp[3];
	char dec_m_tmp[2];
	char dec_s_tmp[10];

	char ra_err_tmp[10];
	char dec_err_tmp[10];

	char ra_dec_corr_tmp[7];

	char mjd_mean_tmp[7];
	char mjd_first_tmp[7];
	char mjd_last_tmp[7];

	char nb_sess_tmp[6];
	char nb_del_tmp[6];

	src_non_vcs_record temp_record;
	SRC_NON_VCS_Data.clear();

	fs.open(filename_ICRF_NON_VCS, fstream::in);

	if (fs.is_open())
	{
		while (!fs.eof())
		{
			fs.getline(data, 512); // telescope name and id
			if (strncmp(data, "ICRF J", 6) == 0)					/* search for the first line of data */
			{
				strncpy(temp_record.src_name_j2000, data + 5, 16);	/* source name J2000 */
				strncpy(temp_record.src_name_1950, data + 23, 8);	/* source name B1950 */
				temp_record.inf = data[33];

				strncpy(ra_h_tmp, data + 36, 2);					/* right ascension hh */
				sscanf(ra_h_tmp, "%lf", &temp_record.ra_h);
				strncpy(ra_m_tmp, data + 39, 2);					/* right ascension mm */
				sscanf(ra_m_tmp, "%lf", &temp_record.ra_m);
				strncpy(ra_s_tmp, data + 42, 11);					/* right ascension ss */
				sscanf(ra_s_tmp, "%lf", &temp_record.ra_s);

				strncpy(dec_d_tmp, data + 55, 3);					/* declination degrees */
				sscanf(dec_d_tmp, "%lf", &temp_record.dec_d);
				strncpy(dec_m_tmp, data + 59, 2);					/* declination mm */
				sscanf(dec_m_tmp, "%lf", &temp_record.dec_m);
				strncpy(dec_s_tmp, data + 62, 10);					/* declination ss */
				sscanf(dec_s_tmp, "%lf", &temp_record.dec_s);

				strncpy(ra_err_tmp, data + 75, 10);					/* right ascension err ss */
				sscanf(ra_err_tmp, "%lf", &temp_record.ra_err);

				strncpy(dec_err_tmp, data + 85, 9);					/* declination err ss */
				sscanf(dec_err_tmp, "%lf", &temp_record.dec_err);

				strncpy(ra_dec_corr_tmp, data + 96, 7);			/* RA-DEC correction */
				sscanf(ra_dec_corr_tmp, "%lf", &temp_record.ra_dec_corr);

				strncpy(mjd_mean_tmp, data + 104, 7);			/* MJD mean */
				sscanf(mjd_mean_tmp, "%lf", &temp_record.MJD_mean);

				strncpy(mjd_first_tmp, data + 112, 7);			/* MJD first */
				sscanf(mjd_first_tmp, "%lf", &temp_record.MJD_first);

				strncpy(mjd_last_tmp, data + 120, 7);			/* MJD last */
				sscanf(mjd_last_tmp, "%lf", &temp_record.MJD_last);

				strncpy(nb_sess_tmp, data + 128, 6);			/* Nb sessions */
				sscanf(nb_sess_tmp, "%lf", &temp_record.Nb_sess);

				strncpy(nb_del_tmp, data + 135, 6);			/* Nb del */
				sscanf(nb_del_tmp, "%lf", &temp_record.Nb_del);

				SRC_NON_VCS_Data.push_back(temp_record);
			}
		}
		return 1;
	}
	else return -1;

	return 0;
}
/* --------------------------------------------------------- */
/* Read source catalogue ICRF        				ICRF.cat */
/* Format:													 */
/* --------------------------------------------------------- */
int ReadICRF(vector<src_icrf_record> &SRC_ICRF_Data, char filename_ICRF[256])
{
	fstream fs;
	char data[512];
	char ra_h_tmp[2];
	char ra_m_tmp[2];
	char ra_s_tmp[11];

	char dec_d_tmp[3];
	char dec_m_tmp[2];
	char dec_s_tmp[10];

	char ra_err_tmp[10];
	char dec_err_tmp[10];

	char ra_dec_corr_tmp[7];

	char mjd_mean_tmp[7];
	char mjd_first_tmp[7];
	char mjd_last_tmp[7];

	char nb_sess_tmp[5];
	char nb_del_tmp[6];
	char nb_rat_tmp[6];

	src_icrf_record temp_record;
	SRC_ICRF_Data.clear();

	fs.open(filename_ICRF, fstream::in);

	if (fs.is_open())
	{
		while (!fs.eof())
		{
			fs.getline(data, 512); // telescope name and id
			if (strncmp(data, "ICRF J", 6) == 0)					/* search for the first line of data */
			{
				strncpy(temp_record.src_name_j2000, data + 5, 16);	/* source name J2000 */
				strncpy(temp_record.src_name_1950, data + 25, 8);	/* source name B1950 */
				temp_record.inf = data[35];
				
				strncpy(ra_h_tmp, data + 41, 2);					/* right ascension hh */
				sscanf(ra_h_tmp, "%lf", &temp_record.ra_h);
				strncpy(ra_m_tmp, data + 43, 2);					/* right ascension mm */
				sscanf(ra_m_tmp, "%lf", &temp_record.ra_m);
				strncpy(ra_s_tmp, data + 45, 11);					/* right ascension ss */
				sscanf(ra_s_tmp, "%lf", &temp_record.ra_s);

				strncpy(dec_d_tmp, data + 61, 3);					/* declination degrees */
				sscanf(dec_d_tmp, "%lf", &temp_record.dec_d);
				strncpy(dec_m_tmp, data + 65, 2);					/* declination mm */
				sscanf(dec_m_tmp, "%lf", &temp_record.dec_m);
				strncpy(dec_s_tmp, data + 67, 10);					/* declination ss */
				sscanf(dec_s_tmp, "%lf", &temp_record.dec_s);

				strncpy(ra_err_tmp, data + 83, 10);					/* right ascension err ss */
				sscanf(ra_err_tmp, "%lf", &temp_record.ra_err);

				strncpy(dec_err_tmp, data + 98, 9);					/* declination err ss */
				sscanf(dec_err_tmp, "%lf", &temp_record.dec_err);

				strncpy(ra_dec_corr_tmp, data + 108, 7);			/* RA-DEC correction */
				sscanf(ra_dec_corr_tmp, "%lf", &temp_record.ra_dec_corr);
				
				strncpy(mjd_mean_tmp, data + 118, 7);			/* MJD mean */
				sscanf(mjd_mean_tmp, "%lf", &temp_record.MJD_mean);

				strncpy(mjd_first_tmp, data + 127, 7);			/* MJD first */
				sscanf(mjd_first_tmp, "%lf", &temp_record.MJD_first);

				strncpy(mjd_last_tmp, data + 136, 7);			/* MJD last */
				sscanf(mjd_last_tmp, "%lf", &temp_record.MJD_last);

				strncpy(nb_sess_tmp, data + 144, 5);			/* Nb sessions */
				sscanf(nb_sess_tmp, "%lf", &temp_record.Nb_sess);

				strncpy(nb_del_tmp, data + 150, 6);			/* Nb del */
				sscanf(nb_del_tmp, "%lf", &temp_record.Nb_del);

				strncpy(nb_rat_tmp, data + 157, 6);			/* Nb rat */
				sscanf(nb_rat_tmp, "%lf", &temp_record.Nb_rat);

				SRC_ICRF_Data.push_back(temp_record);
			}
		}
		return 1;
	}
	else return -1;

	return 0;
}

/* --------------------------------------------------------- */
/* Read telescopes coordinates and velocities		VTRF.cat */
/* Format:													 */
/* --------------------------------------------------------- */
int ReadVTRF(vector<ant_vtrf_record> &ANT_VTRF_Data, char filename_VTRF[256])
{
	fstream fs;
	char data[1024];
	char temp_data[1024];

	char tmp_domes_nb[9];

	ant_vtrf_record temp_record;
	ANT_VTRF_Data.clear();

	fs.open(filename_VTRF, fstream::in);

	if (fs.is_open())
	{
		while (!fs.eof())
		{
			fs.getline(data, 1024); // telescope name and id
			if (strncmp(data, "#", 1) != 0)						/* skip commented lines */
			{
				strncpy(temp_record.domes_nb, data, 9);			/* domes nb telescope ID */
				strncpy(temp_record.ant_name, data + 10, 9);	/* antenna name */
				strncpy(temp_record.tech_id, data + 32, 4);		/* tech id */

				strncpy(temp_data, data + 36, sizeof(data));
				sscanf(temp_data, "%lf %lf %lf %lf %lf %lf", &temp_record.x, &temp_record.y, &temp_record.z, &temp_record.sx, &temp_record.sy, &temp_record.sz);
				fs.getline(data, 1024);
				if (strncmp(data, "#", 1) != 0)						/* skip commented lines */
				{
					strncpy(tmp_domes_nb, data, 9);					/* next line domes nb telescope ID to compare*/
					if (strcmp(temp_record.domes_nb, tmp_domes_nb))
					{
						strncpy(temp_data, data + 9, sizeof(data));
						sscanf(temp_data, "%lf %lf %lf %lf %lf %lf", &temp_record.vx, &temp_record.vy, &temp_record.vz, &temp_record.svx, &temp_record.svy, &temp_record.svz);
						ANT_VTRF_Data.push_back(temp_record);
					}
				}
			}
		}
		return 1;
	}
	else return -1;

	return 0;
}

/* --------------------------------------------------------- */
/* Read telescopes coordinates and velocities		ITRF.cat */
/* Format:													 */
/* --------------------------------------------------------- */
int ReadITRF(vector<ant_itrf_record> &ANT_ITRF_Data, char filename_ITRF[256])
{
	fstream fs;
	char data[1024];
	char temp_data[1024];

	ant_itrf_record temp_record;
	ANT_ITRF_Data.clear();

	fs.open(filename_ITRF, fstream::in);

	if (fs.is_open())
	{
		while (!fs.eof())
		{
			fs.getline(data, 1024); // telescope name and id
			
			sscanf(data, "%s %s %*s %s %*d %*s %s %lf", temp_record.domes_nb, temp_record.ant_name, temp_record.tech_id, temp_record.mount_type, &temp_record.axis_offset);
			ANT_ITRF_Data.push_back(temp_record);
		}
		return 1;
	}
	else return -1;

	return 0;
}