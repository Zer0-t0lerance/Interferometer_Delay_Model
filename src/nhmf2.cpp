// nhmf2.cpp
#include "functions.h"

namespace ariadna {

void nhmf2(double epoch, double latitude, double height, double elev, Eigen::Vector2d& hmf) {
    // Convert height to kilometers
    double hs_km = height / cnst::NHMF_HEIGHT_SCALE;

    // Calculate day of year
    double doy = epoch - cnst::NHMF_DOY_REF;
    double lat_deg = latitude / cnst::CDEGRAD;
//    printf("LAT_DEG = %le\n", lat_deg);
    if (lat_deg < 0.0) {
        doy += cnst::NHMF_SEASON_OFFSET;
    }
    double doy_atm = doy - cnst::NHMF_PHASE_OFFSET;
    double doyr_atm = doy_atm * cnst::TWOPI / 365.25;

//    printf("NHMF DOY = %.16le\n", doy);
//    printf("NHMF DOY_ATM = %.16le\n", doy_atm);
//    printf("NHMF DOYR_ATM = %.16le\n", doyr_atm);

//    printf("NHMF EPOCH = %.16le\n", epoch);
//    printf("NHMF ELEV = %.16le\n", elev);

//    printf("NHMF HS_KM = %.16le\n", hs_km);

//    printf("NHMF LAT = %.16le\n", latitude);
//    printf("NHMF LAT_DEG = %.16le\n", lat_deg);

    double cost = std::cos(doyr_atm);

    // Interpolate coefficients
    double a, b, c, dl;
    double l = std::abs(lat_deg);
//    printf("ABS LAT_DEG = %le\n", l);
    if (l <= cnst::NHMF_LATITUDES[0]) {
        a = cnst::NHMF_ABC_AVG[0][0];
        b = cnst::NHMF_ABC_AVG[0][1];
        c = cnst::NHMF_ABC_AVG[0][2];
    } else if (l >= cnst::NHMF_LATITUDES[4]) {
        a = cnst::NHMF_ABC_AVG[4][0];
        b = cnst::NHMF_ABC_AVG[4][1];
        c = cnst::NHMF_ABC_AVG[4][2];
    } else {
        for (int i = 0; i < 4; ++i) {
            if (l > cnst::NHMF_LATITUDES[i] && l <= cnst::NHMF_LATITUDES[i + 1]) {
                dl = (l - cnst::NHMF_LATITUDES[i]) / (cnst::NHMF_LATITUDES[i + 1] - cnst::NHMF_LATITUDES[i]);
                double daavg = cnst::NHMF_ABC_AVG[i + 1][0] - cnst::NHMF_ABC_AVG[i][0];
                double daamp = cnst::NHMF_ABC_AMP[i + 1][0] - cnst::NHMF_ABC_AMP[i][0];
                double dbavg = cnst::NHMF_ABC_AVG[i + 1][1] - cnst::NHMF_ABC_AVG[i][1];
                double dbamp = cnst::NHMF_ABC_AMP[i + 1][1] - cnst::NHMF_ABC_AMP[i][1];
                double dcavg = cnst::NHMF_ABC_AVG[i + 1][2] - cnst::NHMF_ABC_AVG[i][2];
                double dcamp = cnst::NHMF_ABC_AMP[i + 1][2] - cnst::NHMF_ABC_AMP[i][2];
                a = cnst::NHMF_ABC_AVG[i][0] + dl * daavg - (cnst::NHMF_ABC_AMP[i][0] + dl * daamp) * cost;
                b = cnst::NHMF_ABC_AVG[i][1] + dl * dbavg - (cnst::NHMF_ABC_AMP[i][1] + dl * dbamp) * cost;
                c = cnst::NHMF_ABC_AVG[i][2] + dl * dcavg - (cnst::NHMF_ABC_AMP[i][2] + dl * dcamp) * cost;
                break;
            }
        }
    }

//    printf("NHMF a = %.16le\n", a);
//    printf("NHMF b = %.16le\n", b);
//    printf("NHMF c = %.16le\n", c);

    // Compute mapping function
    double sine = std::sin(elev);
    double cose = std::cos(elev);
    double beta = b / (sine + c);
    double gamma = a / (sine + beta);
    double topcon = 1.0 + a / (1.0 + b / (1.0 + c));

    hmf(0) = topcon / (sine + gamma);
    hmf(1) = -topcon / ((sine + gamma) * (sine + gamma)) *
     (cose - a / ((sine + beta) * (sine + beta)) * cose *
      (1.0 - b / ((sine + c) * (sine + c))));

    // Apply height correction
    beta = cnst::NHMF_B_HT / (sine + cnst::NHMF_C_HT);
    gamma = cnst::NHMF_A_HT / (sine + beta);
    topcon = (1.0 + cnst::NHMF_A_HT / (1.0 + cnst::NHMF_B_HT / (1.0 + cnst::NHMF_C_HT)));

//    printf("NHMF SINE = %.16le\n", sine);
//    printf("NHMF COSE = %.16le\n", cose);
//    printf("NHMF BETA = %.16le\n", beta);
//    printf("NHMF GAMMA = %.16le\n", gamma);
//    printf("NHMF TOPCON = %.16le\n", topcon);

    double ht_corr_coef = 1.0 / sine - topcon / (sine + gamma);
    double ht_corr = ht_corr_coef * hs_km;
    hmf(0) += ht_corr;

//  Calculate the component of d hmf/ d el due to ht_corr.

    double dhmf_ht_del = -(topcon * cose / ((sine + gamma) * (sine + gamma))) * 
                     (1.0 - cnst::NHMF_A_HT / ((sine + beta) * (sine + beta)) * 
                     (1.0 - cnst::NHMF_B_HT / ((sine + cnst::NHMF_C_HT) * (sine + cnst::NHMF_C_HT))));

    double dht_corr_coef_del = -cose/(sine*sine) - dhmf_ht_del;
    double dht_corr_del = dht_corr_coef_del * hs_km;
    hmf(1) += dht_corr_del;

//    printf("NHMF ht_corr_coef = %.16le\n", ht_corr_coef);
//    printf("NHMF ht_corr = %.16le\n", ht_corr);
//    printf("NHMF dhmf_ht_del = %.16le\n", dhmf_ht_del);
//    printf("NHMF dht_corr_coef_del = %.16le\n", dht_corr_coef_del);
//    printf("NHMF dht_corr_del = %.16le\n", dht_corr_del);
//    printf("NHMF hmf(0) = %.16le\n", hmf(0));
//    printf("NHMF hmf(1) = %.16le\n\n", hmf(1));
}

} // namespace ariadna