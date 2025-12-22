// nwmf2.cpp
#include "functions.h"

namespace ariadna {

void nwmf2(double latitude, double elev, Eigen::Vector2d& wmf) {
    double lat_deg = latitude / cnst::CDEGRAD;
    double l = std::abs(lat_deg);

//    printf("NWMF LAT = %.16f\n", latitude);
//    printf("NWMF ELEV = %.16f\n", elev);

    // Interpolate coefficients
    double a, b, c, dl;
    if (l <= cnst::NHMF_LATITUDES[0]) {
        a = cnst::NWMF_ABC_W2P0[0][0];
        b = cnst::NWMF_ABC_W2P0[0][1];
        c = cnst::NWMF_ABC_W2P0[0][2];
    } else if (l >= cnst::NHMF_LATITUDES[4]) {
        a = cnst::NWMF_ABC_W2P0[4][0];
        b = cnst::NWMF_ABC_W2P0[4][1];
        c = cnst::NWMF_ABC_W2P0[4][2];
    } else {
        for (int i = 0; i < 4; ++i) {
            if (l > cnst::NHMF_LATITUDES[i] && l <= cnst::NHMF_LATITUDES[i + 1]) {
                dl = (l - cnst::NHMF_LATITUDES[i]) / (cnst::NHMF_LATITUDES[i + 1] - cnst::NHMF_LATITUDES[i]);
                a = cnst::NWMF_ABC_W2P0[i][0] + dl * (cnst::NWMF_ABC_W2P0[i + 1][0] - cnst::NWMF_ABC_W2P0[i][0]);
                b = cnst::NWMF_ABC_W2P0[i][1] + dl * (cnst::NWMF_ABC_W2P0[i + 1][1] - cnst::NWMF_ABC_W2P0[i][1]);
                c = cnst::NWMF_ABC_W2P0[i][2] + dl * (cnst::NWMF_ABC_W2P0[i + 1][2] - cnst::NWMF_ABC_W2P0[i][2]);
                break;
            }
        }
    }

//    printf("NWMF a = %.16f\n", a);
//    printf("NWMF b = %.16f\n", b);
//    printf("NWMF c = %.16f\n", c);
//    printf("NWMF dl = %.16f\n", dl);    

    // Compute mapping function
    double sine = std::sin(elev);
    double cose = std::cos(elev);
    double beta = b / (sine + c);
    double gamma = a / (sine + beta);
    double topcon = 1.0 + a / (1.0 + b / (1.0 + c));
    wmf(0) = topcon / (sine + gamma);
    wmf(1) = -topcon / ((sine + gamma) * (sine + gamma)) * (cose - a / ((sine + beta) * (sine + beta)) * cose * (1.0 - b / ((sine + c) * (sine + c))));

//    printf("NWMF SINE = %.16le\n", sine);
//    printf("NWMF COSE = %.16le\n", cose);
//    printf("NWMF BETA = %.16le\n", beta);
//    printf("NWMF GAMMA = %.16le\n", gamma);
//    printf("NWMF TOPCON = %.16le\n", topcon);

//    printf("wmf(0) = %.16le\n", wmf(0));
//    printf("wmf(1) = %.16le\n\n", wmf(1));
}

} // namespace ariadna