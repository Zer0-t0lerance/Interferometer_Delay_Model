#include "functions.h"

namespace ariadna {

double find_val(double mjd) {
    if (mjd < 41317.0) return 0.0; // До 1972 года
    int count = 0;
    for (int i = 0; i < 28; i++) {
        if (mjd >= cnst::LEAP_SECOND_MJD_INTERVALS[i]) count++;
        else break;
    }
    return (cnst::LEAP_SECONDS_BASE + count);
}

void nsec(double mjd, double& idelt) {
    idelt = find_val(mjd);
}

} // namespace ariadna