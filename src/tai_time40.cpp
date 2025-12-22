#include "functions.h"

namespace ariadna {

void tai_time(double mjd, double UTC, double &TAI, double &TT) {
	double idelt;

	nsec(mjd, idelt);
//	printf("idelt=%lf\n", idelt);

	TAI = UTC + double(idelt) / 86400.0;
	TT = TAI + 32.184 / 86400.0;
}
} // namespace ariadna