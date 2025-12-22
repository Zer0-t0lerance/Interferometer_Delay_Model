#include "../src/functions.h"

int main(int argc, char **argv) {
    double el_rad = 0.674516446350307;
    double temp_k = 30.427 + cnst::KELVIN_OFFSET;
    double humid_f = 60.973 / cnst::PERCENT_TO_FRACTION;
    double press_hg = 755.087095978288;

    double result = ariadna::sbend(el_rad, temp_k, humid_f, press_hg);

    printf("SBend result: %f\n", result);
    
    double el_rad2 = 0.693637398629522;
    double temp_k2 = 25.884 + cnst::KELVIN_OFFSET;
    double humid_f2 = 43.395 / cnst::PERCENT_TO_FRACTION;
    double press_hg2 = 645.977873180360;

    double result2 = ariadna::sbend(el_rad2, temp_k2, humid_f2, press_hg2);

    printf("SBend result: %f\n", result2);
}