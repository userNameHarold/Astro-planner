#ifndef UniversalLambert_dF_h
#define UniversalLambert_dF_h
#include <math.h>
#include "stumpffc2.cpp"
#include "stumpffc3.cpp"
#include "UniversalLambertsProblem_y.cpp"

double UniversalLambert_dF(double z, double r1, double r2, double A) {
    double yz, c2, c3;

    if (z == 0) {
        double y0 = UniversalLambertsProblem_y(0, r1, r2, A);
        return (sqrt(2)) / 40 * pow(y0, 1.5) + (A / 8) * (sqrt(y0) + A * sqrt(1 / (2 * y0)));
    } else {
        yz = UniversalLambertsProblem_y(z, r1, r2, A);
        c2 = stumpffc2(z);
        c3 = stumpffc3(z);

        return pow(yz / c2, 1.5) * ((1 / (2 * z)) * (c2 - (3 * c3 / (2 * c2))) + (3 * pow(c3, 2)) / (4 * c2)) +
               (A / 8) * ((3 * c3 / c2) * sqrt(yz) + A * sqrt(c2 / yz));
    }
}
#endif

