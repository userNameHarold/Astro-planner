#ifndef UniversalLambertsProblem_y_h
#define UniversalLambertsProblem_y_h
#include <cmath>
#include "stumpffc2.cpp"
#include "stumpffc3.cpp"

double UniversalLambertsProblem_y(double z, double r1, double r2, double A)
{
    double c2 = stumpffc2(z);
    double c3 = stumpffc3(z);
    double yz = r1 + r2 + A * ((z * c3 - 1) / sqrt(c2));
    return yz;
}
#endif
