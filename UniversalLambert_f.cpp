#ifndef UniversalLambert_F_h
#define UniversalLambert_F_h
#include <cmath>
#include "stumpffc2.cpp"
#include "stumpffc3.cpp"
#include "UniversalLambertsProblem_y.cpp"

double UniversalLambert_F(double z, double r1, double r2, double A, double TOF, double mu = 1.327e11) {
double yz = UniversalLambertsProblem_y(z, r1, r2, A);
double c2 = stumpffc2(z);
double c3 = stumpffc3(z);

return c3 * pow(yz/c2, 3./2) + A*sqrt(yz) - sqrt(mu) * TOF;
}

#endif