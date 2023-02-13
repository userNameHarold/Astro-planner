#ifndef LAGRANGE_FDOT_H
#define LAGRANGE_FDOT_H
#include <cmath>
#include "stumpffc2.h"
#include "stumpffc3.h"
#include "UniversalLambertsProblem_y.h"

using namespace std;

double langrangian_fdot(double z, double r1, double r2, double A, double mu)
{
    double answer = (sqrt(mu)/(r1*r2)) * sqrt(UniversalLambertsProblem_y(z, r1, r2, A) / stumpffc2(z)) * (z * stumpffc3(z) - 1);
    return answer;
}


#endif

