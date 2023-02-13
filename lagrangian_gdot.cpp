#ifndef LAGRANGIAN_G_DOT_H
#define LAGRANGIAN_G_DOT_H
#include "UniversalLambertsProblem_y.cpp"

double lagrangian_gdot(double z, double r1, double r2, double A)
{
    return 1 - (UniversalLambertsProblem_y(z, r1, r2, A) / r2);
}

#endif

