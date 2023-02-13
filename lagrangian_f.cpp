#ifndef LAGRANGIAN_F_H
#define LAGRANGIAN_F_H
#include "UniversalLambertsProblem_y.cpp"

double lagrangian_f(double z, double r1, double r2, double A)
{
    return 1 - (UniversalLambertsProblem_y(z, r1, r2, A) / r1);
}

#endif

