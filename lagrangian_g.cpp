#ifndef LAGRANGIAN_G_H
#define LAGRANGIAN_G_H
#include <cmath>

double lagrangian_g(double z, double r1, double r2, double A, double mu) {
  return A * sqrt(UniversalLambertsProblem_y(z, r1, r2, A) / mu);
}

#endif

