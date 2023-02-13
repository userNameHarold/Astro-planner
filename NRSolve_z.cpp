#ifndef NRSolve_z_h
#define NRSolve_z_h
#include <iostream>
#include <cmath>
#include <chrono>
#include "UniversalLambert_f.cpp"
#include "UniversalLambert_dF.cpp"


double NRSolve_z(double r1, double r2, double A, double TOF, double z0 = 1.0, double absTol = pow(10,-10), double relTol = pow(10,-10), double mu = 1.327 * pow(10,11), int max_recursion_depth = 500, bool use_timer = false){
	
	/* r1 is the magitude of the vector that points from the barycenter of the system to the departure body
	   r2 is the magnitude of the vector that points from the barycenter of the system to the arrival body
   	   A is a constant calculated for the universal variable solution to Lambert's problem
	   TOF is the time of flight between the 2 bodies, in seconds
 	   z0 is the initial guess of the system
	   the absolute and relative tollerances define how much accuracy we need for the calculation of Z
	   mu is the gravitational paramter of the body both the arrival and departure point are orbiting. Default value is mu_sol
	   max_recursion_depth defines the number of iterations to try to converge before aborting. Z is garaunteed to converge eventually, but could use exceptional resources is bad input values are given
	   use_timer is a debug feature to measue execution time, defaults to off */

 
    double z[max_recursion_depth];
    int i = 1;
	double F;
	double dF;
	double relTolVal;
    double absTolVal;
	std::chrono::high_resolution_clock::time_point start;
	std::chrono::high_resolution_clock::time_point end;
    if (use_timer) {

	 start = std::chrono::high_resolution_clock::now();
    }

    z[0] = z0;
    while (i < max_recursion_depth) {
        F = UniversalLambert_F(z[i-1], r1, r2, A, TOF, mu);
        dF = UniversalLambert_dF(z[i-1], r1, r2, A);
        z[i] = z[i-1] - (F) / (dF);

        relTolVal = std::abs((z[i] - z[i-1]) / z[i-1]);
        absTolVal = std::abs(z[i] - z[i-1]);

        if ((relTolVal < relTol) && (absTolVal < absTol)) {
            break;
        }
        ++i;
    }

    if (i == max_recursion_depth) {
        return NAN;
    }
    if (use_timer) {
        end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        std::cout << "Time taken by function: " << duration << " microseconds" << std::endl;
		std::cout << "Z converged to " << z[i] << " on iteration: " << i <<  std::endl;
    }
    return z[i];
}
#endif


