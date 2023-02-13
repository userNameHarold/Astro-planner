#ifndef Universal_lambert_vector_h
#define Universal_lambert_vector_h
#include <cmath>
#include <vector>
#include <iostream>
#include "UniversalLambertsProblem_y.cpp"
#include "stumpffc2.cpp"
#include "stumpffc3.cpp"
#include "lagrangian_f.cpp"
#include "lagrangian_g.cpp"
#include "lagrangian_gdot.cpp"
#include "NRSolve_z.cpp"

using namespace std;

vector<double> Universal_lambert_vector(vector<double> r1vec, vector<double> r2vec, double TOF, int returnSelect = 1, double mu = 1.327* pow(10,11)) {
	/*
	Universal_lambert_vector takes arguments of 2, 3-Dimensional position vectors, time of flight (seconds), an optional return selector, and optionally the 
	gravitational contant of the body the 2 vectors are orbiting. mu defaults to mu_sol
	return selector has 2 cases:
	returnSelect = 1
		returns the vector v1, the total velocity needed to initiate the transfer
	returnSelect = 2
		returns the vector v2, the total velocity needed to end the transfer at the target body
	
	*/
	
    double r1 = sqrt(r1vec[0] * r1vec[0] + r1vec[1] * r1vec[1] + r1vec[2] * r1vec[2]);
    double r2 = sqrt(r2vec[0] * r2vec[0] + r2vec[1] * r2vec[1] + r2vec[2] * r2vec[2]);

    double dotProduct = (r1vec[0] * r2vec[0] + r1vec[1] * r2vec[1] + r1vec[2] * r2vec[2]);
    double dTheta = acos(dotProduct / (r1 * r2));

    double typeCheck = r1 * r2 * sin(dTheta);
    if (typeCheck < 0) {
        dTheta = 2 * M_PI - dTheta;
    }

    double A = sin(dTheta)*(sqrt((r1*r2)/(1 - cos(dTheta))));
	//cout<<"A = " << A <<endl;
	double z = NRSolve_z(r1, r2, A, TOF);
	//cout<<"Z = " << z <<endl;
	
	double f = lagrangian_f(z,r1,r2,A);
	//cout<< "F = "<<f <<endl;
	double g = lagrangian_g(z,r1,r2,A,mu);
	//cout<< "G =" <<g <<endl;
	double gdot = lagrangian_gdot(z,r1,r2,A);
	//cout<< "Gdot = " <<gdot <<endl;
	

    vector<double> v1(3);
    vector<double> v2(3);
	
	if (returnSelect == 1) {
		for (int i = 0; i < 3; ++i) {
			v1[i] = (1 / g) * (r2vec[i] - f * r1vec[i]);
		}
		return v1;
	}else {
		for (int i = 0; i < 3; ++i) {
			v2[i] = (1 / g) * (gdot * r2vec[i] - r1vec[i]);
		}
		return v2;
	}


}

#endif

