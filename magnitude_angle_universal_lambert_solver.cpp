#ifndef magnitude_angle_universal_lambert_solver_h
#define magnitude_angle_universal_lambert_solver_h
#include <cmath>
#include <iostream>
#include "UniversalLambertsProblem_y.cpp"
#include "stumpffc2.cpp"
#include "stumpffc3.cpp"
#include "lagrangian_f.cpp"
#include "lagrangian_g.cpp"
#include "lagrangian_gdot.cpp"
#include "NRSolve_z.cpp"

using namespace std;

double magnitude_angle_universal_lambert_solver(double TOF, double r1, double theta1, double r2, double theta2 = M_PI, double mu=1.327*pow(10,11)) {
    double dTheta = theta2 - theta1;
    
    if (dTheta == M_PI) {
		//Hohmann transfer needs to be treated differently
        double v1num = abs(sqrt(mu/r1) * (sqrt((2*r2)/(r1+r2)) - 1));
        //double v2num = abs(sqrt(mu/r2) * (1 - sqrt((2*r1)/(r1+r2))));
        return  v1num; // + v2num; // we just want the initial velocity for this model, left in for reference 
    } else {
        double r1vec[2] = {r1*cos(theta1), r1*sin(theta1)};
        double r2vec[2] = {r2*cos(theta2), r2*sin(theta2)};
        
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
        
        double v1num[2] = {(1/g) * (r2vec[0] - f * r1vec[0]), (1/g) * (r2vec[1] - f * r1vec[1])};
        //double v2num[2] = {(1/g) * (gdot * r2vec[0] - r1vec[0]), (1/g) * (gdot * r2vec[1] - r1vec[1])};
		// left in for reference, not needed for current usage
        
        double v = sqrt(pow(v1num[0],2) + pow(v1num[1],2)); // turns vector into scalar
        //v2 = sqrt(pow(v2num[0],2) + pow(v2num[1],2)); // turns vector into scalar, not needed for current usage of function
		
		double vinit[2] = {sqrt(mu/r1) * (-sin(theta1)), sqrt(mu/r1) * cos(theta1)};
		
		
		double deltaV[2] ={ v1num[0] - vinit[0], v1num[1] - vinit[1]};
		double returnDeltaV = sqrt(pow(deltaV[0],2) + pow(deltaV[1],2));
		return returnDeltaV;
		
    }
}



#endif

