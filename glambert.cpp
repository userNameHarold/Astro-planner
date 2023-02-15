#ifndef glambert_h
#define glambert_h

#include <cmath>
#include <vector>
#include <iostream>
#include "tlamb.cpp"
#include "vlamb.cpp"
#include "xlamb.cpp"
#include "d8rt.cpp"
#include "dot.cpp"
#include "cross.cpp"
#include "vectorMag.cpp"
#include "vectorScalarDivision.cpp"
#include "vectorScalarMultiplication.cpp"

/* function glambert derived from glambert.m, written by Dr. Dario Izzo of the European Space Agency (ESA) 
 * Advanced Concepts Team (ACT)
 *
 * References:
 *
 * R. H. Gooding, Technical Report 88027
 * On the Solution of Lambert's Orbital Boundary-Value Problem,
 * Royal Aerospace Establishment, April 1988
 *
 * R. H. Gooding, Technical Memo SPACE 378
 * A Procedure for the Solution of Lambert's Orbital Boundary-Value Problem
 * Royal Aerospace Establishment, March 1991
 *
 * Translated from MATLAB by R. Harold with permission
 */
 
 using namespace std;
 
 void glambert(vector<double> *vi, vector<double> *vf, double cbmu, vector<double> pv1, vector<double> Vv1, vector<double> pv2, vector<double> Vv2, double TOF, int nrev) {
	
	double vr11, vr12, vr21, vr22, vt11, vt12, vt21, vt22;
	int n;

	double r1mag = vectorMag(pv1);
	double r2mag = vectorMag(pv2);

	vector<double> ur1xv1(3);
	ur1xv1	= cross(&pv1, &Vv1);

	ur1xv1 = vectorScalarDivision(ur1xv1, vectorMag(ur1xv1));

	vector<double> ux1(3);
	ux1 = vectorScalarDivision(pv1, r1mag);
	vector<double> ux2(3);
	ux2 = vectorScalarDivision(pv2, r2mag);

	vector<double> uz1(3);
	uz1 = cross(&ux1, &ux2);
	uz1 = vectorScalarDivision(uz1, vectorMag(uz1));
	
	// vector<double> uz2; // note: in MATLAB code, uz2 is created and initialized to uz1, which is never changed after that. if incorrect values, look here

	double theta = dot(&ux1, &ux2);

	if (theta > 1.0){
	 theta = 1.0;
	} else if (theta < -1.0) {
	 theta = -1.0;
	}

	theta = acos(theta);
	
	double angle_to_on = dot(&ur1xv1, &uz1);
	if (angle_to_on > 1.0) {
		angle_to_on = 1.0;
	} else if (angle_to_on < -1.0) {
		angle_to_on = -1.0;
	}
	angle_to_on = acos(angle_to_on);
	
	if ((angle_to_on > 0.5 * M_PI) && (TOF > 0.0)){
		theta = 2.0 * M_PI - theta;
		uz1[0] = -1 * uz1[0];
		uz1[1] = -1 * uz1[1];
		uz1[2] = -1 * uz1[2];
	} else if ((angle_to_on < 0.5 * M_PI) && (TOF < 0.0)){
		theta = 2.0 * M_PI - theta;
		uz1[0] = -1 * uz1[0];
		uz1[1] = -1 * uz1[1];
		uz1[2] = -1 * uz1[2];
	}
	
	
	//uz2 = uz1; in MATLAB code, uz1 is never changed after here: if incorrect values, assign uz2 to uz1 and edit per MATLAB code
	vector<double> uy1(3);
	uy1 = cross(&uz1, &ux1);
	uy1 = vectorScalarDivision(uy1, vectorMag(uy1));
	
	vector<double> uy2(3);
	uy2 = cross(&uz1, &ux2);
	uy2 = vectorScalarDivision(uy2, vectorMag(uy2));
	
	theta = theta + 2.0 * M_PI * abs(nrev);
	
	vlamb(&vr11, &vr12, &vr21, &vr22, &vt11, &vt12, &vt21, &vt22, &n, cbmu, r1mag, r2mag, theta, TOF);
	
	if (abs(nrev) > 0){
		if (n == -1) {
			std::cout<<"No tminimum "<< std::endl;;
			for (int i = 0; i < 3; ++i){
				(*vi)[i] = 0.0;
				(*vf)[i] = 0.0;
			}
			return;
		} else if (n == 0){
			std::cout<<"No solution time"<< std::endl;;
			for (int i = 0; i < 3; ++i){
				(*vi)[i] = 0.0;
				(*vf)[i] = 0.0;
			}
			return;
		}
	}
	
	if((nrev > 0) && (n > 1)){
		for (int i = 0; i < 3; ++i){
			(*vi)[i] = ux1[i] * vr21 + uy1[i] * vt21;
			(*vf)[i] = ux2[i] * vr22 + uy2[i] * vt22;
		}
		return;
	} else{
		for(int i = 0; i < 3; ++i){
			(*vi)[i] = ux1[i] * vr11 + uy1[i] * vt11;
			(*vf)[i] = ux2[i] * vr12 + uy2[i] * vt12;
		}
		return;
	}
 }

	


#endif


