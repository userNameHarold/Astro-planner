#ifndef glambert_h
#define glambert_h


#include <iostream>
#include "tlamb.hpp" // function tlamb
#include "vlamb.hpp" // function vlamb
#include "xlamb.hpp" // function xlamb
#include "d8rt.hpp"  // function d8rt
#include "dvec3Mag.hpp" // for calculating double precision magnitude of vectors, type glm::dvec3
#include "makeUnitDvec3.hpp" // for turning vectors of type glm::dvec3 into unit vectors


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
 
 /* Note: This project relies heavily on the OpenGL Mathmematics (GLM) Library. Out of respect for thier
  * Happy Bunny license, we ask that this software not be used for military purposes.
  *
  * OpenGL Mathematics (GLM) is a header only C++ mathematics library for graphics software based on the 
  * OpenGL Shading Language (GLSL) specifications. Please consider supporting their project: https://glm.g-truc.net/
  */

#include <glm/ext/vector_double3.hpp> // for dvec3 type, 3-component vectors of double precision
#include <glm/geometric.hpp> // dot, cross, etc..
#include <glm/trigonometric.hpp> // radians, cos, etc..
 
 int glambert(glm::dvec3 *vi, glm::dvec3 *vf, double cbmu, glm::dvec3 pv1, glm::dvec3 Vv1, glm::dvec3 pv2, glm::dvec3 Vv2, double TOF, int nrev) {
	
	using namespace glm;
	
	double vr11, vr12, vr21, vr22, vt11, vt12, vt21, vt22 = 0;
	int n = 0;

	dvec3 sv1_1_3 = {pv1[0], Vv1[0], pv1[1]}; // these declarations courtesy of "trust me bro" mathmatics 
	dvec3 sv2_1_3 = {pv2[0], Vv2[0], pv2[1]}; // seriously, I have no idea why these are so, but its necessary
	dvec3 sv1_4_6 = {Vv1[1], pv1[2], Vv1[2]}; // or you wind up taking the cross product of a vector and itself 
											 // and the floating pooint handlers don't seem to catch it.
	double r1mag = dvec3Mag(sv1_1_3);
	double r2mag = dvec3Mag(sv2_1_3);
	

	dvec3 ur1xv1 = cross(sv1_1_3, sv1_4_6);
	
	ur1xv1 =  makeUnitDvec3(ur1xv1);
	
	dvec3 ux1 = makeUnitDvec3(sv1_1_3);
	dvec3 ux2 = makeUnitDvec3(sv2_1_3);

	
	dvec3 uz1 = cross(ux1, ux2);
	uz1 = makeUnitDvec3(uz1);
	

	double theta = dot(ux1, ux2);

	if (theta > 1.0){
	 theta = 1.0;
	} else if (theta < -1.0) {
	 theta = -1.0;
	}

	theta = acos(theta);
	
	double angle_to_on = dot(ur1xv1, uz1);
	
	
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
	
	
	dvec3 uy1 = cross(uz1, ux1);
	uy1 = makeUnitDvec3(uy1);
	

	dvec3 uy2 = cross(uz1, ux2);
	uy2 = makeUnitDvec3(uy2);
	
	theta = theta + 2.0 * M_PI * abs(nrev);
	
	vlamb(&vr11, &vr12, &vr21, &vr22, &vt11, &vt12, &vt21, &vt22, &n, cbmu, r1mag, r2mag, theta, TOF);
	
	
	if (abs(nrev) > 0){
		if (n == -1) {
			std::cout<<"No tminimum "<< std::endl;;
			for (int i = 0; i < 3; ++i){
				(*vi)[i] = 0.0;
				(*vf)[i] = 0.0;
			}
			return 1;
		} else if (n == 0){
			std::cout<<"No solution time"<< std::endl;;
			for (int i = 0; i < 3; ++i){
				(*vi)[i] = 0.0;
				(*vf)[i] = 0.0;
			}
			return 1;
		}
	}
	
	if((nrev > 0) && (n > 1)){
		for (int i = 0; i < 3; ++i){
			(*vi)[i] = ux1[i] * vr21 + uy1[i] * vt21;
			(*vf)[i] = ux2[i] * vr22 + uy2[i] * vt22;
		}
		return 1;
	} else{
		for(int i = 0; i < 3; ++i){
			(*vi)[i] = ux1[i] * vr11 + uy1[i] * vt11;
			(*vf)[i] = ux2[i] * vr12 + uy2[i] * vt12;
		}
		return 1;
	}
	
	
	return -1;
 }

	


#endif


