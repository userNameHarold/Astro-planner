#ifndef vlamb_h
#define vlamb_h

#include <glm/trigonometric.hpp> // radians, cos, etc..
#include "tlamb.hpp" // tlamb function 
#include "xlamb.hpp" //xlamb function


/* function vlamb derived from glambert.m, written by Dr. Dario Izzo of the European Space Agency (ESA) 
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





void vlamb(double *vr11, double *vr12, double *vr21, double *vr22, double *vt11, double *vt12, double *vt21, double *vt22, int *n, double cbmu, double r1, double r2, double thr2, double tdelt){
	const double pi = M_PI;
	double dr, r1r2, r1r2th, csq, c, s, gms, qsqfm1, q, rho, sig = 0;
	double x, x1, x2, qzminx, qzplx, zplqx, t, foo, vt1, vt2, vr1, vr2 = 0;
	int m  = 0;
	
	//ensure pointers for return values dont have bad values already
	*vr11 = 0;
	*vr12 = 0;
	*vr21 = 0;
	*vr22 = 0;
	*vt11 = 0;
	*vt12 = 0;
	*vt21 = 0;
	*vt22 = 0;
	*n = 0;

	// adjust for multi-rev
	while(thr2 > (2.0 * pi)){
		thr2 = thr2 - 2.0 * pi;
		++m;
	}
	
	//set up variables 
	thr2 = thr2 / 2.0;
	dr = r1 - r2;
	r1r2 = r1 * r2;
	r1r2th = 4.0 * r1r2 * pow((sin(thr2)), 2);
	csq = dr * dr + r1r2th;
	c = sqrt(csq);
	s = (r1 + r2 + c) / 2.0;
	gms = sqrt(cbmu * s / 2.0);
	qsqfm1 = c/s;
	q = sqrt(r1r2) * cos(thr2) / s;
	
	if (c != 0.0){
		rho = dr/c;
		sig = r1r2th / csq;
	} else{
		rho = 0.0;
		sig = 1.0;
	}
	
	t = 4.0 * gms * tdelt / (s * s);
	

	xlamb(&x1, &x2, n, m, q, qsqfm1, t);
	
	
	for(int i = 0; i < *n; ++i){
		if (i == 0){
			x = x1;
		} else{
			x = x2;
		}

		tlamb(&qzminx, &qzplx, &zplqx, &foo, m, q, qsqfm1, x, -1); 
		
		
		vt2 = gms * zplqx * sqrt(sig);

		vr1 = gms * (qzminx - qzplx * rho) / r1;

		vt1 = vt2 / r1;

		vr2 = -gms * (qzminx + qzplx * rho) / r2;

		vt2 = vt2 / r2;
		
		if ( i == 0){
			*vr11 = vr1;
			*vt11 = vt1;
			*vr12 = vr2;
			*vt12 = vt2;
		} else if (i == (*n - 1)){
			*vr21 = vr1;
			*vt21 = vt1;
			*vr22 = vr2;
			*vt22 = vt2;
		}
		
	}
}




#endif