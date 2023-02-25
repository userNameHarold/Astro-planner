#ifndef tlamb_h
#define tlamb_h

#include <glm/trigonometric.hpp> // radians, cos, etc..
#include <iostream>

/* function tlamb derived from glambert.m, written by Dr. Dario Izzo of the European Space Agency (ESA) 
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

void tlamb(double *dt, double *d2t, double *d3t, double *t, int m, double q, double qsqfm1, double x, int n){
	bool lm1 = (n == -1);
	bool l1 = (n >= 1);
	bool l2 = (n >= 2);
	bool l3 = (n == 3);
	double sw = 0.4;
	double qsq = q * q;
	double xsq = x * x;
	double u = (1.0 - x) * (1.0 + x);
	double y, z, qx, a, b, aa, bb, g, f, fg1, fg1sq, term, told, twoi1 = 0;
	double u0i, u1i, u2i, u3i, tq, tqsum, ttmold, tterm, tqterm, qz, qz2 = 0;
	int icounter = 0;
	int p = 0;
	
	if (!lm1){
		*dt = 0.0;
		*d2t = 0.0;
		*d3t = 0.0;
	}
	 
	cout<<"lm1 "<<lm1<<endl<<"m "<<m<<endl<<"x "<<x<<endl<<"u "<<u<<endl<<"sw "<<sw<<endl;
	 
	if (lm1 || (m > 0) || (x < 0.0) || (fabs(u) > sw)){
		// direct computation, no series
		
		y = sqrt(fabs(u));
		z = sqrt(qsqfm1 + qsq * xsq);
		qx = q * x;
		
		if (qx <= 0.0){
			a = z - qx;
			b = q * z - x;
		}
		
		if ((qx < 0 ) && lm1){
			aa = qsqfm1 / a;
			bb = qsqfm1 * (qsq * u - xsq) / b;
		}
		
		if (((qx == 0.0) && lm1 ) || (qx > 0.0)){
			aa = z + qx;
			bb = q * z + x;
		}
		
		if (qx > 0){
			a = qsqfm1 / aa;
			b = qsqfm1 * (qsq * u - xsq) / bb;
		}
		
		if (lm1){
			*dt = b;
			*d2t = bb;
			*d3t = aa;
		} else {
			if ((qx * u) >= 0.0){
				g = x * z + q * u;
			} else {
				g = (xsq  - qsq * u) / (x * z - q * u);
			}
			
			f = a * y;
			
			if (x <= 1){
				*t = m * M_PI + atan2(f,g);
			} else {
				if ( f > sw) {
					*t = log(f + g);
				} else {
					fg1 = f / (g + 1.0);
					term = 2.0 * fg1;
					fg1sq = fg1 * fg1;
					*t = term;
					twoi1 = 1.0;
					told = 0.0;
					
					while(*t != told){
						twoi1 += 2.0;
						term = term * fg1sq;
						told = *t;
						*t = *t + term / twoi1;
					}
				}
			}
			
			*t = 2.0 * (*t / y + b) / u;
			
			if (l1 && (z != 0.0)){
				qz = q / z;
				qz2 = qz * qz;
				qz = qz * qz2;
				*dt = (3.0 * x * *t - 4.0 * (a + qx * qsqfm1) / z) / u;
				
				if (l2){
					*d2t = ( 3.0 * *t + 5.0 * x * *dt + 4.0 * qz * qsqfm1) / u;
				}
				if (l3) {
					*d3t = (8.0 * *dt + 7.0 * x * *d2t - 12.0 * qz * qz2 * x * qsqfm1) / u;
				}
			}
		}
	} else {
		// compute by series
		u0i = 1.0;
		
		if (l1){
			u1i = 1.0;
		}
		if (l2){
			u2i = 1.0;
		}
		if (l3){
			u3i = 1.0;
		}
		
		term = 4.0;
		tq = q * qsqfm1;
		
		if ( q < 0.5){
			tqsum = 1.0 - q * qsq;
		} else {
			tqsum = (1.0 / (1.0 + q) + q) * qsqfm1;
		}
		
		ttmold = term / 3.0;
		*t = ttmold * tqsum;
		
		while((icounter < n ) || (*t != told)){
			++icounter;
			p = icounter;
			u0i = u0i * u;
			
			if(l1 && (icounter > 1)){
				u1i = u1i * u;
			}
			if(l2 && (icounter > 2)){
				u2i = u2i * u;
			}
			if (l3 && (icounter > 3)){
				u3i = u3i * u;
			}
			
			term = term * (p - 0.5) / p;
			tq = tq * qsq;
			tqsum = tqsum + tq;
			told = *t;
			tterm = term / (2.0 * p + 3.0);
			tqterm = tterm * tqsum;
			*t = *t - u0i * ((1.5 * p + 0.25) * tqterm / (p * p - 0.25) - ttmold * tq);
			ttmold = tterm;
			tqterm = tqterm * p;
			
			if (l1) {
				*dt = *dt + tqterm * u1i;
			}
			if(l2) {
				*d2t = *d2t + tqterm * u2i * (p - 1.0);
			}
			if(l3) {
				*d3t = *d3t + tqterm * u3i * (p - 1.0) * ( p - 2.0);
			}
		}
		
		if(l3) {
			*d3t = 8.0 * x * (1.5 * *d2t - xsq * *d3t);
		}
		if (l2) {
			*d2t = 2.0 * (2.0 * xsq * *d2t - *dt);
		}
		if(l1) {
			*dt = -2.0 * x * *dt;
		}
		
		*t = *t / xsq;
	}
	
	return;
}


#endif