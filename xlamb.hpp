#ifndef xlamb_h
#define xlamb_h

#include <glm/trigonometric.hpp> // radians, cos, etc..
#include "d8rt.hpp"
#include "tlamb.hpp"
#include <iostream> //for cout

/* function xlamb derived from glambert.m, written by Dr. Dario Izzo of the European Space Agency (ESA) 
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



void xlamb(double *x, double *xpl, int *n, int m, double q, double qsqfm1, double tin){
	const double pi = M_PI;
	const double tol = 3*pow(10,-7);
	const double c0 = 1.7;
	const double c1 = 0.5;
	const double c2 = 0.03;
	const double c3 = 0.15;
	const double c41 = 1.0;
	const double c42 = 0.24;
	double t0, dt, d2t, d3t, tdiff, w, xm, xmold, xtest = 0;
	double tmin, tdiffm, d2t2, tdiff0, t = 0;
	bool solnflag = false;
	bool termflag = false;

	*xpl = 0;
	*x = 0;
	*n = 0;
	
	double thr2 = atan2(qsqfm1, 2.0 * q) / pi;
	
	if (m == 0) { // single rev starter from t (@ x = 0)
		*n = 1;
		tlamb(&dt, &d2t, &d3t, &t0, m, q, qsqfm1, 0.0, 0); // keep an eye on this line for errors
		// std::cout<<dt<<d2t<<d3t<<t0<<std::endl;
		tdiff = tin - t0;
		// std::cout<<tin<<" "<<t0<<" "<<tdiff<<std::endl;
		if (tdiff <= 0.0) {
		  *x = t0 * tdiff / (-4.0 * tin);
		} else {
		  *x = -tdiff / (tdiff + 4.0);
		  w = *x + c0 * sqrt(2.0 * (1.0 - thr2));
		  if (w < 0.0) {
			*x = *x - sqrt(d8rt(-w)) * (*x + sqrt(tdiff / (tdiff + 1.5 * t0)));
		  }
		  w = 4.0 / (4.0 + tdiff);
		  *x = *x * (1.0 + *x * (c1 * w - c2 * *x * sqrt(w)));
		}
		// std::cout<<w<<*x<<std::endl;
	} else { // multirevs, get t(min) as basis 
		xm = 1.0 / (1.5 * (m + 0.5) * pi);
		
		if (thr2 < 0.5){
			xm = d8rt(2.0 * thr2) * xm;
		} else if (thr2 > 0.5) {
			xm = (2.0 - d8rt(2.0 - 2.0 * thr2)) * xm;
		} else {
			std::cout<< "glambert does not have code for this eventuallity, thr2 == 0.5 '('xlamb.cpp, line 50')'" << std::endl;
		}
		

		for(int i = 0; i < 12; ++i){
			
			tlamb(&dt, &d2t, &d3t, &tmin, m, q, qsqfm1, xm, 3); // keep an eye on this line for errors
			
			if (d2t == 0.0){
				solnflag = true;
				break;
			}
			
			xmold = xm;
			xm = xm - dt * d2t / (d2t * d2t - dt * d3t / 2.0);
			xtest = fabs(xmold / xm - 1.0);
			
			if (xtest <= tol){
				solnflag = true;
				break;
			}
		}
		
		if (solnflag){
			tdiffm = tin - tmin;
			
			if (tdiffm < 0.0){
				*n = 0;
				termflag = true;
				// this is the case where the provided m-value produces no solution
			} else if (tdiffm == 0.0){
				*x = xm;
				*n = 1;
				termflag = true;
				// this is the case where unique solution has already been found from x(tmin)
			} else {
				*n = 3;
				
				if (d2t == 0.0){
					d2t = 6.0 * m * pi;
				}
				
				*x = sqrt(tdiffm / (d2t / 2.0 + tdiffm / (pow((1.0 - xm), 2))));
				w = xm + *x;
				w = w * 4.0 / (4.0 + tdiffm) + pow((1.0 - w), 2);
				*x = *x * (1.0 - (1.0 + m + c41 * (thr2 - 0.5)) / (1.0 + c3 * m) * *x * (c1 * w + c2 * *x * sqrt(w))) + xm;
				d2t2 = d2t / 2.0;
				
				if (*x >= 1.0){
					*n = 1;
					
					tlamb(&dt, &d2t, &d3t, &t0, m, q, qsqfm1, 0.0, 0); // keep an eye on this line for errors
					tdiff0 = t0 - tmin;
					tdiff = tin - t0;
					
					if (tdiff <= 0.0){
						*x = xm - sqrt(tdiffm / (d2t2 - tdiffm * (d2t2 / tdiff0 - 1.0 / (xm * xm))));
					} else{
						*x = -tdiff / (tdiff + 4.0);
						w = *x + c0 * sqrt(2.0 * (1.0 - thr2));
						
						if (w < 0.0){
							*x = *x - sqrt(d8rt(-w)) * (*x + sqrt(tdiff / (tdiff + 1.5 * t0)));
						}
						
						w = 4.0 / (4.0 + tdiff);
						*x = *x * (1.0 + (1.0 + m + c42 * (thr2 - 0.5)) / (1.0 + c3 * m) * *x * (c1 * w - c2 * *x * sqrt(w)));
						
						if (*x <= -1.0){
							*n = *n - 1;
							
							if (*n == 1){
								*x = *xpl;
							}
						}
					}
				}
			}
		} else {
			*n = -1;
			termflag = true;
		}
	}
	std::cout<<"termflag "<<termflag<<std::endl;
	while( !termflag){
		for(int i = 0; i < 3; ++i){
			tlamb(&dt, &d2t, &d3t, &t, m, q, qsqfm1, *x, 2); // keep an eye on this line for errors
			
			using namespace std; //delete when done debugging
			
			// cout<<"dt "<<dt<<endl;
			// cout<<"d2t "<<d2t<<endl;
			// cout<<"d3t "<<d3t<<endl;
			// cout<<"m "<<m<<endl;
			// cout<<"q "<<q<<endl;
			// cout<<"qsqfm1 "<<qsqfm1<<endl;
			// cout<<"t pre-change "<<t<<endl;
			
			t = tin - t;
			
			// cout<<"t post change "<<t<<endl;
			// cout<<"x before 'if' "<<*x<<endl;
			if (dt != 0.0){
				*x = *x + t * dt / (dt * dt + t * d2t / 2.0);
				// cout<<"x after 'if' "<<*x<<endl;
			}
		}
		
		if ( *n != 3){
			//exit when only one solution, normally when m = 0
			return;
		}
		
		*n = 2;
		*xpl = *x;
		
		// second multi-rev starter
		
		tlamb(&dt, &d2t, &d3t, &t0, m, q, qsqfm1, 0.0, 0); // keep an eye on this line for errors
		
		tdiff0 = t0 - tmin;
		tdiff = tin - t0;
		
		if (tdiff <= 0.0){
			*x = xm - sqrt(tdiffm / (d2t2 - tdiffm * (d2t2 / tdiff0 - 1.0 / (xm * xm))));
		} else {
			*x = -tdiff / (tdiff + 4.0);
			w = *x + c0 * sqrt(2.0 * (1.0 - thr2));
			if(w < 0.0){
				*x = *x - sqrt(d8rt(-w)) * (*x + sqrt(tdiff / (tdiff + 1.5 * t0)));
			}
			w = 4.0 / (4.0 + tdiff);
			*x = *x * (1.0 + (1.0 + m + c42 * (thr2 - 0.5)) / (1.0 + c3 * m) * *x * (c1 * w - c2 * *x * sqrt(w)));
			
			if (*x <= -1.0){
				*n = *n - 1;
				if (*n == 1){
					// no finite soln for *x < xm
					*x = *xpl;
				}
			}
		}
	}
	
	return;
}
	
#endif
