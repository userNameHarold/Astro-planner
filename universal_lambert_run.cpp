#include <cmath>
#include <iostream>
#include "magnitude_angle_universal_lambert_solver.cpp"
#include <fstream>
#include <chrono>
//#include "Gnuplot.h"

using namespace std;

int main() {
    // Define constant variables
    const double TOF = 146 * 86400;   // time of flight in seconds
    const double r1 = 149597800;      // radius of Earth Orbit in meters
    const double r2 = 108204088;      // radius of Venus Orbit in meters
    const double daily_change = ((2* M_PI) / 365);  // daily change of earth's position in radians
    
    // Array to store C3 values
    double c3[365];
	int days[365];
    
	// timers
	std::chrono::high_resolution_clock::time_point start;
	std::chrono::high_resolution_clock::time_point end;
	
	//start timer
	start = std::chrono::high_resolution_clock::now();
	
    // Loop through all the days from -15 to 15
    for (int i = -182; i <= 182; ++i) {
        double A;
        
        // special case for the hohmann transfer
        if (i == 0) {
            A = magnitude_angle_universal_lambert_solver(TOF, r1, 0, r2);
        }
        
        else {
            A = magnitude_angle_universal_lambert_solver(TOF, r1, i * daily_change, r2);
        }
        
        // Store the C3 value in the array
        c3[i + 183] = pow(A, 2);
		//cout << c3[ i + 15]<<endl;
		days[ i + 183] = i;
    }
	//stop timer
	end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
	auto milli = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
	auto average = duration / 365;
	
    std::cout << "Time taken to solve for 365 days: " <<  milli << " milliseconds" << std::endl;
	cout<<"Average time to solve: "<<average<< " microseconds" <<endl;
    
	/*
	Gnuplot plot;
	plot.set_style("lines");
	plot.set_xlabel("Days from optimal Hohmann Alignment");
	plot.set_ylabel("C_{3} (km^2/s^2)");
	plot.set_title("C_{3} needed as a function of days from hohmann alignment");
	plot.plot_xy(days, c3);
    // Plot the C3 values
    // ...
    return 0;
	*/
	return 0;
}
