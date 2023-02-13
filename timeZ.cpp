#include <cmath>
#include <iostream>
#include "NRSolve_z.cpp"

using namespace std;

int main(){
	std::chrono::high_resolution_clock::time_point start;
	std::chrono::high_resolution_clock::time_point end;
	 const double TOF = 146 * 86400;   // time of flight in seconds
    const double r1 = 149597800;      // radius of Earth Orbit in meters
    const double r2 = 108204088;      // radius of Venus Orbit in meters
	const double daily_change = ((2* M_PI) / 365);  // daily change of earth's position in radians
	
	start = std::chrono::high_resolution_clock::now();
	for (int i = -182; i <= 182; ++i){
	double A = sin(i*daily_change)*(sqrt((r1*r2)/(1 - cos(i*daily_change))));
		//cout<<"A = " << A <<endl;
        double z = NRSolve_z(r1, r2, A, TOF);
		//cout<<"Z = " << z <<endl;
	}
	
	end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
	auto milli = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
	auto average = duration / 365;
	
    std::cout << "Time taken to solve Z for 365 days: " <<  milli << " milliseconds" << std::endl;
	cout<<"Average time to solve Z: "<<average<< " microseconds" <<endl;
	
	
	
}