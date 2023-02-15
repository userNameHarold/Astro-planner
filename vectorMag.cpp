#ifndef vectorMag_h
#define vectorMag_h

#include <cmath>
#include <vector>

using namespace std;

double vectorMag(vector<double> v){
	double sum = 0;
	
	for(int i = 0; i < v.size(); ++i){
		sum += v[i] * v[i];
	}
	return sqrt(sum);

}

#endif

