#ifndef vectorScalarDivision_h
#define vectorScalarDivision_h

#include <cmath>
#include <vector>

using namespace std;

vector<double> vectorScalarDivision(vector<double> v, double num){
	int vectorSize = v.size();
	vector<double> result(vectorSize);
	
	for(int i = 0; i < vectorSize; ++i){
		result[i] = v[i] / num;
	}
	return result;

}

#endif