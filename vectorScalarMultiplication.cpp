#ifndef vectorScalarMultiplication_h 
#define vectorScalarMultiplication_h 

#include <cmath> 
#include <vector>

using namespace std; 

vector<double> vectorScalarMultiplication(vector<double> v, double num){
	int vectorSize = v.size();
	vector<double> result(vectorSize);
	
	for(int i = 0; i < vectorSize; ++i){
		result[i] = v[i] * num;
	}
	return result;
}
#endif
