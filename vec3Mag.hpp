#ifndef vec3Mag_h
#define vec3Mag_h


#include <glm/ext/vector_double3.hpp> // for dvec3 type, 3-component vectors of double precision

// returns the magnitude of a vector of type dvec3
double vec3Mag(glm::dvec3 v){
	double sum = 0;
	
	for(int i = 0; i < 3; ++i){
		sum += v[i] * v[i];
	}
	return sqrt(sum);

}

#endif