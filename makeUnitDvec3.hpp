#ifndef makeUnitDvec3_h
#define makeUnitDvec3_h


#include <glm/ext/vector_double3.hpp> // for dvec3 type, 3-component vectors of double precision
#include "dvec3Mag.hpp" // for getting vector magnitude 


// takes a vector of type dvec3 and returns the same vector scaled by its magnitude
glm::dvec3 makeUnitDvec3(glm::dvec3 v){
	double magnitude = dvec3Mag(v);
	
	for(int i = 0; i < 3; ++i){
		v[i]= v[i] / magnitude;
	}
	return v;

}

#endif