#ifndef printdvec3_h
#define printdvec3_h

#include <glm/ext/vector_double3.hpp> // for dvec3 type, 3-component vectors of double precision
#include <iostream> // for cout

// prints the 3 elements of a vector type dvec3
void printdvec3(glm::dvec3 v){
	std::cout<<"{ ";
	for(int i = 0;i<3;++i){
		std::cout<<v[i]<<" ";
	}
	std::cout<<"}"<<std::endl;
}

#endif