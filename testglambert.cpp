#include <vector>
#include <iostream>

#include <glm/ext/vector_double3.hpp> // for dvec3 type, 3-component vectors of double precision

#include "glambert.cpp"




int main(){
	using namespace glm;
	dvec3 vi = {0,0,0};
	dvec3 vf = {0,0,0};
	double mu = 1.327 * pow(10,11);
	dvec3 pv1 = {-77448649, 228632436, 6691087};
	dvec3 Vv1 = {-22.069,-5.607, -0.42428};
	dvec3 pv2 = {-120342944, 87422241, 27891};
	dvec3 Vv2 = {-18.132, -24.123, -0.002482};
	
	double tof = 86400*185;
	int nrev = 0;
	
	int check = glambert(&vi, &vf, mu, pv1, Vv1, pv2, Vv2, tof, nrev);
	if (check == -1) {
		std::cout<<"Glambert has returned -1"<<std::endl;
	}
	
	
	using namespace std;
	//vi[1] = 3.14159;
	//vf[1] = 6.283;
	cout<<endl<<endl<< "GLambert inputs:"<<endl;
	cout<< "Position of Body 1: {" << pv1[0] << " " << pv1[1]<< " " << pv1[2]<< "}"<< endl;
	cout<< "Velocity of Body 1: {" << Vv1[0] << " " << Vv1[1]<< " " << Vv1[2]<< "}"<< endl<<endl;
	cout<< "Position of Body 2: {" << pv2[0] << " " << pv2[1]<< " " << pv2[2]<< "}"<< endl;
	cout<< "Velocity of Body 2: {" << Vv2[0] << " " << Vv2[1]<< " " << Vv2[2]<< "}"<< endl<<endl;
	cout<< "mu = "<<mu<<" km^3/s^2" <<endl;
	cout<< "Time of flight: "<<tof<<" seconds"<<endl<<endl;
	
	cout<<"GLambert has returned "<<check<<endl;
	
	cout<< "vi = {" << vi[0] << " " << vi[1]<< " " << vi[2]<< "}"<< endl;
	cout<< "vf = {" << vf[0] << " " << vf[1]<< " " << vf[2]<< "}"<< endl;
	
	return 1;
}