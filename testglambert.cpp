#include <vector> // for indexing ephemeris data
#include <iostream> // for printing computational results to screen
#include <fstream> // for reading data from csv
#include <sstream> // for reading data from csv


#include <glm/ext/vector_double3.hpp> // for dvec3 type, 3-component vectors of double precision

#include "glambert.cpp"
#include "printdvec3.hpp"




int main(){
	
	int count = 0;
	int count_body_one = 0;
	int count_body_two = 0;
	std::string x_str, y_str, z_str, vx_str, vy_str, vz_str, date_str, JDTDB_str, lt_str, rg_str, rr_str, file1, file2;
	double x, y, z, vx, vy, vz, jd;
	std::string line;
    bool start_reading = false;
    std::vector<glm::dvec3> positions_e;
    std::vector<glm::dvec3> velocities_e;
	std::vector<glm::dvec3> positions_v;
    std::vector<glm::dvec3> velocities_v;
	std::vector<std::string> bodyOneDate;
	std::vector<std::string> bodyTwoDate;
	std::vector<double> JDbodyOne;
	std::vector<double> JDbodyTwo;
	
	// get filenames to open
	std::cout<<"Enter the filename of the first ephemeris data: ";
	std::cin>>file1;
	std::cout<<"Enter the filename of the second ephemeris data: ";
	std::cin>>file2;
	
	// read in first body's data
	std::ifstream infile_one(file1);
    if (!infile_one.is_open()) {
        std::cerr << "Error opening file 1." << std::endl;
        return -1;
    }

   // Read the file line by line
    	
    while (std::getline(infile_one, line)) {
        // Ignore lines before the "$$SOE" line
        if (line.find("$$SOE") != std::string::npos) {
            start_reading = true;
            continue;
        }
		
		// ignore lines after "$$EOE" and stop reading
		if (line.find("$$EOE") != std::string::npos) {
			start_reading = false;
			break;
        }

        // Parse the line if we're past the "$$SOE" line
        if (start_reading) {
			// accept csv data and store each entry as a string
            std::istringstream iss(line);
			std::getline(iss, JDTDB_str, ',');
			std::getline(iss, date_str, ',');
			std::getline(iss, x_str, ',');
			std::getline(iss, y_str, ',');
			std::getline(iss, z_str, ',');
			std::getline(iss, vx_str, ',');
			std::getline(iss, vy_str, ',');
			std::getline(iss, vz_str, ',');
			std::getline(iss, lt_str, ',');
			std::getline(iss, rg_str, ',');
			std::getline(iss, rr_str, ',');
			
			// convert the strings to doubles
			jd = std::stod(JDTDB_str);
			x = std::stod(x_str);
			y = std::stod(y_str);
			z = std::stod(z_str);
			vx = std::stod(vx_str);
			vy = std::stod(vy_str);
			vz = std::stod(vz_str);
			
			// assemble the doubles into dvec3 arrays
            positions_e.emplace_back(x, y, z);
            velocities_e.emplace_back(vx, vy, vz);
			bodyOneDate.emplace_back(date_str);
			JDbodyOne.emplace_back(jd);
			++count_body_one;
        }
    }
    infile_one.close();
	
	
	// read in second body's data
	
	
	std::ifstream infile_two(file2);
    if (!infile_two.is_open()) {
        std::cerr << "Error opening file 2." << std::endl;
        return -1;
    }
	
	    while (std::getline(infile_two, line)) {
        // Ignore lines before the "$$SOE" line
        if (line.find("$$SOE") != std::string::npos) {
            start_reading = true;
            continue;
        }
		
		// ignore lines after "$$EOE" and stop reading
		if (line.find("$$EOE") != std::string::npos) {
			start_reading = false;
			break;
        }

        // Parse the line if we're past the "$$SOE" line
        if (start_reading) {
			// accept csv data and store each entry as a string
            std::istringstream iss(line);
			std::getline(iss, JDTDB_str, ',');
			std::getline(iss, date_str, ',');
			std::getline(iss, x_str, ',');
			std::getline(iss, y_str, ',');
			std::getline(iss, z_str, ',');
			std::getline(iss, vx_str, ',');
			std::getline(iss, vy_str, ',');
			std::getline(iss, vz_str, ',');
			std::getline(iss, lt_str, ',');
			std::getline(iss, rg_str, ',');
			std::getline(iss, rr_str, ',');
			
			// convert values of interest to doubles
			jd = std::stod(JDTDB_str);
			x = std::stod(x_str);
			y = std::stod(y_str);
			z = std::stod(z_str);
			vx = std::stod(vx_str);
			vy = std::stod(vy_str);
			vz = std::stod(vz_str);
			
			// assemble dvec3 arrays
            positions_v.emplace_back(x, y, z);
            velocities_v.emplace_back(vx, vy, vz);
			bodyTwoDate.emplace_back(date_str);
			JDbodyTwo.emplace_back(jd);
			++count_body_two;
        }
    }
    infile_two.close();
	
	
	if (count_body_one != count_body_two){
		std::cout<<"You have different numbers of datum's for the two bodies."<<std::endl<<std::endl;
		if (count_body_one < count_body_two){
			count = count_body_one;
		} else {
			count = count_body_two;
		}
	} else {
		count = count_body_one;
	};
	
	glm::dvec3 vi = {0,0,0};
	glm::dvec3 vf = {0,0,0};
	double mu = 1.327 * pow(10,11);
	
	double tof = 86400*185;
	int nrev = 0;
	int check = -1; // errors if not reassigned properly
	
	
	
	for (int i = 0; i < (count); ++i){
		
		tof = (JDbodyTwo[i] - JDbodyOne[i]) * 86400; // turn difference in dates into time of flight by multiplying by seconds per day
		if (tof == 0){
			std::cout<<"Instantaneous travel not possible."<<std::endl;
			tof = 1;
		}
		
		check = glambert(&vi, &vf, mu, positions_e[i], velocities_e[i], positions_v[i], velocities_v[i], tof, nrev);
		if (check == -1) {
			std::cout<<"Glambert has returned -1"<<std::endl;
		} else {
			std::cout<<"Departing on "<<bodyOneDate[i]<<std::endl<<"Arriving on "<<bodyTwoDate[i]<<std::endl;
			std::cout<<"vi = ";
			printdvec3(vi);
			std::cout<<"vf = ";
			printdvec3(vf);
			std::cout<<std::endl;
		}
	}
	
	
	return check;
}
