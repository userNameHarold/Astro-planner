#include <vector> // for indexing ephemeris data
#include <iostream> // for printing computational results to screen
#include <fstream> // for reading data from csv
#include <sstream> // for reading data from csv
#include <cfloat> // for DBL_MAX

#include <glm/ext/vector_double3.hpp> // for dvec3 type, 3-component vectors of double precision

#include "glambert.hpp" // function GLambert
#include "dvec3Mag.hpp" // for DeltaV calculations




int main(){
	
	int count = 0;
	int loopCount = 0;
	int count_body_one = 0;
	int count_body_two = 0;
	std::string x_str, y_str, z_str, vx_str, vy_str, vz_str, date_str, JDTDB_str, lt_str, rg_str, rr_str, foo, file1, file2;
	double x, y, z, vx, vy, vz, jd;
	std::string line;
    bool start_reading = false;
	bool foundName = false;
    std::vector<glm::dvec3> positions_b1;
    std::vector<glm::dvec3> velocities_b1;
	std::vector<glm::dvec3> positions_b2;
    std::vector<glm::dvec3> velocities_b2;
	std::vector<std::string> bodyOneDate;
	std::vector<std::string> bodyTwoDate;
	std::vector<double> JDbodyOne;
	std::vector<double> JDbodyTwo;
	
	
	std::string body1name;
	std::string body2name;
	std::string target = "Target body name: ";
	size_t pos = 0;
	size_t end_pos = 0;
	
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

		// get the name of the planet in question
		if (!foundName){
			pos = line.find(target);
			if (pos != std::string::npos){
				pos += target.size();
				end_pos = line.find("(", pos);
				body1name = line.substr(pos, end_pos - pos);
				foundName = true;
			}
		}

		
        // Ignore lines before the "$$SOE" line
        if (line.find("$$SOE") != std::string::npos) {
            start_reading = true;
            continue;
        }
		
		// ignore lines after "$$EOE" and stop reading
		if (line.find("$$EOE") != std::string::npos) {

			start_reading = false; // set up bools for next read
			foundName = false; 
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

            positions_b1.emplace_back(x, y, z);
            velocities_b1.emplace_back(vx, vy, vz);
			bodyOneDate.emplace_back(date_str);
			JDbodyOne.emplace_back(jd);
			++count_body_one;
        }

		
		++loopCount;
    }
    infile_one.close();
	loopCount = 0;
	
	
	// read in second body's data
	
	
	std::ifstream infile_two(file2);
    if (!infile_two.is_open()) {
        std::cerr << "Error opening file 2." << std::endl;
        return -1;
    }
	

	while (std::getline(infile_two, line)) {
		
		// get name of second body
		if (!foundName){
		pos = line.find(target);
		if (pos != std::string::npos){
			pos += target.size();
			end_pos = line.find("(", pos);
			body2name = line.substr(pos, end_pos - pos);
			foundName = true;
			}
		}
		
		
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

			positions_b2.emplace_back(x, y, z);
			velocities_b2.emplace_back(vx, vy, vz);
			bodyTwoDate.emplace_back(date_str);
			JDbodyTwo.emplace_back(jd);
			++count_body_two;
		}
    }
    infile_two.close();
	
	std::string outFileName = body1name + "To " + body2name + "Transfer.txt";
	
// the following is an outdated feature from a rpevious version. It is left in as comments for the situation where a 1:1 related between depature and arrival dates is needed	
	
	// if (count_body_one != count_body_two){
		// std::cout<<"You have different numbers of datum's for the two bodies."<<std::endl<<std::endl;
		// if (count_body_one < count_body_two){
			// count = count_body_one;
		// } else {
			// count = count_body_two;
		// }
	// } else {
		// count = count_body_one;
	// };
	
	glm::dvec3 vi = {0,0,0};
	glm::dvec3 vf = {0,0,0};
	glm::dvec3 deltaV1 = {0,0,0};
	glm::dvec3 deltaV2 = {0,0,0};
	double mu = 1.327 * pow(10,11);
	double DV1 = 0.0;
	double DV2 = 0.0;
	double DVTot = 0.0;
	double maxDV = 0.0;
	std::vector<int> maxDVpos = {0,0};
	double minDV = DBL_MAX;
	std::vector<int> minDVpos = {0,0};
	
	double tof = 86400*185;
	int nrev = 0;
	int check = -1; // errors if not reassigned properly

	
	std::ofstream outFile;
	outFile.open(outFileName, std::ios::app);
	if (!outFile.is_open()) {
        std::cerr << "Error opening "<<outFileName << std::endl;
        return -1;
    }
	outFile<<"Transfer from "<<body1name<<" departing "<<bodyOneDate[0]<< " to "<<bodyOneDate[count_body_one - 1]<<std::endl;
	outFile<<"To "<<body2name<<" arriving "<<bodyTwoDate[0]<< " to "<< bodyTwoDate[count_body_two - 1]<<std::endl<<std::endl;
	
	
	for (int i = 0; i < count_body_one; ++i){
		
		for (int j = 0; j < count_body_two; ++j){
		
			tof = (JDbodyTwo[j] - JDbodyOne[i]) * 86400; // turn difference in dates into time of flight by multiplying by seconds per day
			if (tof == 0){
				//std::cout<<"Instantaneous travel not possible."<<std::endl;
				tof = 1;
			}
			
			if (JDbodyTwo[j] < JDbodyOne[i]){
				// arrival date before departure date, this is not possible
			} else {
				
				check = glambert(&vi, &vf, mu, positions_b1[i], velocities_b1[i], positions_b2[j], velocities_b2[j], tof, nrev);
				if (check == -1) {
					std::cout<<"Glambert has returned -1"<<std::endl;
				} else {
					
					// calculate delta V requirements
					for (int k = 0; k < 3; ++k){
						deltaV1[k] = fabs(vi[k] - velocities_b1[i][k]);
						deltaV2[k] = fabs(vf[k] - velocities_b2[j][k]);
					}
					DV1 = dvec3Mag(deltaV1);
					DV2 = dvec3Mag(deltaV2);
					DVTot = DV1 + DV2;
					
					
					if (DVTot < minDV){
						minDV = DVTot;
						minDVpos[0] = i;
						minDVpos[1] = j;
					}
					if (DVTot > maxDV){
						maxDV = DVTot;
						maxDVpos[0] = i;
						maxDVpos[1] = j;
					}
						
					outFile<<std::endl<<"Departing on "<<bodyOneDate[i]<<std::endl<<"Arriving on  "<<bodyTwoDate[j]<<std::endl<<std::endl;
					outFile<<"      vi         = {"<<vi[0]<<" "<<vi[1]<<" "<<vi[2]<<"} km/s"<<std::endl;
					outFile<<"      delta V1   = {"<<deltaV1[0]<<" "<<deltaV1[1]<<" "<<deltaV1[2]<<"} km/s"<<std::endl;
					outFile<<"      |delta V1| = "<< DV1<<"km/s"<<std::endl;
					outFile<<"      vf         = {"<<vf[0]<<" "<<vf[1]<<" "<<vf[2]<<"} km/s"<<std::endl;
					outFile<<"      delta V2   = {"<<deltaV2[0]<<" "<<deltaV2[1]<<" "<<deltaV2[2]<<"} km/s"<<std::endl;
					outFile<<"      |delta v2| = "<<DV2<<"km/s"<<std::endl;
					outFile<<"For a total of "<<DVTot<<"km/s of delta V"<<std::endl<<std::endl;
					
				}
			}
		}
	}
	
	outFile<<std::endl<<"The best transfer takes "<<minDV<<"km/s of delta V, and occurs when leaving "<<body1name<<"on "<<bodyOneDate[minDVpos[0]]<<" and arriving at "<<body2name<<"on "<<bodyTwoDate[minDVpos[1]]<<std::endl;
	outFile<<"The worst transfer takes "<<maxDV<<"km/s of dela V, and occurs when leaving "<<body1name<<"on "<<bodyOneDate[maxDVpos[0]]<<" and arriving at "<<body2name<<"on "<<bodyTwoDate[maxDVpos[1]]<<std::endl;
	outFile.close();
	
	
	/* for (int i = 0; i < (count); ++i){
		
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
	} */
	std::cout<<std::endl<<"The best transfer takes "<<minDV<<"km/s of delta V, and occurs when leaving "<<body1name<<"on "<<bodyOneDate[minDVpos[0]]<<" and arriving at "<<body2name<<"on "<<bodyTwoDate[minDVpos[1]]<<std::endl;
	std::cout<<"The worst transfer takes "<<maxDV<<"km/s of dela V, and occurs when leaving "<<body1name<<"on "<<bodyOneDate[maxDVpos[0]]<<" and arriving at "<<body2name<<"on "<<bodyTwoDate[maxDVpos[1]]<<std::endl;

	
	
	return check;
}
