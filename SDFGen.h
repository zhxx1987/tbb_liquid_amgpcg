#ifndef SDFGEN_H
#define SDFGEN_H
#include "makelevelset3.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <limits>
class SDFGenerator
{
public:
	SDFGenerator()
	{
	}
	~SDFGenerator(){}
	int ni, nj, nk;

	Array3f sdf;
	float dx;
	
	float getPhi(Vec3f pos)
	{
		return interpolate_value(pos/dx,sdf);
	}
	void objToSDF(char* filename, Vec3f min_box, Vec3f max_box)
	{
		int padding = 2;

		std::ifstream infile(filename);
		if(!infile) {
			std::cerr << "Failed to open. Terminating.\n";
			exit(-1);
		}

		int ignored_lines = 0;
		std::string line;
		std::vector<Vec3f> vertList;
		std::vector<Vec3ui> faceList;
		while(!infile.eof()) {
			std::getline(infile, line);
			if(line.substr(0,1) == std::string("v")) {
				std::stringstream data(line);
				char c;
				Vec3f point;
				data >> c >> point[0] >> point[1] >> point[2];
				vertList.push_back(point);
				update_minmax(point, min_box, max_box);
			}
			else if(line.substr(0,1) == std::string("f")) {
				std::stringstream data(line);
				char c;
				int v0,v1,v2;
				data >> c >> v0 >> v1 >> v2;
				faceList.push_back(Vec3ui(v0-1,v1-1,v2-1));
			}
			else {
				++ignored_lines; 
			}
		}
		infile.close();

		if(ignored_lines > 0)
			std::cout << "Warning: " << ignored_lines << " lines were ignored since they did not contain faces or vertices.\n";

		std::cout << "Read in " << vertList.size() << " vertices and " << faceList.size() << " faces." << std::endl;

		//Add padding around the box.
		Vec3f unit(1,1,1);
		/*min_box -= padding*dx*unit;
		max_box += padding*dx*unit;*/
		Vec3ui sizes = Vec3ui((max_box - min_box)/dx);

		std::cout << "Bound box size: (" << min_box << ") to (" << max_box << ") with dimensions " << sizes << "." << std::endl;
		std::cout << "Computing signed distance field.\n";
		
		ni = sizes[0];
		nj = sizes[1];
		nk = sizes[2];
		make_level_set3(faceList, vertList, min_box, dx, sizes[0], sizes[1], sizes[2], sdf);


	}
};


#endif