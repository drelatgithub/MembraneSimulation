/*
Loading an already defined mesh file into the data structure.
*/

#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<sys/stat.h>

#include"surface_mesh.h"
#include"simulation_process.h"

bool mesh_init(std::vector<MS::vertex*> &vertices) {
	bool success = false;

	char *position_file = "position.txt";
	char *neighbors_file = "neighbors.txt";

	struct stat buffer;
	if (stat(position_file,&buffer)==0 && stat(neighbors_file, &buffer) == 0) { // File exists
		std::cout << "Saved mesh file found. Trying to build from file...\n";

		std::ifstream position_in;
		position_in.open(position_file);
		std::ifstream neighbors_in;
		neighbors_in.open(neighbors_file);

		std::string line;

		// getting positions
		while (std::getline(position_in, line)) {
			std::stringstream ss(line);
			double x, y, z;
			if (ss >> x >> y >> z) {
				MS::point_3 *pt = new MS::point_3(x, y, z);
				MS::vertex *new_vertex = new MS::vertex(pt);
				vertices.push_back(new_vertex);
			}
		}

		// getting neighbors (and count vertices)
		int num = 0;
		while (std::getline(neighbors_in, line)) {
			std::stringstream ss(line);
			int n;
			while (ss >> n) {
				vertices[num]->n.push_back(vertices[n]);
				vertices[num]->dump_data_vectors();
			}
			vertices[num]->gen_next_prev_n();
			num++;
		}

		std::cout << "Number of vertices: " << num << std::endl;

		position_in.close();
		neighbors_in.close();

		success = true;

	}
	else {
		std::cout << "At least one file needed for mesh data is not found." << std::endl;
	}

	return success;
}