/*
Loading an already defined mesh file into the data structure.
*/

#include<fstream>
#include<sstream>
#include<string>
#include<sys/stat.h>

#include"common.h"
#include"math_public.h"
#include"surface_mesh.h"
#include"simulation_process.h"

bool mesh_init(std::vector<MS::vertex*> &vertices, std::vector<MS::facet*> &facets) {
	bool success = false;

	char *position_file = "position.txt";
	char *neighbors_file = "neighbors.txt";

	struct stat buffer;
	if (stat(position_file, &buffer) == 0 && stat(neighbors_file, &buffer) == 0) { // File exists
		LOG(INFO) << "Saved mesh file found. Trying to generate from file...";

		std::ifstream position_in;
		position_in.open(position_file);
		std::ifstream neighbors_in;
		neighbors_in.open(neighbors_file);

		std::string line;

		int num_vertices, num_edges, num_facets;

		// getting positions
		while (std::getline(position_in, line)) {
			std::stringstream ss(line);
			double x, y, z;
			if (ss >> x >> y >> z) {
				math_public::Vec3 *pt = new math_public::Vec3(x, y, z);
				MS::vertex *new_vertex = new MS::vertex(pt);
				vertices.push_back(new_vertex);
			}
		}

		// getting neighbors (and count vertices)
		num_vertices = 0;
		num_edges = 0;
		while (std::getline(neighbors_in, line)) {
			std::stringstream ss(line);
			int n;
			while (ss >> n) {
				vertices[num_vertices]->n.push_back(vertices[n]);
				vertices[num_vertices]->dump_data_vectors();
				num_edges++;
			}
			vertices[num_vertices]->gen_next_prev_n();
			num_vertices++;
		}
		num_edges /= 2;

		LOG(INFO) << "Number of vertices: " << num_vertices << "; Number of edges: " << num_edges;
		int predicted_num_facets = num_edges - num_vertices + 2; // Euler characteristic is 2

		LOG(INFO) << "Registering facets...";
		facets.reserve(predicted_num_facets);
		num_facets = 0;
		for (int i = 0; i < num_vertices; i++) {
			for (int j = 0; j < vertices[i]->neighbors; j++) {
				// Propose a facet
				MS::facet *f = new MS::facet(vertices[i], vertices[i]->n[j], vertices[i]->nn[j]);

				// Check whether this facet has already existed
				bool exist = false;
				for (int k = 0; k < facets.size(); k++) {
					if (*facets[k] == *f) {
						exist = true;
						delete f;
						break;
					}
				}
				if (!exist) {
					facets.push_back(f); // Leave the facet in heap
					f->ind[0] = f->v[0]->neighbor_indices_map[f->v[1]];
					f->ind[1] = f->v[1]->neighbor_indices_map[f->v[2]];
					f->ind[2] = f->v[2]->neighbor_indices_map[f->v[0]];
					num_facets++;
				}
				else {
					LOG(WARNING) << "Check your code! This should never happen.";
				}
			}
		}
		LOG(DEBUG) << "Predicted number of facets: " << predicted_num_facets << "; Number of facets: " << num_facets;
		if (predicted_num_facets != num_facets)
			LOG(WARNING) << "The number of facets (" << num_facets << ") is not as expected (" << predicted_num_facets << ").";

		position_in.close();
		neighbors_in.close();

		success = true;

	}
	else {
		LOG(ERROR) << "At least one file needed for mesh data is not found.";
	}

	return success;
}