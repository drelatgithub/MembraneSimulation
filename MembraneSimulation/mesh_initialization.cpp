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

bool mesh_init(MS::surface_mesh &sm) {
	auto &vertices = sm.vertices;
	auto &facets = sm.facets;
	auto &edges = sm.edges;

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

		LOG(INFO) << "Registering edges and facets...";
		facets.reserve(predicted_num_facets);
		edges.reserve(num_edges);
		num_facets = 0;
		for (int i = 0; i < num_vertices; i++) {
			for (int j = 0; j < vertices[i]->neighbors; j++) {
				// Propose a facet
				MS::facet *f = new MS::facet(vertices[i], vertices[i]->n[j], vertices[i]->nn[j]);
				// Propose an edge
				MS::edge *e = new MS::edge(vertices[i], vertices[i]->n[j]);

				// Check whether this facet has already existed
				bool facet_exist = false, edge_exist = false;
				int k_facet, k_edge;
				for (k_facet = 0; k_facet < facets.size(); k_facet++) {
					if (*facets[k_facet] == *f) {
						facet_exist = true;
						break;
					}
				}
				if (!facet_exist) {
					facets.push_back(f); // Leave the facet in heap
					f->update_geo();
					vertices[i]->f.push_back(f);
					num_facets++;
				}
				else {
					vertices[i]->f.push_back(facets[k_facet]);
					delete f;
				}

				for (k_edge = 0; k_edge < edges.size(); k_edge++) {
					if (*edges[k_edge] == *e) {
						edge_exist = true;
						break;
					}
				}
				if (!edge_exist) {
					edges.push_back(e); // Leave the edge in heap
					//e->update_geo();
					vertices[i]->e.push_back(e);
				}
				else {
					vertices[i]->e.push_back(edges[k_edge]);
					delete e;
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