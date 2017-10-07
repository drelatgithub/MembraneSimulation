#include"surface_mesh.h"
using namespace MS;

void MS::surface_mesh::initialize() {
	size_t N = vertices.size();
	size_t N_f = facets.size();

	for (size_t i = 0; i < N_f; ++i) {
		facets[i]->update_geo();
	}
	for (size_t i = 0; i < N; ++i) {
		vertices[i]->count_neighbors();
		if (i == 0) {
			int a = 1;
		}
		vertices[i]->update_geo();
		vertices[i]->make_initial();
	}
}