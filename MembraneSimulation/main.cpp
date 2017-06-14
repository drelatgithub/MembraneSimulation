#define _CRT_SECURE_NO_WARNINGS

#include"common.h"
#include"mesh_initialization.h"
#include"surface_mesh.h"
#include"simulation_process.h"

int main() {
	logger::Logger::default_init("simulation.log");

	std::vector<MS::vertex*> vertices;
	std::vector<MS::facet*> facets;

	if (mesh_init(vertices, facets)) {

		// starting simulation
		LOG(INFO) << "Simulation starting...";
		MS::simulation_start(vertices, facets);

	}


	system("pause");

	return 0;
}