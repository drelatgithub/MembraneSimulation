#define _CRT_SECURE_NO_WARNINGS

#include<iostream>

#include"mesh_initialization.h"
#include"surface_mesh.h"
#include"simulation_process.h"

int main() {
	std::vector<MS::vertex*> vertices;
	std::vector<MS::facet*> facets;

	if (mesh_init(vertices, facets)) {

		// starting simulation
		std::cout << "Simulation starting...\n";
		MS::simulation_start(vertices, facets);

	}


	system("pause");

	return 0;
}