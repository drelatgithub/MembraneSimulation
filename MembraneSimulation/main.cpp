#define _CRT_SECURE_NO_WARNINGS

#include<iostream>

#include"mesh_initialization.h"
#include"surface_mesh.h"
#include"simulation_process.h"

int main() {
	std::vector<MS::vertex*> vertices;

	mesh_init(vertices);

	// starting simulation
	std::cout << "Simulation starting...\n";
	MS::simulation_start(vertices);
	


	system("pause");

	return 0;
}