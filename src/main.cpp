#define _CRT_SECURE_NO_WARNINGS

#include"common.h"
#include"mesh_initialization.h"
#include"Mechanics/surface_mesh.h"
#include"Mechanics/surface_mesh_tip.h"
#include"simulation_process.h"
#include"test.h"

int main() {
	logger::Logger::default_init("simulation.log");

	test::run_all_tests();

	MS::surface_mesh sm;
	std::vector<MS::filament_tip*> tips;

	if (mesh_init(sm)) {

		// starting simulation
		LOG(INFO) << "Simulation starting...";
		MS::simulation_start(sm, tips);

	}


	system("pause");

	return 0;
}