#pragma once

#include<vector>

#include"Mechanics/surface_mesh.h"
#include"Mechanics/surface_mesh_tip.h"

namespace MS{
	int simulation_start(MS::surface_mesh &sm, std::vector<MS::filament_tip*> &tips);
}