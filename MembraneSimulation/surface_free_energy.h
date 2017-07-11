/*
	This file contains a definition of a point, as would be simulating a filament tip.
	Currently the point has absolute position, and no counter-interaction is imposed on the point.
*/

#pragma once

#include"math_public.h"
#include"surface_mesh.h"

namespace MS {

	extern math_public::Vec3 *po;
	double update_len(double param); // should be removed when considering real polymers

	extern std::vector<facet*> po_neighbor;
}
