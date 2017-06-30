#pragma once

#include"math_public.h"
#include"surface_mesh.h"

namespace MS {
	/********************************************
		Free energies
	********************************************/
	double h_curv_h(vertex * v); // Free energy - Curvature - Mean
	double h_curv_g(vertex * v); // Free energy - Curvature - Gaussian
	double h_area(vertex *v); // Free energy - Surface tesnsion and pressure
	//double h_osmotic(vertex * v); // Free energy - Osmotic pressure
	double h_point_interact_v(math_public::Vec3 *p, vertex * v); // Free energy - Potential energy

	/********************************************
		Partial derivatives of free energies
	********************************************/
	math_public::Vec3 d_h_curv_h(vertex *v);
	math_public::Vec3 d_h_curv_h(vertex *v);
	math_public::Vec3 d_h_curv_g(vertex *v);
	math_public::Vec3 d_h_area(vertex *v);
	//double d_h_
	math_public::Vec3 d_h_potential(vertex * v);

	double h_all(vertex * v);
	math_public::Vec3 d_h_all(vertex *v);

	double update_len(double param); // should be removed when considering real polymers
}
