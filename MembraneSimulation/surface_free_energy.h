#pragma once

#include"surface_mesh.h"

#define QUADRATIC_SURFACE_ENERGY false

namespace MS {
	/********************************************
		Free energies
	********************************************/
	double h_curv_h(vertex * v); // Free energy - Curvature - Mean
	double h_curv_g(vertex * v); // Free energy - Curvature - Gaussian
	// surface energies use either quadratic potential or tension + pressure depending on QUADRATIC_SURFACE_ENERGY
	double h_tension(vertex * v); // Free energy - Surface tension
	double h_pressure(vertex *v); // Free energy - Surface pressure
	double h_surface_quadratic(vertex *v); // Free energy - Surface tension and pressure near the lowest energy
	//double h_osmotic(vertex * v); // Free energy - Osmotic pressure
	double h_point_interact_v(math_public::Vec3 *p, vertex * v); // Free energy - Potential energy

	/********************************************
		Partial derivatives of free energies
	********************************************/
	double d_h_curv_h(vertex *v, int c_index);
	double d_h_curv_g(vertex *v, int c_index);
	double d_h_tension(vertex *v, int c_index);
	double d_h_pressure(vertex *v, int c_index);
	double d_h_surface_quadratic(vertex *v, int c_index);
	//double d_h_
	double d_h_potential(vertex * v, int c_index);

	double h_all(vertex * v);
	double d_h_all(vertex *v, int c_index);

	double update_len(double param); // should be removed when considering real polymers
}
