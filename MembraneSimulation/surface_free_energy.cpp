#include"surface_free_energy.h"
#include"surface_mesh.h"

// TODO: Find experimental results for constants
const double k_c = 1.0; // Bending modulus
const double k_g = 1.0; // Saddle-splay modulus
const double c_0 = 0.0; // Spontaneous curvature
const double gamma = 0.01; // Surface tension

double MS::h_curv_h(vertex * v) {
	return 2 * k_c*(v->curv_h - c_0)*(v->curv_h - c_0);
}
double MS::d_h_curv_h(vertex *v, int c_index) {
	switch (c_index) {
	case 0:
		return 4 * k_c*(v->curv_h - c_0)*v->dx_curv_h;
	case 1:
		return 4 * k_c*(v->curv_h - c_0)*v->dy_curv_h;
	case 2:
		return 4 * k_c*(v->curv_h - c_0)*v->dz_curv_h;
	}
}

double MS::h_curv_g(vertex * v) {
	return k_g * v->curv_g;
}
double MS::d_h_curv_g(vertex *v, int c_index) {
	switch (c_index) {
	case 0:
		return k_g*v->dx_curv_g;
	case 1:
		return k_g*v->dy_curv_g;
	case 2:
		return k_g*v->dz_curv_g;
	}
}

double MS::h_tension(vertex * v) {
	return gamma * (v->area - v->area0);
}
double MS::d_h_tension(vertex *v, int c_index) {
	switch (c_index) {
	case 0:
		return gamma*v->dx_area;
	case 1:
		return gamma*v->dy_area;
	case 2:
		return gamma*v->dz_area;
	}
}

double MS::h_all(vertex * v) {
	return h_tension(v);
	return h_curv_h(v) + h_curv_g(v) + h_tension(v);
}
double MS::d_h_all(vertex * v, int c_index) {
	return d_h_tension(v, c_index);
	return d_h_curv_h(v, c_index)
		+ d_h_curv_g(v, c_index)
		+ d_h_tension(v, c_index);
}