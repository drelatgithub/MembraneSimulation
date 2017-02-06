#include"surface_free_energy.h"
#include"surface_mesh.h"

// TODO: Find experimental results for constants
const double k_c = 0.01; // Bending modulus
const double k_g = 0.01; // Saddle-splay modulus
const double c_0 = 0.0; // Spontaneous curvature
const double gamma = 0.01; // Surface tension

double MS::h_curv_h(vertex * v) {
	return 2 * k_c*(v->curv_h - c_0)*(v->curv_h - c_0);
}
double MS::d_h_curv_h(vertex *v, int c_index) {
	double ans = 0;
	switch (c_index) {
	case 0:
		ans += 4 * k_c*(v->curv_h - c_0)*v->dx_curv_h;
		for each(vertex * n in v->n) {
			ans += 4 * k_c*(n->curv_h - c_0)*n->dxn_curv_h[n->neighbour_indices_map[v]];
		}
		break;
	case 1:
		ans += 4 * k_c*(v->curv_h - c_0)*v->dy_curv_h;
		for each(vertex * n in v->n) {
			ans += 4 * k_c*(n->curv_h - c_0)*n->dyn_curv_h[n->neighbour_indices_map[v]];
		}
		break;
	case 2:
		ans += 4 * k_c*(v->curv_h - c_0)*v->dz_curv_h;
		for each(vertex * n in v->n) {
			ans += 4 * k_c*(n->curv_h - c_0)*n->dzn_curv_h[n->neighbour_indices_map[v]];
		}
		break;
	}
	return ans;
}

double MS::h_curv_g(vertex * v) {
	return k_g * v->curv_g;
}
double MS::d_h_curv_g(vertex *v, int c_index) {
	double ans = 0;
	switch (c_index) {
	case 0:
		ans += k_g*v->dx_curv_g;
		for each(vertex * n in v->n) {
			ans += k_g*n->dxn_curv_g[n->neighbour_indices_map[v]];
		}
		break;
	case 1:
		ans += k_g*v->dy_curv_g;
		for each(vertex * n in v->n) {
			ans += k_g*n->dyn_curv_g[n->neighbour_indices_map[v]];
		}
		break;
	case 2:
		ans += k_g*v->dz_curv_g;
		for each(vertex * n in v->n) {
			ans += k_g*n->dzn_curv_g[n->neighbour_indices_map[v]];
		}
		break;
	}
	return ans;
}

double MS::h_tension(vertex * v) {
	return gamma * (v->area - v->area0);
}
double MS::d_h_tension(vertex *v, int c_index) {
	double ans = 0;
	switch (c_index) {
	case 0:
		ans += gamma*v->dx_area;
		for each(vertex * n in v->n) {
			ans += gamma*n->dxn_area[n->neighbour_indices_map[v]];
		}
		break;
	case 1:
		ans += gamma*v->dy_area;
		for each(vertex * n in v->n) {
			ans += gamma*n->dyn_area[n->neighbour_indices_map[v]];
		}
		break;
	case 2:
		ans += gamma*v->dz_area;
		for each(vertex * n in v->n) {
			ans += gamma*n->dzn_area[n->neighbour_indices_map[v]];
		}
		break;
	}
	return ans;
}

double MS::h_all(vertex * v) {
	//return h_tension(v);
	return h_curv_h(v) + h_curv_g(v) + h_tension(v);
}
double MS::d_h_all(vertex * v, int c_index) {
	//return d_h_tension(v, c_index);
	return d_h_curv_h(v, c_index)
		+ d_h_curv_g(v, c_index)
		+ d_h_tension(v, c_index);
}