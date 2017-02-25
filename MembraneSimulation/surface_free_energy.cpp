#include<math.h>

#include"surface_free_energy.h"
#include"surface_mesh.h"

// TODO: Find experimental results for constants
const double k_c = 2.0; // Bending modulus
const double k_g = -4.0; // Saddle-splay modulus
const double c_0 = 0.0; // Spontaneous curvature
const double gamma = 1.0; // Surface tension
const double k_ps = 0.02; // Surface pressure constant. For balance around area0, k_ps = area0 * gamma
const double h_max_ps = log(1000)*k_ps;


double polymer_len = 0.6;

double MS::h_curv_h(vertex * v) {
	return 2 * k_c*(v->curv_h - c_0)*(v->curv_h - c_0)*v->area;
}
double MS::d_h_curv_h(vertex *v, int c_index) {
	double ans = 0;
	switch (c_index) {
	case 0:
		ans += 4 * k_c*(v->curv_h - c_0)*v->dx_curv_h*v->area + 2 * k_c*(v->curv_h - c_0)*(v->curv_h - c_0)*v->dx_area;
		for each(vertex * n in v->n) {
			int i = n->neighbour_indices_map[v];
			ans += 4 * k_c*(n->curv_h - c_0)*n->dxn_curv_h[i] * v->area + 2 * k_c*(v->curv_h - c_0)*(v->curv_h - c_0)*n->dxn_area[i];
		}
		break;
	case 1:
		ans += 4 * k_c*(v->curv_h - c_0)*v->dy_curv_h*v->area + 2 * k_c*(v->curv_h - c_0)*(v->curv_h - c_0)*v->dy_area;
		for each(vertex * n in v->n) {
			int i = n->neighbour_indices_map[v];
			ans += 4 * k_c*(n->curv_h - c_0)*n->dyn_curv_h[i] * v->area + 2 * k_c*(v->curv_h - c_0)*(v->curv_h - c_0)*n->dyn_area[i];
		}
		break;
	case 2:
		ans += 4 * k_c*(v->curv_h - c_0)*v->dz_curv_h*v->area + 2 * k_c*(v->curv_h - c_0)*(v->curv_h - c_0)*v->dz_area;
		for each(vertex * n in v->n) {
			int i = n->neighbour_indices_map[v];
			ans += 4 * k_c*(n->curv_h - c_0)*n->dzn_curv_h[i] * v->area + 2 * k_c*(v->curv_h - c_0)*(v->curv_h - c_0)*n->dzn_area[i];
		}
		break;
	}
	return ans;
}

double MS::h_curv_g(vertex * v) {
	return k_g * v->curv_g * v->area;
}
double MS::d_h_curv_g(vertex *v, int c_index) {
	double ans = 0;
	switch (c_index) {
	case 0:
		ans += k_g*(v->dx_curv_g*v->area + v->curv_g*v->dx_area);
		for each(vertex * n in v->n) {
			int i = n->neighbour_indices_map[v];
			ans += k_g*(n->dxn_curv_g[i] * v->area + v->curv_g*n->dxn_area[i]);
		}
		break;
	case 1:
		ans += k_g*(v->dy_curv_g*v->area + v->curv_g*v->dy_area);
		for each(vertex * n in v->n) {
			int i = n->neighbour_indices_map[v];
			ans += k_g*(n->dyn_curv_g[i] * v->area + v->curv_g*n->dyn_area[i]);
		}
		break;
	case 2:
		ans += k_g*(v->dz_curv_g*v->area + v->curv_g*v->dz_area);
		for each(vertex * n in v->n) {
			int i = n->neighbour_indices_map[v];
			ans += k_g*(n->dzn_curv_g[i] * v->area + v->curv_g*n->dzn_area[i]);
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

double MS::h_pressure(vertex *v) {
	if (v->area < 0.001*v->area0)return h_max_ps;
	return -k_ps*log(v->area / v->area0);
}
double MS::d_h_pressure(vertex *v, int c_index) {
	double ans = 0;
	if (v->area < 0.001*v->area0)return 0;
	switch (c_index) {
	case 0:
		ans += -k_ps / v->area * v->dx_area;
		for each(vertex * n in v->n) {
			ans += -k_ps / v->area * n->dxn_area[n->neighbour_indices_map[v]];
		}
		break;
	case 1:
		ans += -k_ps / v->area * v->dy_area;
		for each(vertex * n in v->n) {
			ans += -k_ps / v->area * n->dyn_area[n->neighbour_indices_map[v]];
		}
		break;
	case 2:
		ans += -k_ps / v->area * v->dz_area;
		for each(vertex * n in v->n) {
			ans += -k_ps / v->area * n->dzn_area[n->neighbour_indices_map[v]];
		}
		break;
	}
	return ans;
}

double MS::h_all(vertex * v) {
	//return h_tension(v);
	return h_curv_h(v) + h_curv_g(v) + h_tension(v) + h_pressure(v) + h_potential(v);
}
double MS::d_h_all(vertex * v, int c_index) {
	//return d_h_tension(v, c_index);
	return d_h_curv_h(v, c_index)
		+ d_h_curv_g(v, c_index)
		+ d_h_tension(v, c_index)
		+ d_h_pressure(v, c_index)
		+ d_h_potential(v, c_index);
}

double MS::update_len(double param) {
	polymer_len = param;
	return polymer_len;
}
double MS::h_potential(vertex * v) {
	if (v->point->x < polymer_len - 0.6)return exp(10) - 1;
	else if (v->point->x < polymer_len)return exp((polymer_len - v->point->x) / 0.06) - 1;
	else return 0;
}

double MS::d_h_potential(vertex * v, int c_index) {
	if (c_index == 1 || c_index == 2)return 0;
	if (v->point->x >= polymer_len - 0.6 && v->point->x < polymer_len)return -exp((polymer_len - v->point->x) / 0.06) / 0.06;
	else return 0;
}
