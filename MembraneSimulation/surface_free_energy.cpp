#include<math.h>

#include"surface_free_energy.h"
#include"surface_mesh.h"

// TODO: Find experimental results for constants
const double k_c = 1e-19; // Bending modulus
const double k_g = -2*k_c; // Saddle-splay modulus
const double c_0 = 0.0; // Spontaneous curvature
const double gamma = 0.4; // Surface tension


double MS::h_curv_h(vertex * v) {
	return 2 * k_c*(v->curv_h - c_0)*(v->curv_h - c_0)*v->area;
}
double MS::d_h_curv_h(vertex *v, int c_index) {
	double ans = 0;
	switch (c_index) {
	case 0:
		ans += 4 * k_c*(v->curv_h - c_0)*v->dx_curv_h*v->area + 2 * k_c*(v->curv_h - c_0)*(v->curv_h - c_0)*v->dx_area;
		for each(vertex * n in v->n) {
			int i = n->neighbor_indices_map[v];
			ans += 4 * k_c*(n->curv_h - c_0)*n->dxn_curv_h[i] * n->area + 2 * k_c*(n->curv_h - c_0)*(n->curv_h - c_0)*n->dxn_area[i];
		}
		break;
	case 1:
		ans += 4 * k_c*(v->curv_h - c_0)*v->dy_curv_h*v->area + 2 * k_c*(v->curv_h - c_0)*(v->curv_h - c_0)*v->dy_area;
		for each(vertex * n in v->n) {
			int i = n->neighbor_indices_map[v];
			ans += 4 * k_c*(n->curv_h - c_0)*n->dyn_curv_h[i] * n->area + 2 * k_c*(n->curv_h - c_0)*(n->curv_h - c_0)*n->dyn_area[i];
		}
		break;
	case 2:
		ans += 4 * k_c*(v->curv_h - c_0)*v->dz_curv_h*v->area + 2 * k_c*(v->curv_h - c_0)*(v->curv_h - c_0)*v->dz_area;
		for each(vertex * n in v->n) {
			int i = n->neighbor_indices_map[v];
			ans += 4 * k_c*(n->curv_h - c_0)*n->dzn_curv_h[i] * n->area + 2 * k_c*(n->curv_h - c_0)*(n->curv_h - c_0)*n->dzn_area[i];
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
			int i = n->neighbor_indices_map[v];
			ans += k_g*(n->dxn_curv_g[i] * n->area + n->curv_g * n->dxn_area[i]);
		}
		break;
	case 1:
		ans += k_g*(v->dy_curv_g*v->area + v->curv_g*v->dy_area);
		for each(vertex * n in v->n) {
			int i = n->neighbor_indices_map[v];
			ans += k_g*(n->dyn_curv_g[i] * n->area + n->curv_g * n->dyn_area[i]);
		}
		break;
	case 2:
		ans += k_g*(v->dz_curv_g*v->area + v->curv_g*v->dz_area);
		for each(vertex * n in v->n) {
			int i = n->neighbor_indices_map[v];
			ans += k_g*(n->dzn_curv_g[i] * n->area + n->curv_g * n->dzn_area[i]);
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
			ans += gamma*n->dxn_area[n->neighbor_indices_map[v]];
		}
		break;
	case 1:
		ans += gamma*v->dy_area;
		for each(vertex * n in v->n) {
			ans += gamma*n->dyn_area[n->neighbor_indices_map[v]];
		}
		break;
	case 2:
		ans += gamma*v->dz_area;
		for each(vertex * n in v->n) {
			ans += gamma*n->dzn_area[n->neighbor_indices_map[v]];
		}
		break;
	}
	return ans;
}
double MS::h_pressure(vertex *v) {
	if (v->area < 0.0001*v->area0)return gamma * log(10000); // Maximum value
	return -gamma * v->area0 * log(v->area / v->area0);
}
double MS::d_h_pressure(vertex *v, int c_index) {
	double ans = 0;
	if (v->area < 0.0001*v->area0)return 0;
	switch (c_index) {
	case 0:
		ans += -gamma * v->area0 / v->area * v->dx_area;
		for each(vertex * n in v->n) {
			ans += -gamma * n->area0 / n->area * n->dxn_area[n->neighbor_indices_map[v]];
		}
		break;
	case 1:
		ans += -gamma * v->area0 / v->area * v->dy_area;
		for each(vertex * n in v->n) {
			ans += -gamma * n->area0 / n->area * n->dyn_area[n->neighbor_indices_map[v]];
		}
		break;
	case 2:
		ans += -gamma * v->area0 / v->area * v->dz_area;
		for each(vertex * n in v->n) {
			ans += -gamma * n->area0 / n->area * n->dzn_area[n->neighbor_indices_map[v]];
		}
		break;
	}
	return ans;
}
double MS::h_surface_quadratic(vertex *v) {
	return gamma / 2 / v->area0 * (v->area - v->area0) * (v->area - v->area0);
}
double MS::d_h_surface_quadratic(vertex *v, int c_index) {
	double ans = 0;
	switch (c_index) {
	case 0:
		ans += gamma / v->area0 * (v->area - v->area0) * v->dx_area;
		for each(vertex * n in v->n) {
			ans += gamma / n->area0 * (n->area - n->area0) * n->dxn_area[n->neighbor_indices_map[v]];
		}
		break;
	case 1:
		ans += gamma / v->area0 * (v->area - v->area0) * v->dy_area;
		for each(vertex * n in v->n) {
			ans += gamma / n->area0 * (n->area - n->area0) * n->dyn_area[n->neighbor_indices_map[v]];
		}
		break;
	case 2:
		ans += gamma / v->area0 * (v->area - v->area0) * v->dz_area;
		for each(vertex * n in v->n) {
			ans += gamma / n->area0 * (n->area - n->area0) * n->dzn_area[n->neighbor_indices_map[v]];
		}
		break;
	}
	return ans;
}

double MS::h_all(vertex * v) {
	/*
		Energy caused by Gaussian curvature could be neglected because for a closed surface
		it is a constant (Gauss-Bonnet theorem).
	*/
	double h_surf = (QUADRATIC_SURFACE_ENERGY ? h_surface_quadratic(v) : h_tension(v) + h_pressure(v));

	return h_curv_h(v) + h_surf + h_point_interact_v(NULL, v);
}
double MS::d_h_all(vertex * v, int c_index) {
	double d_h_surf = (QUADRATIC_SURFACE_ENERGY ? d_h_surface_quadratic(v, c_index) : d_h_tension(v, c_index) + d_h_pressure(v, c_index));

	return d_h_curv_h(v, c_index)
		+ d_h_surf
		+ d_h_potential(v, c_index);
}

double polymer_len = 0;
math_public::Vec3 *po = new math_public::Vec3(polymer_len, 0, 0);

double MS::update_len(double param) {
	polymer_len = param;
	po->x = polymer_len;
	return polymer_len;
}
double MS::h_point_interact_v(math_public::Vec3 *p, vertex * v) {
	double r = dist(*(v->point), *po);
	double R = 1e-7;
	double ep = 1e-15;
	if (r > R*pow(2, 1.0 / 6))return -ep;
	return 4 * ep*(pow(R / r, 12) - pow(R / r, 6));

}

double MS::d_h_potential(vertex * v, int c_index) {
	double r = dist(*(v->point), *po);
	double R = 1e-7;
	double ep = 1e-15;
	if (r > R*pow(2, 1.0 / 6))return -ep;

	double dr;
	switch (c_index) {
	case 0:
		dr = (v->point->x - po->x) / r;
		break;
	case 1:
		dr = (v->point->y - po->y) / r;
		break;
	case 2:
		dr = (v->point->z - po->z) / r;

	}
	return 4 * ep*(-12 * pow(R / r, 12) + 6 * pow(R / r, 6)) / r * dr;
}
