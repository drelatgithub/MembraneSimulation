#include<math.h>

#include"surface_free_energy.h"
#include"surface_mesh.h"

using namespace math_public;


const double k_c = 1e-19; // Bending modulus
const double k_g = -2*k_c; // Saddle-splay modulus
const double c_0 = 0.0; // Spontaneous curvature
const double gamma = 0.4; // Surface tension


double MS::h_curv_h(vertex * v) {
	return 2 * k_c*(v->curv_h - c_0)*(v->curv_h - c_0)*v->area;
}
Vec3 MS::d_h_curv_h(vertex *v) {
	Vec3 ans;
	ans += 4 * k_c * (v->curv_h - c_0) * v->d_curv_h * v->area + 2 * k_c * (v->curv_h - c_0) * (v->curv_h - c_0) * v->d_area;
	for each(vertex* n in v->n) {
		int i = n->neighbor_indices_map[v];
		ans += 4 * k_c * (n->curv_h - c_0) * n->dn_curv_h[i] * n->area + 2 * k_c * (n->curv_h - c_0) * (n->curv_h - c_0) * n->dn_area[i];
	}
	return ans;
}
double MS::h_curv_g(vertex * v) {
	return k_g * v->curv_g * v->area;
}
Vec3 MS::d_h_curv_g(vertex *v) {
	Vec3 ans;
	ans += k_g * (v->d_curv_g * v->area + v->curv_g * v->d_area);
	for each(vertex* n in v->n) {
		int i = n->neighbor_indices_map[v];
		ans += k_g * (n->dn_curv_g[i] * n->area + n->curv_g * n->dn_area[i]);
	}
	return ans;
}

double MS::h_area(vertex *v) {
	return gamma / 2 / v->area0 * (v->area - v->area0) * (v->area - v->area0);
}
Vec3 MS::d_h_area(vertex *v) {
	Vec3 ans;
	ans += gamma / v->area0 * (v->area - v->area0) * v->d_area;
	for each(vertex* n in v->n) {
		ans += gamma / n->area0 * (n->area - n->area0) * n->dn_area[n->neighbor_indices_map[v]];
	}
	return ans;
}

double MS::h_all(vertex * v) {
	/*
		Energy caused by Gaussian curvature could be neglected because for a closed surface
		it is a constant (Gauss-Bonnet theorem).
	*/

	return h_curv_h(v) + h_area(v) + h_point_interact_v(NULL, v);
}
Vec3 MS::d_h_all(vertex * v) {

	return d_h_curv_h(v)
		+ d_h_area(v)
		+ d_h_potential(v);
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

Vec3 MS::d_h_potential(vertex * v) {
	double r = dist(*(v->point), *po);
	double R = 1e-7;
	double ep = 1e-15;
	if (r > R*pow(2, 1.0 / 6))return Vec3();

	Vec3 dr = (*(v->point) - *po) / r;
	return 4 * ep*(-12 * pow(R / r, 12) + 6 * pow(R / r, 6)) / r * dr;
}
