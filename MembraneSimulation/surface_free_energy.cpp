#define _USE_MATH_DEFINES

#include<math.h>

#include"common.h"
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

	return h_curv_h(v) + h_area(v);
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
double MS::h_point_interact_facet(math_public::Vec3 *p, facet *f) {
	// The total integral over the whole plane is
	// H = d0^n * d^(2-n) * f(n) * k
	// where k is the coefficient, d0 is the distance of force, f(n) is a polynomial of force.
	static int po_pwr = 6; // negative power of distance in expression of potential energy
	static double r_0_factor = 1.0 / sqrt(po_pwr-1);
	static double d0 = 1e-8;
	static double energy_coe = 1.0; // in J/m^2
	static double d0_pwr = pow(d0, po_pwr);
	static double pre_calc = d0_pwr * energy_coe;

	/**********************************
	find the foot of perpendicular O on the triangle
	r_O = r0 + alpha v1 + beta v2
	**********************************/
	Vec3 v1 = *(f->v[1]->point) - *(f->v[0]->point), v2 = *(f->v[2]->point) - *(f->v[0]->point);
	Vec3 r12 = v2 - v1;
	Vec3 r0p = *p - *(f->v[0]->point);

	// alpha and beta need to satisfy the perpendicular condition
	// A * (alpha, beta)' = B
	// So (alpha, beta)' = A^(-1) * B
	double dot12 = dot(v1, v2);
	v1.calc_norm();
	v2.calc_norm();
	r12.calc_norm();
	double det_A = v1.norm2 * v2.norm2 - dot12 * dot12;
	double AR11 = v2.norm2 / det_A,
		AR12 = -dot12 / det_A,
		AR21 = -dot12 / det_A,
		AR22 = v1.norm2 / det_A;
	double B1 = r0p.dot(v1), B2 = r0p.dot(v2);

	double alpha = AR11*B1 + AR12*B2,
		beta = AR21*B1 + AR22*B2;

	Vec3 r0O = alpha*v1 + beta*v2;
	Vec3 rO = r0O + *(f->v[0]->point),
		r1O = r0O - v1,
		r2O = r0O - v2;
	r0O.calc_norm();
	r1O.calc_norm();
	r2O.calc_norm();
	if (true) {
		Vec3 rOp = *p - rO;
		LOG(DEBUG) << "rOp dot v1 = " << rOp.dot(v1);
		LOG(DEBUG) << "rOp dot v2 = " << rOp.dot(v2);
	}

	/**********************************
	distance from the point to the triangle
	**********************************/
	double d = (*p - rO).get_norm();

	/**********************************
	distance from point O to all 3 edges
	**********************************/
	double d1 = cross(r0O, v1).get_norm() / v1.norm;
	if (beta < 0)d1 = -d1;
	double d2 = cross(r0O, v2).get_norm() / v2.norm;
	if (alpha < 0)d2 = -d2;
	double d3 = cross(r1O, r12).get_norm() / r12.norm;
	if (alpha - beta > 1)d3 = -d3;

	/**********************************
	find the affecting region and calculate energy
	**********************************/
	double sigma = d * r_0_factor;
	if (sigma >= 0.05*(v1.norm + v2.norm + r12.norm))
		LOG(WARNING) << "Sigma is not significantly smaller than the scale of the facet.";

	double I = pow(d, 2 - po_pwr) * pre_calc;
	double vcoe1 = exp(-r0O.norm2 / (sigma*sigma)),
		vcoe2 = exp(-r1O.norm2 / (sigma*sigma)),
		vcoe3 = exp(-r2O.norm2 / (sigma*sigma));
	double ecoe1 = 0.5 + 0.5*tanh(d1 / sigma),
		ecoe2 = 0.5 + 0.5*tanh(d2 / sigma),
		ecoe3 = 0.5 + 0.5*tanh(d3 / sigma);
	double en = vcoe1*(f->v[0]->theta[f->ind[0]] / (2 * M_PI))
		+ vcoe2*(f->v[1]->theta[f->ind[1]] / (2 * M_PI))
		+ vcoe3*(f->v[2]->theta[f->ind[2]] / (2 * M_PI))
		+ (1 - vcoe1)*(1 - vcoe2)*(1 - vcoe3)*ecoe1*ecoe2*ecoe3;
	en *= I;

	return en;
	
}

Vec3 MS::d_h_potential(vertex * v) {
	double r = dist(*(v->point), *po);
	double R = 1e-7;
	double ep = 1e-15;
	if (r > R*pow(2, 1.0 / 6))return Vec3();

	Vec3 dr = (*(v->point) - *po) / r;
	return 4 * ep*(-12 * pow(R / r, 12) + 6 * pow(R / r, 6)) / r * dr;
}
