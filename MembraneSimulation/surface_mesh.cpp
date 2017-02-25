#define _USE_MATH_DEFINES
#include<math.h>
#include<iostream>

#include"surface_mesh.h"

#define USE_VONOROI_CELL true

MS::point_3::point_3(double nx, double ny, double nz) {
	x = nx, y = ny, z = nz;
}

MS::vertex::vertex(MS::point_3 *npoint) {
	point = npoint;
	point_last = new point_3(0,0,0);
}

MS::vertex::~vertex() {
	delete point_last;
}

int MS::vertex::count_neighbours() {
	neighbours = n.size();
	for (int i = 0; i < neighbours; i++) {
		neighbour_indices_map[n[i]] = i;
	}
	return neighbours;
}

int MS::vertex::dump_data_vectors() {
	// theta
	theta.push_back(0), dx_theta.push_back(0), dy_theta.push_back(0), dz_theta.push_back(0);
	dxn_theta.push_back(0), dyn_theta.push_back(0), dzn_theta.push_back(0);
	dxnn_theta.push_back(0), dynn_theta.push_back(0), dznn_theta.push_back(0);
	sin_theta.push_back(0), dx_sin_theta.push_back(0), dy_sin_theta.push_back(0), dz_sin_theta.push_back(0);
	dxn_sin_theta.push_back(0), dyn_sin_theta.push_back(0), dzn_sin_theta.push_back(0);
	dxnn_sin_theta.push_back(0), dynn_sin_theta.push_back(0), dznn_sin_theta.push_back(0);
	// theta2
	theta2.push_back(0), dx_theta2.push_back(0), dy_theta2.push_back(0), dz_theta2.push_back(0);
	dxn_theta2.push_back(0), dyn_theta2.push_back(0), dzn_theta2.push_back(0);
	dxnp_theta2.push_back(0), dynp_theta2.push_back(0), dznp_theta2.push_back(0);
	cot_theta2.push_back(0), dx_cot_theta2.push_back(0), dy_cot_theta2.push_back(0), dz_cot_theta2.push_back(0);
	dxn_cot_theta2.push_back(0), dyn_cot_theta2.push_back(0), dzn_cot_theta2.push_back(0);
	dxnp_cot_theta2.push_back(0), dynp_cot_theta2.push_back(0), dznp_cot_theta2.push_back(0);
	// theta3
	theta3.push_back(0), dx_theta3.push_back(0), dy_theta3.push_back(0), dz_theta3.push_back(0);
	dxn_theta3.push_back(0), dyn_theta3.push_back(0), dzn_theta3.push_back(0);
	dxnn_theta3.push_back(0), dynn_theta3.push_back(0), dznn_theta3.push_back(0);
	cot_theta3.push_back(0), dx_cot_theta3.push_back(0), dy_cot_theta3.push_back(0), dz_cot_theta3.push_back(0);
	dxn_cot_theta3.push_back(0), dyn_cot_theta3.push_back(0), dzn_cot_theta3.push_back(0);
	dxnn_cot_theta3.push_back(0), dynn_cot_theta3.push_back(0), dznn_cot_theta3.push_back(0);
	// distances
	r_p_n.push_back(0), dx_r_p_n.push_back(0), dy_r_p_n.push_back(0), dz_r_p_n.push_back(0), dxn_r_p_n.push_back(0), dyn_r_p_n.push_back(0), dzn_r_p_n.push_back(0);
	r_p_n_prev.push_back(0), dx_r_p_n_prev.push_back(0), dy_r_p_n_prev.push_back(0), dz_r_p_n_prev.push_back(0), dxnp_r_p_n_prev.push_back(0), dynp_r_p_n_prev.push_back(0), dznp_r_p_n_prev.push_back(0);
	r_p_n_next.push_back(0), dx_r_p_n_next.push_back(0), dy_r_p_n_next.push_back(0), dz_r_p_n_next.push_back(0), dxnn_r_p_n_next.push_back(0), dynn_r_p_n_next.push_back(0), dznn_r_p_n_next.push_back(0);

	// Derivatives around a vertex
	dxn_area.push_back(0), dyn_area.push_back(0), dzn_area.push_back(0);
	dxn_curv_h.push_back(0), dyn_curv_h.push_back(0), dzn_curv_h.push_back(0);
	dxn_curv_g.push_back(0), dyn_curv_g.push_back(0), dzn_curv_g.push_back(0);

	return 0;
}

/*
int MS::vertex::fill_vectors_with_zeroes() {
	for (int i = 0, len = n.size(); i < len; i++) {
		theta[i] = 0; dx_theta[i] = dy_theta[i] = dz_theta[i] = 0;
		sin_theta[i] = 0; dx_sin_theta[i] = dy_sin_theta[i] = dz_sin_theta[i] = 0;
		theta2[i] = 0; dx_theta2[i] = dy_theta2[i] = dz_theta2[i] = 0;
		cot_theta2[i] = 0; dx_cot_theta2[i] = dy_cot_theta2[i] = dz_cot_theta2[i] = 0;
		theta3[i] = 0; dx_theta3[i] = dy_theta3[i] = dz_theta3[i] = 0;
		cot_theta3[i] = 0; dx_cot_theta3[i] = dy_cot_theta3[i] = dz_cot_theta3[i] = 0;
		r_p_n[i] = 0; dx_r_p_n[i] = dy_r_p_n[i] = dz_r_p_n[i] = 0;
		r_p_n_prev[i] = 0; dx_r_p_n_prev[i] = dy_r_p_n_prev[i] = dz_r_p_n_prev[i] = 0;
		r_p_n_next[i] = 0; dx_r_p_n_next[i] = dy_r_p_n_next[i] = dz_r_p_n_next[i] = 0;
		return 0;
	}
}
*/


double distance2(const MS::point_3 *p, const MS::point_3 *np) {
	return (p->x - np->x)*(p->x - np->x) + (p->y - np->y)*(p->y - np->y) + (p->z - np->z)*(p->z - np->z);
}
double distance(const MS::point_3 *p, const MS::point_3 *np) {
	return sqrt(distance2(p, np));
}


void MS::vertex::calc_angle() {
	for (int i = 0; i < neighbours; i++) {
		// Distances
		r_p_n[i] = distance(point, n[i]->point);
		dx_r_p_n[i] = (point->x - n[i]->point->x) / r_p_n[i];
		dy_r_p_n[i] = (point->y - n[i]->point->y) / r_p_n[i];
		dz_r_p_n[i] = (point->z - n[i]->point->z) / r_p_n[i];
		dxn_r_p_n[i] = (n[i]->point->x - point->x) / r_p_n[i];
		dyn_r_p_n[i] = (n[i]->point->y - point->y) / r_p_n[i];
		dzn_r_p_n[i] = (n[i]->point->z - point->z) / r_p_n[i];
		r_p_n_prev[i] = distance(point, n_prev[i]->point);
		dx_r_p_n_prev[i] = (point->x - n_prev[i]->point->x) / r_p_n_prev[i];
		dy_r_p_n_prev[i] = (point->y - n_prev[i]->point->y) / r_p_n_prev[i];
		dz_r_p_n_prev[i] = (point->z - n_prev[i]->point->z) / r_p_n_prev[i];
		dxnp_r_p_n_prev[i] = (n_prev[i]->point->x - point->x) / r_p_n_prev[i];
		dynp_r_p_n_prev[i] = (n_prev[i]->point->y - point->y) / r_p_n_prev[i];
		dznp_r_p_n_prev[i] = (n_prev[i]->point->z - point->z) / r_p_n_prev[i];
		r_p_n_next[i] = distance(point, n_next[i]->point);
		dx_r_p_n_next[i] = (point->x - n_next[i]->point->x) / r_p_n_next[i];
		dy_r_p_n_next[i] = (point->y - n_next[i]->point->y) / r_p_n_next[i];
		dz_r_p_n_next[i] = (point->z - n_next[i]->point->z) / r_p_n_next[i];
		dxnn_r_p_n_next[i] = (n_next[i]->point->x - point->x) / r_p_n_next[i];
		dynn_r_p_n_next[i] = (n_next[i]->point->y - point->y) / r_p_n_next[i];
		dznn_r_p_n_next[i] = (n_next[i]->point->z - point->z) / r_p_n_next[i];
		//std::cout << r_p_n[i];

		// Find theta
		double inner_product = (n_next[i]->point->x - point->x)*(n[i]->point->x - point->x) + (n_next[i]->point->y - point->y)*(n[i]->point->y - point->y) + (n_next[i]->point->z - point->z)*(n[i]->point->z - point->z);
		double dx_inner_product = 2 * point->x - n_next[i]->point->x - n[i]->point->x;
		double dy_inner_product = 2 * point->y - n_next[i]->point->y - n[i]->point->y;
		double dz_inner_product = 2 * point->z - n_next[i]->point->z - n[i]->point->z;
		double dxn_inner_product = n_next[i]->point->x - point->x;
		double dyn_inner_product = n_next[i]->point->y - point->y;
		double dzn_inner_product = n_next[i]->point->z - point->z;
		double dxnn_inner_product = n[i]->point->x - point->x;
		double dynn_inner_product = n[i]->point->y - point->y;
		double dznn_inner_product = n[i]->point->z - point->z;
		double cos_theta = inner_product / (r_p_n[i] * r_p_n_next[i]);
		double dx_cos_theta = (r_p_n[i] * r_p_n_next[i] * dx_inner_product - inner_product*(r_p_n[i] * dx_r_p_n_next[i] + dx_r_p_n[i] * r_p_n_next[i])) / (r_p_n[i] * r_p_n[i] * r_p_n_next[i] * r_p_n_next[i]);
		double dy_cos_theta = (r_p_n[i] * r_p_n_next[i] * dy_inner_product - inner_product*(r_p_n[i] * dy_r_p_n_next[i] + dy_r_p_n[i] * r_p_n_next[i])) / (r_p_n[i] * r_p_n[i] * r_p_n_next[i] * r_p_n_next[i]);
		double dz_cos_theta = (r_p_n[i] * r_p_n_next[i] * dz_inner_product - inner_product*(r_p_n[i] * dz_r_p_n_next[i] + dz_r_p_n[i] * r_p_n_next[i])) / (r_p_n[i] * r_p_n[i] * r_p_n_next[i] * r_p_n_next[i]);
		double dxn_cos_theta = (r_p_n[i] * dxn_inner_product - inner_product*dxn_r_p_n[i]) / (r_p_n[i] * r_p_n[i] * r_p_n_next[i]);
		double dyn_cos_theta = (r_p_n[i] * dyn_inner_product - inner_product*dyn_r_p_n[i]) / (r_p_n[i] * r_p_n[i] * r_p_n_next[i]);
		double dzn_cos_theta = (r_p_n[i] * dzn_inner_product - inner_product*dzn_r_p_n[i]) / (r_p_n[i] * r_p_n[i] * r_p_n_next[i]);
		double dxnn_cos_theta = (r_p_n_next[i] * dxnn_inner_product - inner_product*dxnn_r_p_n_next[i]) / (r_p_n_next[i] * r_p_n_next[i] * r_p_n[i]);
		double dynn_cos_theta = (r_p_n_next[i] * dynn_inner_product - inner_product*dynn_r_p_n_next[i]) / (r_p_n_next[i] * r_p_n_next[i] * r_p_n[i]);
		double dznn_cos_theta = (r_p_n_next[i] * dznn_inner_product - inner_product*dznn_r_p_n_next[i]) / (r_p_n_next[i] * r_p_n_next[i] * r_p_n[i]);
		sin_theta[i] = sqrt(1 - cos_theta*cos_theta);
		theta[i] = acos(cos_theta);
		dx_theta[i] = -dx_cos_theta / sin_theta[i];
		dy_theta[i] = -dy_cos_theta / sin_theta[i];
		dz_theta[i] = -dz_cos_theta / sin_theta[i];
		dxn_theta[i] = -dxn_cos_theta / sin_theta[i];
		dyn_theta[i] = -dyn_cos_theta / sin_theta[i];
		dzn_theta[i] = -dzn_cos_theta / sin_theta[i];
		dxnn_theta[i] = -dxnn_cos_theta / sin_theta[i];
		dynn_theta[i] = -dynn_cos_theta / sin_theta[i];
		dznn_theta[i] = -dznn_cos_theta / sin_theta[i];
		dx_sin_theta[i] = cos_theta*dx_theta[i];
		dy_sin_theta[i] = cos_theta*dy_theta[i];
		dz_sin_theta[i] = cos_theta*dz_theta[i];
		dxn_sin_theta[i] = cos_theta*dxn_theta[i];
		dyn_sin_theta[i] = cos_theta*dyn_theta[i];
		dzn_sin_theta[i] = cos_theta*dzn_theta[i];
		dxnn_sin_theta[i] = cos_theta*dxnn_theta[i];
		dynn_sin_theta[i] = cos_theta*dynn_theta[i];
		dznn_sin_theta[i] = cos_theta*dznn_theta[i];
		//std::cout << theta[i] << '\t';

		// Find theta2
		double r_n_n_prev = distance(n[i]->point, n_prev[i]->point); // derivative with regard to p is 0
		double dxn_r_n_n_prev = (n[i]->point->x - n_prev[i]->point->x) / r_n_n_prev;
		double dyn_r_n_n_prev = (n[i]->point->y - n_prev[i]->point->y) / r_n_n_prev;
		double dzn_r_n_n_prev = (n[i]->point->z - n_prev[i]->point->z) / r_n_n_prev;
		double dxnp_r_n_n_prev = (n_prev[i]->point->x - n[i]->point->x) / r_n_n_prev;
		double dynp_r_n_n_prev = (n_prev[i]->point->y - n[i]->point->y) / r_n_n_prev;
		double dznp_r_n_n_prev = (n_prev[i]->point->z - n[i]->point->z) / r_n_n_prev;
		double inner_product2 = (n[i]->point->x - n_prev[i]->point->x)*(point->x - n_prev[i]->point->x) + (n[i]->point->y - n_prev[i]->point->y)*(point->y - n_prev[i]->point->y) + (n[i]->point->z - n_prev[i]->point->z)*(point->z - n_prev[i]->point->z);
		double dx_inner_product2 = n[i]->point->x - n_prev[i]->point->x;
		double dy_inner_product2 = n[i]->point->y - n_prev[i]->point->y;
		double dz_inner_product2 = n[i]->point->z - n_prev[i]->point->z;
		double dxn_inner_product2 = point->x - n_prev[i]->point->x;
		double dyn_inner_product2 = point->y - n_prev[i]->point->y;
		double dzn_inner_product2 = point->z - n_prev[i]->point->z;
		double dxnp_inner_product2 = 2 * n_prev[i]->point->x - point->x - n[i]->point->x;
		double dynp_inner_product2 = 2 * n_prev[i]->point->y - point->y - n[i]->point->y;
		double dznp_inner_product2 = 2 * n_prev[i]->point->z - point->z - n[i]->point->z;
		double cos_theta2 = inner_product2 / (r_p_n_prev[i] * r_n_n_prev);
		double dx_cos_theta2 = (r_p_n_prev[i] * dx_inner_product2 - inner_product2*(dx_r_p_n_prev[i])) / (r_p_n_prev[i] * r_p_n_prev[i] * r_n_n_prev);
		double dy_cos_theta2 = (r_p_n_prev[i] * dy_inner_product2 - inner_product2*(dy_r_p_n_prev[i])) / (r_p_n_prev[i] * r_p_n_prev[i] * r_n_n_prev);
		double dz_cos_theta2 = (r_p_n_prev[i] * dz_inner_product2 - inner_product2*(dz_r_p_n_prev[i])) / (r_p_n_prev[i] * r_p_n_prev[i] * r_n_n_prev);
		double dxn_cos_theta2 = (r_n_n_prev*dxn_inner_product2 - inner_product2*dxn_r_n_n_prev) / (r_n_n_prev*r_n_n_prev*r_p_n_prev[i]);
		double dyn_cos_theta2 = (r_n_n_prev*dyn_inner_product2 - inner_product2*dyn_r_n_n_prev) / (r_n_n_prev*r_n_n_prev*r_p_n_prev[i]);
		double dzn_cos_theta2 = (r_n_n_prev*dzn_inner_product2 - inner_product2*dzn_r_n_n_prev) / (r_n_n_prev*r_n_n_prev*r_p_n_prev[i]);
		double dxnp_cos_theta2 = (r_p_n_prev[i] * r_n_n_prev*dxnp_inner_product2 - inner_product2*(dxnp_r_p_n_prev[i] * r_n_n_prev + r_p_n_prev[i] * dxnp_r_n_n_prev)) / (r_p_n_prev[i] * r_p_n_prev[i] * r_n_n_prev*r_n_n_prev);
		double dynp_cos_theta2 = (r_p_n_prev[i] * r_n_n_prev*dynp_inner_product2 - inner_product2*(dynp_r_p_n_prev[i] * r_n_n_prev + r_p_n_prev[i] * dynp_r_n_n_prev)) / (r_p_n_prev[i] * r_p_n_prev[i] * r_n_n_prev*r_n_n_prev);
		double dznp_cos_theta2 = (r_p_n_prev[i] * r_n_n_prev*dznp_inner_product2 - inner_product2*(dznp_r_p_n_prev[i] * r_n_n_prev + r_p_n_prev[i] * dznp_r_n_n_prev)) / (r_p_n_prev[i] * r_p_n_prev[i] * r_n_n_prev*r_n_n_prev);
		double sin_theta2 = sqrt(1 - cos_theta2*cos_theta2);
		theta2[i] = acos(cos_theta2);
		dx_theta2[i] = -dx_cos_theta2 / sin_theta2;
		dy_theta2[i] = -dy_cos_theta2 / sin_theta2;
		dz_theta2[i] = -dz_cos_theta2 / sin_theta2;
		dxn_theta2[i] = -dxn_cos_theta2 / sin_theta2;
		dyn_theta2[i] = -dyn_cos_theta2 / sin_theta2;
		dzn_theta2[i] = -dzn_cos_theta2 / sin_theta2;
		dxnp_theta2[i] = -dxnp_cos_theta2 / sin_theta2;
		dynp_theta2[i] = -dynp_cos_theta2 / sin_theta2;
		dznp_theta2[i] = -dznp_cos_theta2 / sin_theta2;
		cot_theta2[i] = cos_theta2 / sin_theta2;
		dx_cot_theta2[i] = -dx_theta2[i] / (sin_theta2*sin_theta2);
		dy_cot_theta2[i] = -dy_theta2[i] / (sin_theta2*sin_theta2);
		dz_cot_theta2[i] = -dz_theta2[i] / (sin_theta2*sin_theta2);
		dxn_cot_theta2[i] = -dxn_theta2[i] / (sin_theta2*sin_theta2);
		dyn_cot_theta2[i] = -dyn_theta2[i] / (sin_theta2*sin_theta2);
		dzn_cot_theta2[i] = -dzn_theta2[i] / (sin_theta2*sin_theta2);
		dxnp_cot_theta2[i] = -dxnp_theta2[i] / (sin_theta2*sin_theta2);
		dynp_cot_theta2[i] = -dynp_theta2[i] / (sin_theta2*sin_theta2);
		dznp_cot_theta2[i] = -dznp_theta2[i] / (sin_theta2*sin_theta2);

		// Find theta3
		double r_n_n_next = distance(n[i]->point, n_next[i]->point); // derivative with regard to p is 0
		double dxn_r_n_n_next = (n[i]->point->x - n_next[i]->point->x) / r_n_n_next;
		double dyn_r_n_n_next = (n[i]->point->y - n_next[i]->point->y) / r_n_n_next;
		double dzn_r_n_n_next = (n[i]->point->z - n_next[i]->point->z) / r_n_n_next;
		double dxnn_r_n_n_next = (n_next[i]->point->x - n[i]->point->x) / r_n_n_next;
		double dynn_r_n_n_next = (n_next[i]->point->y - n[i]->point->y) / r_n_n_next;
		double dznn_r_n_n_next = (n_next[i]->point->z - n[i]->point->z) / r_n_n_next;
		double inner_product3 = (n[i]->point->x - n_next[i]->point->x)*(point->x - n_next[i]->point->x) + (n[i]->point->y - n_next[i]->point->y)*(point->y - n_next[i]->point->y) + (n[i]->point->z - n_next[i]->point->z)*(point->z - n_next[i]->point->z);
		double dx_inner_product3 = n[i]->point->x - n_next[i]->point->x;
		double dy_inner_product3 = n[i]->point->y - n_next[i]->point->y;
		double dz_inner_product3 = n[i]->point->z - n_next[i]->point->z;
		double dxn_inner_product3 = point->x - n_next[i]->point->x;
		double dyn_inner_product3 = point->y - n_next[i]->point->y;
		double dzn_inner_product3 = point->z - n_next[i]->point->z;
		double dxnn_inner_product3 = 2 * n_next[i]->point->x - point->x - n[i]->point->x;
		double dynn_inner_product3 = 2 * n_next[i]->point->y - point->y - n[i]->point->y;
		double dznn_inner_product3 = 2 * n_next[i]->point->z - point->z - n[i]->point->z;
		double cos_theta3 = inner_product3 / (r_p_n_next[i] * r_n_n_next);
		double dx_cos_theta3 = (r_p_n_next[i] * dx_inner_product3 - inner_product3*(dx_r_p_n_next[i])) / (r_p_n_next[i] * r_p_n_next[i] * r_n_n_next);
		double dy_cos_theta3 = (r_p_n_next[i] * dy_inner_product3 - inner_product3*(dy_r_p_n_next[i])) / (r_p_n_next[i] * r_p_n_next[i] * r_n_n_next);
		double dz_cos_theta3 = (r_p_n_next[i] * dz_inner_product3 - inner_product3*(dz_r_p_n_next[i])) / (r_p_n_next[i] * r_p_n_next[i] * r_n_n_next);
		double dxn_cos_theta3 = (r_n_n_next*dxn_inner_product3 - inner_product3*dxn_r_n_n_next) / (r_n_n_next*r_n_n_next*r_p_n_next[i]);
		double dyn_cos_theta3 = (r_n_n_next*dyn_inner_product3 - inner_product3*dyn_r_n_n_next) / (r_n_n_next*r_n_n_next*r_p_n_next[i]);
		double dzn_cos_theta3 = (r_n_n_next*dzn_inner_product3 - inner_product3*dzn_r_n_n_next) / (r_n_n_next*r_n_n_next*r_p_n_next[i]);
		double dxnn_cos_theta3 = (r_p_n_next[i] * r_n_n_next*dxnn_inner_product3 - inner_product3*(dxnn_r_p_n_next[i] * r_n_n_next + r_p_n_next[i] * dxnn_r_n_n_next)) / (r_p_n_next[i] * r_p_n_next[i] * r_n_n_next*r_n_n_next);
		double dynn_cos_theta3 = (r_p_n_next[i] * r_n_n_next*dynn_inner_product3 - inner_product3*(dynn_r_p_n_next[i] * r_n_n_next + r_p_n_next[i] * dynn_r_n_n_next)) / (r_p_n_next[i] * r_p_n_next[i] * r_n_n_next*r_n_n_next);
		double dznn_cos_theta3 = (r_p_n_next[i] * r_n_n_next*dznn_inner_product3 - inner_product3*(dznn_r_p_n_next[i] * r_n_n_next + r_p_n_next[i] * dznn_r_n_n_next)) / (r_p_n_next[i] * r_p_n_next[i] * r_n_n_next*r_n_n_next);
		double sin_theta3 = sqrt(1 - cos_theta3*cos_theta3);
		theta3[i] = acos(cos_theta3);
		dx_theta3[i] = -dx_cos_theta3 / sin_theta3;
		dy_theta3[i] = -dy_cos_theta3 / sin_theta3;
		dz_theta3[i] = -dz_cos_theta3 / sin_theta3;
		dxn_theta3[i] = -dxn_cos_theta3 / sin_theta3;
		dyn_theta3[i] = -dyn_cos_theta3 / sin_theta3;
		dzn_theta3[i] = -dzn_cos_theta3 / sin_theta3;
		dxnn_theta3[i] = -dxnn_cos_theta3 / sin_theta3;
		dynn_theta3[i] = -dynn_cos_theta3 / sin_theta3;
		dznn_theta3[i] = -dznn_cos_theta3 / sin_theta3;
		cot_theta3[i] = cos_theta3 / sin_theta3;
		dx_cot_theta3[i] = -dx_theta3[i] / (sin_theta3*sin_theta3);
		dy_cot_theta3[i] = -dy_theta3[i] / (sin_theta3*sin_theta3);
		dz_cot_theta3[i] = -dz_theta3[i] / (sin_theta3*sin_theta3);
		dxn_cot_theta3[i] = -dxn_theta3[i] / (sin_theta3*sin_theta3);
		dyn_cot_theta3[i] = -dyn_theta3[i] / (sin_theta3*sin_theta3);
		dzn_cot_theta3[i] = -dzn_theta3[i] / (sin_theta3*sin_theta3);
		dxnn_cot_theta3[i] = -dxnn_theta3[i] / (sin_theta3*sin_theta3);
		dynn_cot_theta3[i] = -dynn_theta3[i] / (sin_theta3*sin_theta3);
		dznn_cot_theta3[i] = -dznn_theta3[i] / (sin_theta3*sin_theta3);
	}
}

double MS::vertex::calc_area() {
	/*****************************************************************************
	Must be used after angles are calculated.

	Use only Vonoroi area even for obtuse triangles.

	The area around vertex i is:
	A = 1/8 * \sum_{v_j is neighbour of v_i} (cot \alpha_ij + cot \beta_ij ) * |r_i - r_j|^2
	*****************************************************************************/
	if (USE_VONOROI_CELL) {
		area = 0;
		dx_area = dy_area = dz_area = 0;
		for (int i = 0; i < neighbours; i++) {
			dxn_area[i] = dyn_area[i] = dzn_area[i] = 0;
		}
		for (int i = 0; i < neighbours; i++) {
			double dis2 = r_p_n[i] * r_p_n[i];
			area += (cot_theta2[i] + cot_theta3[i])*dis2;
			dx_area += (dx_cot_theta2[i] + dx_cot_theta3[i])*dis2 + (cot_theta2[i] + cot_theta3[i]) * 2 * r_p_n[i] * dx_r_p_n[i];
			dy_area += (dy_cot_theta2[i] + dy_cot_theta3[i])*dis2 + (cot_theta2[i] + cot_theta3[i]) * 2 * r_p_n[i] * dy_r_p_n[i];
			dz_area += (dz_cot_theta2[i] + dz_cot_theta3[i])*dis2 + (cot_theta2[i] + cot_theta3[i]) * 2 * r_p_n[i] * dz_r_p_n[i];
			dxn_area[i] += (dxn_cot_theta2[i] + dxn_cot_theta3[i])*dis2 + (cot_theta2[i] + cot_theta3[i]) * 2 * r_p_n[i] * dxn_r_p_n[i];
			dyn_area[i] += (dyn_cot_theta2[i] + dyn_cot_theta3[i])*dis2 + (cot_theta2[i] + cot_theta3[i]) * 2 * r_p_n[i] * dyn_r_p_n[i];
			dzn_area[i] += (dzn_cot_theta2[i] + dzn_cot_theta3[i])*dis2 + (cot_theta2[i] + cot_theta3[i]) * 2 * r_p_n[i] * dzn_r_p_n[i];
			int i_n = neighbour_indices_map.at(n_next[i]), i_p = neighbour_indices_map.at(n_prev[i]);
			dxn_area[i_n] += (dxnn_cot_theta3[i])*dis2;
			dyn_area[i_n] += (dynn_cot_theta3[i])*dis2;
			dzn_area[i_n] += (dznn_cot_theta3[i])*dis2;
			dxn_area[i_p] += (dxnp_cot_theta2[i])*dis2;
			dyn_area[i_p] += (dynp_cot_theta2[i])*dis2;
			dzn_area[i_p] += (dznp_cot_theta2[i])*dis2;
		}
		area /= 8; dx_area /= 8; dy_area /= 8; dz_area /= 8;
		for (int i = 0; i < neighbours; i++) {
			dxn_area[i] /= 8; dyn_area[i] /= 8; dzn_area[i] /= 8;
		}
		return area;
	}
	else {
		return 0;
	}
}

double MS::vertex::calc_curv_h() {
	/*****************************************************************************
	Must be used after the angles and the area is calculated.
	This function could also calculate the normalized normal vector.

	Use only Vonoroi area even for obtuse triangles.

	Locally, Laplace-Beltrami operator (vector) K = 2 * k_h * n, where n is the normalized normal vector.
	And K = 1/(2A) * \sum_{v_j is neighbour of v_i} (cot \alpha_ij + cot \beta_ij ) * (r_i - r_j)	
	*****************************************************************************/
	if (USE_VONOROI_CELL) {
		double K_x = 0, K_y = 0, K_z = 0; // x, y, z component of K
		double dx_K_x = 0, dy_K_x = 0, dz_K_x = 0, dx_K_y = 0, dy_K_y = 0, dz_K_y = 0, dx_K_z = 0, dy_K_z = 0, dz_K_z = 0;
		std::vector<double> dxn_K_x(neighbours), dyn_K_x(neighbours), dzn_K_x(neighbours);
		std::vector<double> dxn_K_y(neighbours), dyn_K_y(neighbours), dzn_K_y(neighbours);
		std::vector<double> dxn_K_z(neighbours), dyn_K_z(neighbours), dzn_K_z(neighbours);
		for (int i = 0; i < neighbours; i++) {
			double diff_x = point->x - n[i]->point->x, diff_y = point->y - n[i]->point->y, diff_z = point->z - n[i]->point->z;
			K_x += (cot_theta2[i] + cot_theta3[i])*diff_x;
			K_y += (cot_theta2[i] + cot_theta3[i])*diff_y;
			K_z += (cot_theta2[i] + cot_theta3[i])*diff_z;
			dx_K_x += (dx_cot_theta2[i] + dx_cot_theta3[i])*diff_x + (cot_theta2[i] + cot_theta3[i]);
			dy_K_x += (dy_cot_theta2[i] + dy_cot_theta3[i])*diff_x;
			dz_K_x += (dz_cot_theta2[i] + dz_cot_theta3[i])*diff_x;
			dx_K_y += (dx_cot_theta2[i] + dx_cot_theta3[i])*diff_y;
			dy_K_y += (dy_cot_theta2[i] + dy_cot_theta3[i])*diff_y + (cot_theta2[i] + cot_theta3[i]);
			dz_K_y += (dz_cot_theta2[i] + dz_cot_theta3[i])*diff_y;
			dx_K_z += (dx_cot_theta2[i] + dx_cot_theta3[i])*diff_z;
			dy_K_z += (dy_cot_theta2[i] + dy_cot_theta3[i])*diff_z;
			dz_K_z += (dz_cot_theta2[i] + dz_cot_theta3[i])*diff_z + (cot_theta2[i] + cot_theta3[i]);
			dxn_K_x[i] += (dxn_cot_theta2[i] + dxn_cot_theta3[i])*diff_x - (cot_theta2[i] + cot_theta3[i]);
			dyn_K_x[i] += (dyn_cot_theta2[i] + dyn_cot_theta3[i])*diff_x;
			dzn_K_x[i] += (dzn_cot_theta2[i] + dzn_cot_theta3[i])*diff_x;
			dxn_K_y[i] += (dxn_cot_theta2[i] + dxn_cot_theta3[i])*diff_y;
			dyn_K_y[i] += (dyn_cot_theta2[i] + dyn_cot_theta3[i])*diff_y - (cot_theta2[i] + cot_theta3[i]);
			dzn_K_y[i] += (dzn_cot_theta2[i] + dzn_cot_theta3[i])*diff_y;
			dxn_K_z[i] += (dxn_cot_theta2[i] + dxn_cot_theta3[i])*diff_z;
			dyn_K_z[i] += (dyn_cot_theta2[i] + dyn_cot_theta3[i])*diff_z;
			dzn_K_z[i] += (dzn_cot_theta2[i] + dzn_cot_theta3[i])*diff_z - (cot_theta2[i] + cot_theta3[i]);
			int i_n = neighbour_indices_map[n_next[i]], i_p = neighbour_indices_map[n_prev[i]];
			dxn_K_x[i_n] += dxnn_cot_theta3[i] * diff_x;
			dyn_K_x[i_n] += dynn_cot_theta3[i] * diff_x;
			dzn_K_x[i_n] += dznn_cot_theta3[i] * diff_x;
			dxn_K_y[i_n] += dxnn_cot_theta3[i] * diff_y;
			dyn_K_y[i_n] += dynn_cot_theta3[i] * diff_y;
			dzn_K_y[i_n] += dznn_cot_theta3[i] * diff_y;
			dxn_K_z[i_n] += dxnn_cot_theta3[i] * diff_z;
			dyn_K_z[i_n] += dynn_cot_theta3[i] * diff_z;
			dzn_K_z[i_n] += dznn_cot_theta3[i] * diff_z;
			dxn_K_x[i_p] += dxnp_cot_theta2[i] * diff_x;
			dyn_K_x[i_p] += dynp_cot_theta2[i] * diff_x;
			dzn_K_x[i_p] += dznp_cot_theta2[i] * diff_x;
			dxn_K_y[i_p] += dxnp_cot_theta2[i] * diff_y;
			dyn_K_y[i_p] += dynp_cot_theta2[i] * diff_y;
			dzn_K_y[i_p] += dznp_cot_theta2[i] * diff_y;
			dxn_K_z[i_p] += dxnp_cot_theta2[i] * diff_z;
			dyn_K_z[i_p] += dynp_cot_theta2[i] * diff_z;
			dzn_K_z[i_p] += dznp_cot_theta2[i] * diff_z;
		}
		dx_K_x = (area*dx_K_x - K_x*dx_area) / (2 * area*area);
		dy_K_x = (area*dy_K_x - K_x*dy_area) / (2 * area*area);
		dz_K_x = (area*dz_K_x - K_x*dz_area) / (2 * area*area);
		dx_K_y = (area*dx_K_y - K_y*dx_area) / (2 * area*area);
		dy_K_y = (area*dy_K_y - K_y*dy_area) / (2 * area*area);
		dz_K_y = (area*dz_K_y - K_y*dz_area) / (2 * area*area);
		dx_K_z = (area*dx_K_z - K_z*dx_area) / (2 * area*area);
		dy_K_z = (area*dy_K_z - K_z*dy_area) / (2 * area*area);
		dz_K_z = (area*dz_K_z - K_z*dz_area) / (2 * area*area);
		for (int i = 0; i < neighbours; i++) {
			dxn_K_x[i] = (area*dxn_K_x[i] - K_x*dxn_area[i]) / (2 * area*area);
			dyn_K_x[i] = (area*dyn_K_x[i] - K_x*dyn_area[i]) / (2 * area*area);
			dzn_K_x[i] = (area*dzn_K_x[i] - K_x*dzn_area[i]) / (2 * area*area);
			dxn_K_y[i] = (area*dxn_K_y[i] - K_y*dxn_area[i]) / (2 * area*area);
			dyn_K_y[i] = (area*dyn_K_y[i] - K_y*dyn_area[i]) / (2 * area*area);
			dzn_K_y[i] = (area*dzn_K_y[i] - K_y*dzn_area[i]) / (2 * area*area);
			dxn_K_z[i] = (area*dxn_K_z[i] - K_z*dxn_area[i]) / (2 * area*area);
			dyn_K_z[i] = (area*dyn_K_z[i] - K_z*dyn_area[i]) / (2 * area*area);
			dzn_K_z[i] = (area*dzn_K_z[i] - K_z*dzn_area[i]) / (2 * area*area);
		}
		K_x /= 2 * area; K_y /= 2 * area; K_z /= 2 * area;
		curv_h = sqrt(K_x*K_x + K_y*K_y + K_z*K_z) / 2;
		dx_curv_h = (K_x*dx_K_x + K_y*dx_K_y + K_z*dx_K_z) / (4 * curv_h);
		dy_curv_h = (K_x*dy_K_x + K_y*dy_K_y + K_z*dy_K_z) / (4 * curv_h);
		dz_curv_h = (K_x*dz_K_x + K_y*dz_K_y + K_z*dz_K_z) / (4 * curv_h);
		for (int i = 0; i < neighbours; i++) {
			dxn_curv_h[i] = (K_x*dxn_K_x[i] + K_y*dxn_K_y[i] + K_z*dxn_K_z[i]) / (4 * curv_h);
			dyn_curv_h[i] = (K_x*dyn_K_x[i] + K_y*dyn_K_y[i] + K_z*dyn_K_z[i]) / (4 * curv_h);
			dzn_curv_h[i] = (K_x*dzn_K_x[i] + K_y*dzn_K_y[i] + K_z*dzn_K_z[i]) / (4 * curv_h);
		}
		n_x = K_x / (2 * curv_h);
		n_y = K_y / (2 * curv_h);
		n_z = K_z / (2 * curv_h);
		return curv_h;
	}
	else {
		return 0;
	}
}

double MS::vertex::calc_curv_g() {
	/*****************************************************************************
	Must be used after the angles and the area is calculated.

	Use only Vonoroi area even for obtuse triangles.

	Locally, k_g * A = 2\pi - \sum_{angle j around v_i} theta_j
	*****************************************************************************/
	if (USE_VONOROI_CELL) {
		double a = 2 * M_PI, dx_a = 0, dy_a = 0, dz_a = 0;
		std::vector<double> dxn_a(neighbours), dyn_a(neighbours), dzn_a(neighbours);
		for (int i = 0; i < neighbours; i++) {
			a -= theta[i];
			dx_a -= dx_theta[i];
			dy_a -= dy_theta[i];
			dz_a -= dz_theta[i];
			dxn_a[i] -= dxn_theta[i];
			dyn_a[i] -= dyn_theta[i];
			dzn_a[i] -= dzn_theta[i];
			int i_n = neighbour_indices_map[n_next[i]];
			dxn_a[i_n] -= dxnn_theta[i];
			dyn_a[i_n] -= dynn_theta[i];
			dzn_a[i_n] -= dznn_theta[i];
		}
		dx_curv_g = (area*dx_a - a*dx_area) / (area*area);
		dy_curv_g = (area*dy_a - a*dy_area) / (area*area);
		dz_curv_g = (area*dz_a - a*dz_area) / (area*area);
		for (int i = 0; i < neighbours; i++) {
			dxn_curv_g[i] = (area*dxn_a[i] - a*dxn_area[i]) / (area*area);
			dyn_curv_g[i] = (area*dyn_a[i] - a*dyn_area[i]) / (area*area);
			dzn_curv_g[i] = (area*dzn_a[i] - a*dzn_area[i]) / (area*area);
		}
		curv_g = a / area;
		return curv_g;
	}
	else {
		return 0;
	}
}

void MS::vertex::update_geo() {
	calc_angle();
	calc_area();
	calc_curv_h();
	calc_curv_g();
}

void MS::vertex::make_initial() {
	area0 = area;
}

void MS::vertex::make_last() {
	point_last->x = point->x;
	point_last->y = point->y;
	point_last->z = point->z;
}