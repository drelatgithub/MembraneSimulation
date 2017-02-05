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

int MS::vertex::dump_data_vectors() {
	theta.push_back(0), dx_theta.push_back(0), dy_theta.push_back(0), dz_theta.push_back(0);
	sin_theta.push_back(0), dx_sin_theta.push_back(0), dy_sin_theta.push_back(0), dz_sin_theta.push_back(0);
	theta2.push_back(0), dx_theta2.push_back(0), dy_theta2.push_back(0), dz_theta2.push_back(0);
	cot_theta2.push_back(0), dx_cot_theta2.push_back(0), dy_cot_theta2.push_back(0), dz_cot_theta2.push_back(0);
	theta3.push_back(0), dx_theta3.push_back(0), dy_theta3.push_back(0), dz_theta3.push_back(0);
	cot_theta3.push_back(0), dx_cot_theta3.push_back(0), dy_cot_theta3.push_back(0), dz_cot_theta3.push_back(0);
	r_p_n.push_back(0), dx_r_p_n.push_back(0), dy_r_p_n.push_back(0), dz_r_p_n.push_back(0);
	r_p_n_prev.push_back(0), dx_r_p_n_prev.push_back(0), dy_r_p_n_prev.push_back(0), dz_r_p_n_prev.push_back(0);
	r_p_n_next.push_back(0), dx_r_p_n_next.push_back(0), dy_r_p_n_next.push_back(0), dz_r_p_n_next.push_back(0);
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
	for (int i = 0, len = n.size(); i < len; i++) {
		// Distances
		r_p_n[i] = distance(point, n[i]);
		dx_r_p_n[i] = (point->x - n[i]->x) / r_p_n[i];
		dy_r_p_n[i] = (point->y - n[i]->y) / r_p_n[i];
		dz_r_p_n[i] = (point->z - n[i]->z) / r_p_n[i];
		r_p_n_prev[i] = distance(point, n_prev[i]);
		dx_r_p_n_prev[i] = (point->x - n_prev[i]->x) / r_p_n_prev[i];
		dy_r_p_n_prev[i] = (point->y - n_prev[i]->y) / r_p_n_prev[i];
		dz_r_p_n_prev[i] = (point->z - n_prev[i]->z) / r_p_n_prev[i];
		r_p_n_next[i] = distance(point, n_next[i]);
		dx_r_p_n_next[i] = (point->x - n_next[i]->x) / r_p_n_next[i];
		dy_r_p_n_next[i] = (point->y - n_next[i]->y) / r_p_n_next[i];
		dz_r_p_n_next[i] = (point->z - n_next[i]->z) / r_p_n_next[i];
		//std::cout << r_p_n[i];

		// Find theta
		double inner_product = (n_next[i]->x - point->x)*(n[i]->x - point->x) + (n_next[i]->y - point->y)*(n[i]->y - point->y) + (n_next[i]->z - point->z)*(n[i]->z - point->z);
		double dx_inner_product = 2 * point->x - n_next[i]->x - n[i]->x;
		double dy_inner_product = 2 * point->y - n_next[i]->y - n[i]->y;
		double dz_inner_product = 2 * point->z - n_next[i]->z - n[i]->z;
		double cos_theta = inner_product / (r_p_n[i] * r_p_n_next[i]);
		double dx_cos_theta = (r_p_n[i] * r_p_n_next[i] * dx_inner_product - inner_product*(r_p_n[i] * dx_r_p_n_next[i] + dx_r_p_n[i] * r_p_n_next[i])) / (r_p_n[i] * r_p_n[i] * r_p_n_next[i] * r_p_n_next[i]);
		double dy_cos_theta = (r_p_n[i] * r_p_n_next[i] * dy_inner_product - inner_product*(r_p_n[i] * dy_r_p_n_next[i] + dy_r_p_n[i] * r_p_n_next[i])) / (r_p_n[i] * r_p_n[i] * r_p_n_next[i] * r_p_n_next[i]);
		double dz_cos_theta = (r_p_n[i] * r_p_n_next[i] * dz_inner_product - inner_product*(r_p_n[i] * dz_r_p_n_next[i] + dz_r_p_n[i] * r_p_n_next[i])) / (r_p_n[i] * r_p_n[i] * r_p_n_next[i] * r_p_n_next[i]);
		sin_theta[i] = sqrt(1 - cos_theta*cos_theta);
		theta[i] = acos(cos_theta);
		dx_theta[i] = -dx_cos_theta / sin_theta[i];
		dy_theta[i] = -dy_cos_theta / sin_theta[i];
		dz_theta[i] = -dz_cos_theta / sin_theta[i];
		dx_sin_theta[i] = cos_theta*dx_theta[i];
		dy_sin_theta[i] = cos_theta*dy_theta[i];
		dz_sin_theta[i] = cos_theta*dz_theta[i];
		//std::cout << theta[i] << '\t';

		// Find theta2
		double r_n_n_prev = distance(n[i], n_prev[i]); // derivative is 0
		double inner_product2 = (n[i]->x - n_prev[i]->x)*(point->x - n_prev[i]->x) + (n[i]->y - n_prev[i]->y)*(point->y - n_prev[i]->y) + (n[i]->z - n_prev[i]->z)*(point->z - n_prev[i]->z);
		double dx_inner_product2 = n[i]->x - n_prev[i]->x;
		double dy_inner_product2 = n[i]->y - n_prev[i]->y;
		double dz_inner_product2 = n[i]->z - n_prev[i]->z;
		double cos_theta2 = inner_product2 / (r_p_n_prev[i] * r_n_n_prev);
		double dx_cos_theta2 = (r_p_n_prev[i] * r_n_n_prev * dx_inner_product2 - inner_product2*(dx_r_p_n_prev[i] * r_n_n_prev)) / (r_p_n_prev[i] * r_p_n_prev[i] * r_n_n_prev * r_n_n_prev);
		double dy_cos_theta2 = (r_p_n_prev[i] * r_n_n_prev * dy_inner_product2 - inner_product2*(dy_r_p_n_prev[i] * r_n_n_prev)) / (r_p_n_prev[i] * r_p_n_prev[i] * r_n_n_prev * r_n_n_prev);
		double dz_cos_theta2 = (r_p_n_prev[i] * r_n_n_prev * dz_inner_product2 - inner_product2*(dz_r_p_n_prev[i] * r_n_n_prev)) / (r_p_n_prev[i] * r_p_n_prev[i] * r_n_n_prev * r_n_n_prev);
		double sin_theta2 = sqrt(1 - cos_theta2*cos_theta2);
		theta2[i] = acos(cos_theta2);
		dx_theta2[i] = -dx_cos_theta2 / sin_theta2;
		dy_theta2[i] = -dy_cos_theta2 / sin_theta2;
		dz_theta2[i] = -dz_cos_theta2 / sin_theta2;
		cot_theta2[i] = cos_theta2 / sin_theta2;
		dx_cot_theta2[i] = dx_theta2[i] / (sin_theta2*sin_theta2);
		dy_cot_theta2[i] = dy_theta2[i] / (sin_theta2*sin_theta2);
		dz_cot_theta2[i] = dz_theta2[i] / (sin_theta2*sin_theta2);

		// Find theta3
		double r_n_n_next = distance(n[i], n_next[i]); // derivative is 0
		double inner_product3 = (n[i]->x - n_next[i]->x)*(point->x - n_next[i]->x) + (n[i]->y - n_next[i]->y)*(point->y - n_next[i]->y) + (n[i]->z - n_next[i]->z)*(point->z - n_next[i]->z);
		double dx_inner_product3 = n[i]->x - n_next[i]->x;
		double dy_inner_product3 = n[i]->y - n_next[i]->y;
		double dz_inner_product3 = n[i]->z - n_next[i]->z;
		double cos_theta3 = inner_product3 / (r_p_n_next[i] * r_n_n_next);
		double dx_cos_theta3 = (r_p_n_next[i] * r_n_n_next * dx_inner_product3 - inner_product3*(dx_r_p_n_next[i] * r_n_n_next)) / (r_p_n_next[i] * r_p_n_next[i] * r_n_n_next * r_n_n_next);
		double dy_cos_theta3 = (r_p_n_next[i] * r_n_n_next * dy_inner_product3 - inner_product3*(dy_r_p_n_next[i] * r_n_n_next)) / (r_p_n_next[i] * r_p_n_next[i] * r_n_n_next * r_n_n_next);
		double dz_cos_theta3 = (r_p_n_next[i] * r_n_n_next * dz_inner_product3 - inner_product3*(dz_r_p_n_next[i] * r_n_n_next)) / (r_p_n_next[i] * r_p_n_next[i] * r_n_n_next * r_n_n_next);
		double sin_theta3 = sqrt(1 - cos_theta3*cos_theta3);
		theta3[i] = acos(cos_theta3);
		dx_theta3[i] = -dx_cos_theta3 / sin_theta3;
		dy_theta3[i] = -dy_cos_theta3 / sin_theta3;
		dz_theta3[i] = -dz_cos_theta3 / sin_theta3;
		cot_theta3[i] = cos_theta3 / sin_theta3;
		dx_cot_theta3[i] = dx_theta3[i] / (sin_theta3*sin_theta3);
		dy_cot_theta3[i] = dy_theta3[i] / (sin_theta3*sin_theta3);
		dz_cot_theta3[i] = dz_theta3[i] / (sin_theta3*sin_theta3);
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
		for (int i = 0, len = n.size(); i < len; i++) {
			area += (cot_theta2[i] + cot_theta3[i])*r_p_n[i] * r_p_n[i];
			dx_area += (dx_cot_theta2[i] + dx_cot_theta3[i])*r_p_n[i] * r_p_n[i] + (cot_theta2[i] + cot_theta3[i]) * 2 * r_p_n[i] * dx_r_p_n[i];
			dy_area += (dy_cot_theta2[i] + dy_cot_theta3[i])*r_p_n[i] * r_p_n[i] + (cot_theta2[i] + cot_theta3[i]) * 2 * r_p_n[i] * dy_r_p_n[i];
			dz_area += (dz_cot_theta2[i] + dz_cot_theta3[i])*r_p_n[i] * r_p_n[i] + (cot_theta2[i] + cot_theta3[i]) * 2 * r_p_n[i] * dz_r_p_n[i];
		}
		area /= 8; dx_area /= 8; dy_area /= 8; dz_area /= 8;
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
		for (int i = 0, len = n.size(); i < len; i++) {
			K_x += (cot_theta2[i] + cot_theta3[i])*(point->x - n[i]->x);
			K_y += (cot_theta2[i] + cot_theta3[i])*(point->y - n[i]->y);
			K_z += (cot_theta2[i] + cot_theta3[i])*(point->z - n[i]->z);
			dx_K_x += (dx_cot_theta2[i] + dx_cot_theta3[i])*(point->x - n[i]->x) + (cot_theta2[i] + cot_theta3[i]);
			dy_K_x += (dy_cot_theta2[i] + dy_cot_theta3[i])*(point->x - n[i]->x);
			dz_K_x += (dz_cot_theta2[i] + dz_cot_theta3[i])*(point->x - n[i]->x);
			dx_K_y += (dx_cot_theta2[i] + dx_cot_theta3[i])*(point->y - n[i]->y);
			dy_K_y += (dy_cot_theta2[i] + dy_cot_theta3[i])*(point->y - n[i]->y) + (cot_theta2[i] + cot_theta3[i]);
			dz_K_y += (dz_cot_theta2[i] + dz_cot_theta3[i])*(point->y - n[i]->y);
			dx_K_z += (dx_cot_theta2[i] + dx_cot_theta3[i])*(point->z - n[i]->z);
			dy_K_z += (dy_cot_theta2[i] + dy_cot_theta3[i])*(point->z - n[i]->z);
			dz_K_z += (dz_cot_theta2[i] + dz_cot_theta3[i])*(point->z - n[i]->z) + (cot_theta2[i] + cot_theta3[i]);
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
		K_x /= 2 * area; K_y /= 2 * area; K_z /= 2 * area;
		curv_h = sqrt(K_x*K_x + K_y*K_y + K_z*K_z) / 2;
		dx_curv_h = (K_x*dx_K_x + K_y*dx_K_y + K_z*dx_K_z) / (4 * curv_h);
		dy_curv_h = (K_x*dy_K_x + K_y*dy_K_y + K_z*dy_K_z) / (4 * curv_h);
		dx_curv_h = (K_x*dz_K_x + K_y*dz_K_y + K_z*dz_K_z) / (4 * curv_h);
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
		for (int i = 0, len = n.size(); i < len; i++) {
			a -= theta[i];
			dx_a += dx_theta[i];
			dy_a += dy_theta[i];
			dz_a += dz_theta[i];
		}
		dx_curv_g = (area*dx_a - a*dx_area) / (area*area);
		dy_curv_g = (area*dy_a - a*dy_area) / (area*area);
		dz_curv_g = (area*dz_a - a*dz_area) / (area*area);
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
	area0 = area-0.0001;
}

void MS::vertex::make_last() {
	point_last->x = point->x;
	point_last->y = point->y;
	point_last->z = point->z;
}