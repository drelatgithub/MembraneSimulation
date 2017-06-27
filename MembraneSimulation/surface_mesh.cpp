#define _USE_MATH_DEFINES

#include"common.h"
#include"surface_mesh.h"

using namespace MS;
using namespace math_public;

#define USE_VONOROI_CELL true


vertex::vertex(Vec3 *npoint) {
	point = npoint;
	point_last = new Vec3(0,0,0);
}

vertex::~vertex() {
	delete point_last;
}

int vertex::count_neighbors() {
	neighbors = n.size();
	for (int i = 0; i < neighbors; i++) {
		neighbor_indices_map[n[i]] = i;
	}
	return neighbors;
}

int vertex::gen_next_prev_n() {
	count_neighbors();
	for (int i = 0; i < neighbors; i++) {
		nn.push_back(n[i < neighbors - 1 ? i + 1 : 0]);
		np.push_back(n[i > 0 ? i - 1 : neighbors - 1]);
	}
	return 0;
}

int vertex::dump_data_vectors() {
	// theta
	theta.push_back(0), sin_theta.push_back(0);
	d_theta.push_back(Vec3()), dn_theta.push_back(Vec3()), dnn_theta.push_back(Vec3());
	d_sin_theta.push_back(Vec3()), dn_sin_theta.push_back(Vec3()), dnn_sin_theta.push_back(Vec3());
	// theta2
	theta2.push_back(0), cot_theta2.push_back(0);
	d_theta2.push_back(Vec3()), dn_theta2.push_back(Vec3()), dnp_theta2.push_back(Vec3());
	d_cot_theta2.push_back(Vec3()), dn_cot_theta2.push_back(Vec3()), dnp_cot_theta2.push_back(Vec3());
	// theta3
	theta3.push_back(0), cot_theta3.push_back(0);
	d_theta3.push_back(Vec3()), dn_theta3.push_back(Vec3()), dnn_theta3.push_back(Vec3());
	d_cot_theta3.push_back(Vec3()), dn_cot_theta3.push_back(Vec3()), dnn_cot_theta3.push_back(Vec3());
	// distances
	r_p_n.push_back(0), d_r_p_n.push_back(Vec3()), dn_r_p_n.push_back(Vec3());
	r_p_np.push_back(0), d_r_p_np.push_back(Vec3()), dnp_r_p_np.push_back(Vec3());
	r_p_nn.push_back(0), d_r_p_nn.push_back(Vec3()), dnn_r_p_nn.push_back(Vec3());

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


void vertex::calc_angle() {
	for (int i = 0; i < neighbors; i++) {
		// Distances
		r_p_n[i] = dist(*point, *(n[i]->point));
		d_r_p_n[i] = (*point - *(n[i]->point)) / r_p_n[i];
		dn_r_p_n[i] = (*(n[i]->point) - *point) / r_p_n[i];

		r_p_np[i] = dist(*point, *(np[i]->point));
		d_r_p_np[i] = (*point - *(np[i]->point)) / r_p_np[i];
		dnp_r_p_np[i] = (*(np[i]->point) - *point) / r_p_np[i];

		r_p_nn[i] = dist(*point, *(nn[i]->point));
		d_r_p_nn[i]= (*point - *(nn[i]->point)) / r_p_nn[i];
		dnn_r_p_nn[i] = (*(nn[i]->point) - *point) / r_p_nn[i];

		// Find theta
		double inner_product = dot(*(nn[i]->point) - *point, *(n[i]->point) - *point);
		Vec3 d_inner_product = 2 * (*point) - *(nn[i]->point) - *(n[i]->point);
		Vec3 dn_inner_product = *(nn[i]->point) - *point;
		Vec3 dnn_inner_product = *(n[i]->point) - *point;

		double cos_theta = inner_product / (r_p_n[i] * r_p_nn[i]);
		Vec3 d_cos_theta = (r_p_n[i] * r_p_nn[i] * d_inner_product - inner_product*(r_p_n[i] * d_r_p_nn[i] + d_r_p_n[i] * r_p_nn[i])) / (r_p_n[i] * r_p_n[i] * r_p_nn[i] * r_p_nn[i]);
		Vec3 dn_cos_theta = (r_p_n[i] * dn_inner_product - inner_product*dn_r_p_n[i]) / (r_p_n[i] * r_p_n[i] * r_p_nn[i]);
		Vec3 dnn_cos_theta = (r_p_nn[i] * dnn_inner_product - inner_product*dnn_r_p_nn[i]) / (r_p_nn[i] * r_p_nn[i] * r_p_n[i]);

		sin_theta[i] = sqrt(1 - cos_theta*cos_theta);
		theta[i] = acos(cos_theta);
		d_theta[i] = -d_cos_theta / sin_theta[i];
		dn_theta[i] = -dn_cos_theta / sin_theta[i];
		dnn_theta[i] = -dnn_cos_theta / sin_theta[i];
		d_sin_theta[i] = cos_theta*d_theta[i];
		dn_sin_theta[i] = cos_theta*dn_theta[i];
		dnn_sin_theta[i] = cos_theta*dnn_theta[i];

		// Find theta2
		double r_n_np = dist(*(n[i]->point), *(np[i]->point)); // derivative with regard to p is 0
		Vec3 dn_r_n_np = (*(n[i]->point) - *(np[i]->point)) / r_n_np;
		Vec3 dnp_r_n_np = (*(np[i]->point) - *(n[i]->point)) / r_n_np;
		double inner_product2 = dot(*(n[i]->point) - *(np[i]->point), *point - *(np[i]->point));
		Vec3 d_inner_product2 = *(n[i]->point) - *(np[i]->point);
		Vec3 dn_inner_product2 = *point - *(np[i]->point);
		Vec3 dnp_inner_product2 = 2 * *(np[i]->point) - *point - *(n[i]->point);

		double cos_theta2 = inner_product2 / (r_p_np[i] * r_n_np);
		Vec3 d_cos_theta2 = (r_p_np[i] * d_inner_product2 - inner_product2*d_r_p_np[i]) / (r_p_np[i] * r_p_np[i] * r_n_np);
		Vec3 dn_cos_theta2 = (r_n_np * dn_inner_product2 - inner_product2*dn_r_n_np) / (r_n_np * r_n_np * r_p_np[i]);
		Vec3 dnp_cos_theta2= (r_p_np[i] * r_n_np * dnp_inner_product2 - inner_product2*(dnp_r_p_np[i] * r_n_np + r_p_np[i] * dnp_r_n_np)) / (r_p_np[i] * r_p_np[i] * r_n_np*r_n_np);

		double sin_theta2 = sqrt(1 - cos_theta2*cos_theta2);
		theta2[i] = acos(cos_theta2);
		d_theta2[i] = -d_cos_theta2 / sin_theta2;
		dn_theta2[i] = -dn_cos_theta2 / sin_theta2;
		dnp_theta2[i] = -dnp_cos_theta2 / sin_theta2;

		cot_theta2[i] = cos_theta2 / sin_theta2;
		d_cot_theta2[i] = -d_theta2[i] / (sin_theta2*sin_theta2);
		dn_cot_theta2[i] = -dn_theta2[i] / (sin_theta2*sin_theta2);
		dnp_cot_theta2[i] = -dnp_theta2[i] / (sin_theta2*sin_theta2);

		// Find theta3
		double r_n_nn = dist(*(n[i]->point), *(nn[i]->point)); // derivative with regard to p is 0
		Vec3 dn_r_n_nn = (*(n[i]->point) - *(nn[i]->point)) / r_n_nn;
		Vec3 dnn_r_n_nn = (*(nn[i]->point) - *(n[i]->point)) / r_n_nn;
		double inner_product3 = dot(*(n[i]->point) - *(nn[i]->point), *point - *(nn[i]->point));
		Vec3 d_inner_product3 = *(n[i]->point) - *(nn[i]->point);
		Vec3 dn_inner_product3 = *point - *(nn[i]->point);
		Vec3 dnn_inner_product3 = 2 * *(nn[i]->point) - *point - *(n[i]->point);

		double cos_theta3 = inner_product3 / (r_p_nn[i] * r_n_nn);
		Vec3 d_cos_theta3 = (r_p_nn[i] * d_inner_product3 - inner_product3*d_r_p_nn[i]) / (r_p_nn[i] * r_p_nn[i] * r_n_nn);
		Vec3 dn_cos_theta3 = (r_n_nn * dn_inner_product3 - inner_product3*dn_r_n_nn) / (r_n_nn * r_n_nn * r_p_nn[i]);
		Vec3 dnn_cos_theta3 = (r_p_nn[i] * r_n_nn * dnn_inner_product3 - inner_product3*(dnn_r_p_nn[i] * r_n_nn + r_p_nn[i] * dnn_r_n_nn)) / (r_p_nn[i] * r_p_nn[i] * r_n_nn*r_n_nn);

		double sin_theta3 = sqrt(1 - cos_theta3*cos_theta3);
		theta3[i] = acos(cos_theta3);
		d_theta3[i] = -d_cos_theta3 / sin_theta3;
		dn_theta3[i] = -dn_cos_theta3 / sin_theta3;
		dnn_theta3[i] = -dnn_cos_theta3 / sin_theta3;

		cot_theta3[i] = cos_theta3 / sin_theta3;
		d_cot_theta3[i] = -d_theta3[i] / (sin_theta3*sin_theta3);
		dn_cot_theta3[i] = -dn_theta3[i] / (sin_theta3*sin_theta3);
		dnn_cot_theta3[i] = -dnn_theta3[i] / (sin_theta3*sin_theta3);

	}
}

double vertex::calc_area() {
	/*****************************************************************************
	Must be used after angles are calculated.

	Use only Vonoroi area even for obtuse triangles.

	The area around vertex i is:
	A = 1/8 * \sum_{v_j is neighbour of v_i} (cot \alpha_ij + cot \beta_ij ) * |r_i - r_j|^2
	*****************************************************************************/
	if (USE_VONOROI_CELL) {
		area = 0;
		dx_area = dy_area = dz_area = 0;
		for (int i = 0; i < neighbors; i++) {
			dxn_area[i] = dyn_area[i] = dzn_area[i] = 0;
		}
		for (int i = 0; i < neighbors; i++) {
			double dis2 = r_p_n[i] * r_p_n[i];
			area += (cot_theta2[i] + cot_theta3[i])*dis2;
			dx_area += (dx_cot_theta2[i] + dx_cot_theta3[i])*dis2 + (cot_theta2[i] + cot_theta3[i]) * 2 * r_p_n[i] * dx_r_p_n[i];
			dy_area += (dy_cot_theta2[i] + dy_cot_theta3[i])*dis2 + (cot_theta2[i] + cot_theta3[i]) * 2 * r_p_n[i] * dy_r_p_n[i];
			dz_area += (dz_cot_theta2[i] + dz_cot_theta3[i])*dis2 + (cot_theta2[i] + cot_theta3[i]) * 2 * r_p_n[i] * dz_r_p_n[i];
			dxn_area[i] += (dxn_cot_theta2[i] + dxn_cot_theta3[i])*dis2 + (cot_theta2[i] + cot_theta3[i]) * 2 * r_p_n[i] * dxn_r_p_n[i];
			dyn_area[i] += (dyn_cot_theta2[i] + dyn_cot_theta3[i])*dis2 + (cot_theta2[i] + cot_theta3[i]) * 2 * r_p_n[i] * dyn_r_p_n[i];
			dzn_area[i] += (dzn_cot_theta2[i] + dzn_cot_theta3[i])*dis2 + (cot_theta2[i] + cot_theta3[i]) * 2 * r_p_n[i] * dzn_r_p_n[i];
			int i_n = neighbor_indices_map.at(n_next[i]), i_p = neighbor_indices_map.at(n_prev[i]);
			dxn_area[i_n] += (dxnn_cot_theta3[i])*dis2;
			dyn_area[i_n] += (dynn_cot_theta3[i])*dis2;
			dzn_area[i_n] += (dznn_cot_theta3[i])*dis2;
			dxn_area[i_p] += (dxnp_cot_theta2[i])*dis2;
			dyn_area[i_p] += (dynp_cot_theta2[i])*dis2;
			dzn_area[i_p] += (dznp_cot_theta2[i])*dis2;
		}
		area /= 8; dx_area /= 8; dy_area /= 8; dz_area /= 8;
		for (int i = 0; i < neighbors; i++) {
			dxn_area[i] /= 8; dyn_area[i] /= 8; dzn_area[i] /= 8;
		}
		return area;
	}
	else {
		return 0;
	}
}

double vertex::calc_curv_h() {
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
		std::vector<double> dxn_K_x(neighbors), dyn_K_x(neighbors), dzn_K_x(neighbors);
		std::vector<double> dxn_K_y(neighbors), dyn_K_y(neighbors), dzn_K_y(neighbors);
		std::vector<double> dxn_K_z(neighbors), dyn_K_z(neighbors), dzn_K_z(neighbors);
		for (int i = 0; i < neighbors; i++) {
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
			int i_n = neighbor_indices_map[n_next[i]], i_p = neighbor_indices_map[n_prev[i]];
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
		for (int i = 0; i < neighbors; i++) {
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
		for (int i = 0; i < neighbors; i++) {
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

double vertex::calc_curv_g() {
	/*****************************************************************************
	Must be used after the angles and the area is calculated.

	Use only Vonoroi area even for obtuse triangles.

	Locally, k_g * A = 2\pi - \sum_{angle j around v_i} theta_j
	*****************************************************************************/
	if (USE_VONOROI_CELL) {
		double a = 2 * M_PI, dx_a = 0, dy_a = 0, dz_a = 0;
		std::vector<double> dxn_a(neighbors), dyn_a(neighbors), dzn_a(neighbors);
		for (int i = 0; i < neighbors; i++) {
			a -= theta[i];
			dx_a -= dx_theta[i];
			dy_a -= dy_theta[i];
			dz_a -= dz_theta[i];
			dxn_a[i] -= dxn_theta[i];
			dyn_a[i] -= dyn_theta[i];
			dzn_a[i] -= dzn_theta[i];
			int i_n = neighbor_indices_map[n_next[i]];
			dxn_a[i_n] -= dxnn_theta[i];
			dyn_a[i_n] -= dynn_theta[i];
			dzn_a[i_n] -= dznn_theta[i];
		}
		dx_curv_g = (area*dx_a - a*dx_area) / (area*area);
		dy_curv_g = (area*dy_a - a*dy_area) / (area*area);
		dz_curv_g = (area*dz_a - a*dz_area) / (area*area);
		for (int i = 0; i < neighbors; i++) {
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

void vertex::update_geo() {
	calc_angle();
	calc_area();
	calc_curv_h();
	//calc_curv_g();
}

void vertex::make_initial() {
	area0 = area;
}

void vertex::make_last() {
	point_last->x = point->x;
	point_last->y = point->y;
	point_last->z = point->z;
}


bool facet::operator==(const facet& operand) {
	int first_index = 0;
	while (first_index < 3) {
		if (v[0] == operand.v[first_index]) { // found a common vertex
			for (int i = 1; i < 3; i++) {
				if (v[i] != operand.v[math_public::loop_add(first_index, i, 3)])
					return false;
			}
			return true;
		}
		++first_index;
	}

	return false;
}