#pragma once

/**********************************************************

Definition of point and vertex structure, and some common functions associated.

**********************************************************/

#include<vector>
#include<map>

namespace MS {
	class point_3 {
	public:
		double x, y, z;
		point_3(double nx, double ny, double nz);
	};

	double distance2(const point_3 *p, const point_3 *np);
	double distance(const point_3 *p, const point_3 *np);


	class vertex {
	public:
		point_3 *point;
		std::vector<vertex*> n, n_prev, n_next;

		vertex(point_3 *npoint);
		~vertex();

		int neighbours;
		std::map<vertex*, int> neighbour_indices_map;
		int count_neighbours();

		// Geometry
		// Local properties other than coordinates may also depend on neighbouring points.
		// Theta is the angle between p->n and p->n_next
		std::vector<double> theta, dx_theta, dy_theta, dz_theta, sin_theta, dx_sin_theta, dy_sin_theta, dz_sin_theta,
			dxn_theta, dyn_theta, dzn_theta, dxn_sin_theta, dyn_sin_theta, dzn_sin_theta,
			dxnn_theta, dynn_theta, dznn_theta, dxnn_sin_theta, dynn_sin_theta, dznn_sin_theta;
		// Theta2 is the angle between n_prev->p and n_prev->n
		std::vector<double> theta2, dx_theta2, dy_theta2, dz_theta2, cot_theta2, dx_cot_theta2, dy_cot_theta2, dz_cot_theta2,
			dxn_theta2, dyn_theta2, dzn_theta2, dxn_cot_theta2, dyn_cot_theta2, dzn_cot_theta2,
			dxnp_theta2, dynp_theta2, dznp_theta2, dxnp_cot_theta2, dynp_cot_theta2, dznp_cot_theta2;
		// Theta3 is the angle between n_next->p and n_next->n
		std::vector<double> theta3, dx_theta3, dy_theta3, dz_theta3, cot_theta3, dx_cot_theta3, dy_cot_theta3, dz_cot_theta3,
			dxn_theta3, dyn_theta3, dzn_theta3, dxn_cot_theta3, dyn_cot_theta3, dzn_cot_theta3,
			dxnn_theta3, dynn_theta3, dznn_theta3, dxnn_cot_theta3, dynn_cot_theta3, dznn_cot_theta3;
		// Distances
		std::vector<double> r_p_n, dx_r_p_n, dy_r_p_n, dz_r_p_n, dxn_r_p_n, dyn_r_p_n, dzn_r_p_n;
		std::vector<double> r_p_n_prev, dx_r_p_n_prev, dy_r_p_n_prev, dz_r_p_n_prev, dxnp_r_p_n_prev, dynp_r_p_n_prev, dznp_r_p_n_prev;
		std::vector<double> r_p_n_next, dx_r_p_n_next, dy_r_p_n_next, dz_r_p_n_next, dxnn_r_p_n_next, dynn_r_p_n_next, dznn_r_p_n_next;

		double area, dx_area, dy_area, dz_area;
		std::vector<double> dxn_area, dyn_area, dzn_area;
		double curv_h, dx_curv_h, dy_curv_h, dz_curv_h;
		std::vector<double> dxn_curv_h, dyn_curv_h, dzn_curv_h;
		double curv_g, dx_curv_g, dy_curv_g, dz_curv_g;
		std::vector<double>dxn_curv_g, dyn_curv_g, dzn_curv_g;
		double n_x, n_y, n_z; // Components of normalized normal vector

		int dump_data_vectors();

		void calc_angle();
		double calc_area();
		double calc_curv_h();
		double calc_curv_g();

		void update_geo();

		double area0;
		point_3 *point_last;
		void make_initial(); // Making the current geometry the initial geometry
		void make_last(); // Recording some of the geometry as the last time geometry
	};

}