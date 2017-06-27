#pragma once

/**********************************************************

Definition of point and vertex structure, and some common functions associated.

**********************************************************/

#include<vector>
#include<map>

#include"common.h"
#include"math_public.h"

namespace MS {

	class vertex {
	public:
		math_public::Vec3 *point;
		std::vector<vertex*> n, np, nn; // neighbor, previous neighbor, next neighbor

		vertex(math_public::Vec3 *npoint);
		~vertex();

		int neighbors;
		std::map<vertex*, int> neighbor_indices_map;
		int count_neighbors();
		int gen_next_prev_n();

		// Geometry
		// Local properties other than coordinates may also depend on neighbouring points.
		// Theta is the angle between p->n and p->nn
		std::vector<double> theta, sin_theta;
		std::vector<math_public::Vec3>d_theta, d_sin_theta, dn_theta, dn_sin_theta, dnn_theta, dnn_sin_theta;
		// Theta2 is the angle between np->p and np->n
		std::vector<double> theta2, dx_theta2, dy_theta2, dz_theta2, cot_theta2, dx_cot_theta2, dy_cot_theta2, dz_cot_theta2,
			dxn_theta2, dyn_theta2, dzn_theta2, dxn_cot_theta2, dyn_cot_theta2, dzn_cot_theta2,
			dxnp_theta2, dynp_theta2, dznp_theta2, dxnp_cot_theta2, dynp_cot_theta2, dznp_cot_theta2;
		std::vector<math_public::Vec3>d_theta2, d_cot_theta2, dn_theta2, dn_cot_theta2, dnp_theta2, dnp_cot_theta2;
		// Theta3 is the angle between nn->p and nn->n
		std::vector<double> theta3, dx_theta3, dy_theta3, dz_theta3, cot_theta3, dx_cot_theta3, dy_cot_theta3, dz_cot_theta3,
			dxn_theta3, dyn_theta3, dzn_theta3, dxn_cot_theta3, dyn_cot_theta3, dzn_cot_theta3,
			dxnn_theta3, dynn_theta3, dznn_theta3, dxnn_cot_theta3, dynn_cot_theta3, dznn_cot_theta3;
		std::vector<math_public::Vec3>d_theta3, d_cot_theta3, dn_theta3, dn_cot_theta3, dnn_theta3, dnn_cot_theta3;
		// Distances
		std::vector<double> r_p_n;
		std::vector<math_public::Vec3> d_r_p_n, dn_r_p_n;
		std::vector<double> r_p_np;
		std::vector<math_public::Vec3> d_r_p_np, dnp_r_p_np;
		std::vector<double> r_p_nn;
		std::vector<math_public::Vec3> d_r_p_nn, dnn_r_p_nn;

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
		math_public::Vec3 *point_last;
		void make_initial(); // Making the current geometry the initial geometry
		void make_last(); // Recording some of the geometry as the last time geometry
	};

	// A facet stores pointers to 3 vertices which form a triangle in the counter-clockwise direction
	class facet {
	public:
		vertex *v[3];

		facet(vertex *v0, vertex *v1, vertex *v2) { v[0] = v0; v[1] = v1; v[2] = v2; }
		bool operator==(const facet& operand);
	};

}