#pragma once

#include<vector>

namespace MS {
	class point_3 {
	public:
		double x, y, z;
		point_3(double nx, double ny, double nz);
	};


	class vertex {
	public:
		point_3 *point;
		std::vector<point_3*> n, n_prev, n_next;

		vertex(point_3 *npoint);
		~vertex();

		// Geometry
		// Local properties other than coordinates may also depend on neighbouring points.
		std::vector<double> theta, dx_theta, dy_theta, dz_theta, sin_theta, dx_sin_theta, dy_sin_theta, dz_sin_theta; // Theta is the angle between p->n and p->n_next
		std::vector<double> theta2, dx_theta2, dy_theta2, dz_theta2, cot_theta2, dx_cot_theta2, dy_cot_theta2, dz_cot_theta2; // Theta2 is the angle between n_prev->p and n_prev->n
		std::vector<double> theta3, dx_theta3, dy_theta3, dz_theta3, cot_theta3, dx_cot_theta3, dy_cot_theta3, dz_cot_theta3; // Theta3 is the angle between n_next->p and n_next->n
		std::vector<double> r_p_n, dx_r_p_n, dy_r_p_n, dz_r_p_n, r_p_n_prev, dx_r_p_n_prev, dy_r_p_n_prev, dz_r_p_n_prev, r_p_n_next, dx_r_p_n_next, dy_r_p_n_next, dz_r_p_n_next;
		int dump_data_vectors();
		double area, dx_area, dy_area, dz_area;
		double curv_h, dx_curv_h, dy_curv_h, dz_curv_h;
		double curv_g, dx_curv_g, dy_curv_g, dz_curv_g;
		double n_x, n_y, n_z; // Components of normalized normal vector
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