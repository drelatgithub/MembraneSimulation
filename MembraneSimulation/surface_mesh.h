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
		inline void release_point() { delete point; }

		int neighbors;
		std::map<vertex*, int> neighbor_indices_map;
		int count_neighbors();
		int gen_next_prev_n();

		/******************************
		Geometry part
		******************************/
		// Local properties other than coordinates may also depend on neighbouring points.
		// Theta is the angle between p->n and p->nn
		std::vector<double> theta, sin_theta;
		std::vector<math_public::Vec3>d_theta, d_sin_theta, dn_theta, dn_sin_theta, dnn_theta, dnn_sin_theta;
		// Theta2 is the angle between np->p and np->n
		std::vector<double> theta2, cot_theta2;
		std::vector<math_public::Vec3>d_theta2, d_cot_theta2, dn_theta2, dn_cot_theta2, dnp_theta2, dnp_cot_theta2;
		// Theta3 is the angle between nn->p and nn->n
		std::vector<double> theta3, cot_theta3;
		std::vector<math_public::Vec3>d_theta3, d_cot_theta3, dn_theta3, dn_cot_theta3, dnn_theta3, dnn_cot_theta3;
		// Distances
		std::vector<double> r_p_n;
		std::vector<math_public::Vec3> d_r_p_n, dn_r_p_n;
		std::vector<double> r_p_np;
		std::vector<math_public::Vec3> d_r_p_np, dnp_r_p_np;
		std::vector<double> r_p_nn;
		std::vector<math_public::Vec3> d_r_p_nn, dnn_r_p_nn;

		double area;
		math_public::Vec3 d_area;
		std::vector<math_public::Vec3> dn_area;
		double curv_h;
		math_public::Vec3 d_curv_h;
		std::vector<math_public::Vec3> dn_curv_h;
		double curv_g;
		math_public::Vec3 d_curv_g;
		std::vector<math_public::Vec3> dn_curv_g;
		
		math_public::Vec3 n_vec; // Normal vector
		// we do not calculate derivative of normal vector here, because it would be very inaccurate when curv_h is close to 0.
		//math_public::Mat3 d_n_vec;
		//std::vector<math_public::Mat3> dn_n_vec;

		int dump_data_vectors(int size=1);

		void calc_angle();
		double calc_area();
		double calc_curv_h();
		double calc_curv_g();

		void update_geo();

		double area0;
		math_public::Vec3 *point_last;
		void make_initial(); // Making the current geometry the initial geometry
		void make_last(); // Recording some of the geometry as the last time geometry


		/******************************
		Energy part
		******************************/
		double H_area, // Energy due to change in area
			H_curv_h, // Energy due to mean curvature
			H_curv_g, // Energy due to Gaussian curvature
			H_int; // Energy due to interaction with objects that are not neighboring vertices.
		double H; // Free energy of this vertex in J. Interaction energy with facet lies in the facet class.

		math_public::Vec3 d_H_area,
			d_H_curv_h,
			d_H_curv_g,
			d_H_int;
		math_public::Vec3 d_H; // Energy derivative on this vertex. in J/m.

		void clear_energy();

		// Calculate energy and derivatives
		void calc_H_area();
		void calc_H_curv_h();
		void calc_H_curv_g();
		inline void calc_H_int() { H_int = 0; d_H_int.set(0, 0, 0); } // This actually serves as cleaning
		void inc_d_H_int(const math_public::Vec3 &d);
		
		inline void sum_energy() {
			/*
				Energy caused by Gaussian curvature could be neglected because for a closed surface
				it is a constant (Gauss-Bonnet theorem).
			*/
			// some of d_H_int might come from facet energies
			H = H_area + H_curv_h + H_int;
			d_H = d_H_area + d_H_curv_h + d_H_int;
		}

		void update_energy();

		/******************************
		Test
		******************************/
		static test::TestCase test_case;

	};

	// A facet stores pointers to 3 vertices which form a triangle in the counter-clockwise direction
	class facet {
	public:
		vertex *v[3];

		/******************************
		Geometry is mostly calculated in vertex class
		******************************/
		int ind[3]; // the neighbor indices of v0-v1, v1-v2, v2-v0.

		facet(vertex *v0, vertex *v1, vertex *v2) { v[0] = v0; v[1] = v1; v[2] = v2; }
		bool operator==(const facet& operand);

		void update_geo();

		/******************************
		Energy part
		******************************/
		double H_int;
		double H;
		// Energy derivative would be on vertices

		inline void sum_energy() {
			H = H_int;
		}

		inline void calc_H_int() { H_int = 0; } // Serves as cleaning
		void inc_H_int(math_public::Vec3 *p);

		void update_energy(math_public::Vec3 *p); // TODO: This should be changed to interacting with all points

		/******************************
		Test
		******************************/
		static test::TestCase test_case;

	};

	class surface_mesh {
	public:
		std::vector<vertex*> vertices;
		std::vector<facet*> facets;
	};

}