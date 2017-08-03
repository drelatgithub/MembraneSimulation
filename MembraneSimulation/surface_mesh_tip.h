/*
	This file contains a definition of a point, as would be simulating a filament tip.
	Currently the point has absolute position, and no counter-interaction is imposed on the point.
*/

#pragma once

#include"math_public.h"
#include"surface_mesh.h"

namespace MS {

	struct tip_facet_interaction {
		const facet* which_facet;
		double d;
		math_public::Vec3 nearest_vec;
		int pos; // 0: in triangle; 1, 2, 3: edge 01, 12, 20; 4, 5, 6: vertex 0, 1, 2
		math_public::Vec3 d_d[3];
		math_public::Vec3 dp_d;
	};
	struct tip_edge_interaction {
		double d;
		math_public::Vec3 nearest_vec;
		math_public::Vec3 d_d[2];
		math_public::Vec3 dp_d;
	};
	class filament_tip {
	public:
		math_public::Vec3 *point;
		std::vector<facet*> n_facets; // neighbor facet list

		filament_tip(math_public::Vec3 *np) :point(np) {}

		void get_neighbor_facets(const surface_mesh& sm);
		tip_facet_interaction get_facet_interaction(const facet& f);
		tip_edge_interaction get_edge_interaction(const edge& e);

		std::vector<tip_facet_interaction> interactions;
		tip_facet_interaction* interaction_winner;
		void recalc_interactions();

		/******************************
		Energy part
		******************************/
		double H;
		math_public::Vec3 d_H; // derivative of energy on THIS tip
		void calc_repulsion();

	};
	extern math_public::Vec3 *po;

	extern std::vector<facet*> po_neighbor;
}
