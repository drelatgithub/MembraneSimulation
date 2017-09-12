/*
	This file contains a definition of a point, as would be simulating a filament tip.
	Currently the point has absolute position, and no counter-interaction is imposed on the point.
*/

#pragma once

#include"math_public.h"
#include"surface_mesh.h"

namespace MS {

	class filament_tip {
	public:
		math_public::Vec3 *point;
		std::vector<facet*> n_facets; // neighbor facet list

		filament_tip(math_public::Vec3 *np) :point(np) {}

		/******************************
		Energy part
		******************************/
		double H;
		math_public::Vec3 d_H; // derivative of energy on THIS tip. Other derivatives go to vertices.
		void calc_repulsion_facet(facet& f);
		void calc_repulsion(surface_mesh& sm);


		/******************************
		Test
		******************************/
		static test::TestCase test_case;


	};


}
