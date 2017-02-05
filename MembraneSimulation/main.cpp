#define _CRT_SECURE_NO_WARNINGS

#include<iostream>

#include "CGAL/Surface_mesh_default_triangulation_3.h"
#include "CGAL/Complex_2_in_triangulation_3.h"
#include "CGAL/make_surface_mesh.h"
#include "CGAL/Implicit_surface_3.h"

#include "CGAL/Polyhedron_3.h"
#include "CGAL/IO/output_surface_facets_to_polyhedron.h"

#include"surface_mesh.h"
#include"simulation_process.h"

// default triangulation for Surface_mesher
typedef CGAL::Surface_mesh_default_triangulation_3 Tr;
// c2t3
typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;

typedef Tr::Geom_traits GT;
typedef GT::Sphere_3 Sphere_3;
typedef GT::Point_3 Point_3;
typedef GT::FT FT;

typedef FT(*Function)(Point_3);

typedef CGAL::Implicit_surface_3<GT, Function> Surface_3;

typedef CGAL::Polyhedron_3<GT> Polyhedron;

FT sphere_function(Point_3 p) {
	const FT x2 = p.x()*p.x(), y2 = p.y()*p.y(), z2 = p.z()*p.z();
	return x2 + y2 + z2 - 2;
}

int main() {
	Tr tr;
	C2t3 c2t3(tr);

	// defining the surface
	Surface_3 surface(sphere_function,             // pointer to function
		Sphere_3(CGAL::ORIGIN, 4.)); // bounding sphere
									 // Note that "2." above is the *squared* radius of the bounding sphere!

	// defining meshing criteria
	CGAL::Surface_mesh_default_criteria_3<Tr> criteria(30.,  // angular bound
														0.1,  // radius bound
														0.1); // distance bound

	// meshing surface
	std::cout << "Meshing surface using CGAL...\n";
	CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Manifold_tag());
	std::cout << "Number of vertices: " << tr.number_of_vertices() << std::endl;

	// converting to polyhedron
	Polyhedron p;
	std::cout << "Converting mesh into polyhedron...\n";
	CGAL::output_surface_facets_to_polyhedron(c2t3, p);

	// storing vertices info
	std::vector<MS::vertex> vertices;
	vertices.reserve(tr.number_of_vertices());
	std::map<Polyhedron::Vertex_iterator, int> vertices_index_map;
	//std::map<MS::point_3*, int> vertices_index_map_point;

	int num = 0;
	for (Polyhedron::Vertex_iterator vit = p.vertices_begin(); vit != p.vertices_end(); vit++) {
		MS::point_3 *pt = new MS::point_3(vit->point().x(), vit->point().y(), vit->point().z());
		MS::vertex *new_vertex = new MS::vertex(pt);
		vertices.push_back(*new_vertex);
		//vertices_index_map_point[pt] = num;
		vertices_index_map[vit] = num++;
	}
	// registering neighbours
	for (Polyhedron::Halfedge_iterator hit = p.halfedges_begin(); hit != p.halfedges_end(); hit++) {
		int this_index = vertices_index_map[hit->vertex()];
		vertices[this_index].n.push_back(vertices[vertices_index_map[hit->opposite()->vertex()]].point);
		vertices[this_index].n_prev.push_back(vertices[vertices_index_map[hit->opposite()->next()->vertex()]].point);
		vertices[this_index].n_next.push_back(vertices[vertices_index_map[hit->next()->vertex()]].point);
		//MS::point_3 * t = vertices[vertices_index_map[hit->next()->vertex()]].point;
		//std::cout << t->x << '\t';
		vertices[this_index].dump_data_vectors();
	}

	// starting simulation
	std::cout << "Simulation starting...\n";
	MS::simulation_start(vertices);
	


	system("pause");

	return 0;
}