#define _CRT_SECURE_NO_WARNINGS

#include<iostream>
#include<fstream>

#include "CGAL/Surface_mesh_default_triangulation_3.h"
#include "CGAL/Complex_2_in_triangulation_3.h"
#include "CGAL/make_surface_mesh.h"
#include "CGAL/Implicit_surface_3.h"

#include "CGAL/Polyhedron_3.h"
#include "CGAL/IO/output_surface_facets_to_polyhedron.h"

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

const FT CELL_RADIUS = 1e-6; // 1 micron
const FT CELL_RADIUS_2 = CELL_RADIUS * CELL_RADIUS;

FT sphere_function(Point_3 p) {
	const FT x2 = p.x()*p.x(), y2 = p.y()*p.y(), z2 = p.z()*p.z();
	return x2 + y2 + z2 - CELL_RADIUS_2;
}

int mesh_init() {
	char *position_file = "position.txt";
	char *neighbors_file = "neighbors.txt";
	char *triangles_file = "triangles.txt";

	Tr tr;
	C2t3 c2t3(tr);

	// defining the surface
	Surface_3 surface(sphere_function,             // pointer to function
		Sphere_3(CGAL::ORIGIN, CELL_RADIUS_2 * 1.5)); // bounding sphere
									 // Note that "2." above is the *squared* radius of the bounding sphere!

									 // defining meshing criteria
	CGAL::Surface_mesh_default_criteria_3<Tr> criteria(30.,  // angular bound
		CELL_RADIUS / 20,  // radius bound
		CELL_RADIUS / 20); // distance bound

			  // meshing surface
	std::cout << "Meshing surface using CGAL...\n";
	CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Manifold_tag());
	std::cout << "Number of vertices: " << tr.number_of_vertices() << std::endl;

	// converting to polyhedron
	Polyhedron p;
	std::cout << "Converting mesh into polyhedron...\n";
	CGAL::output_surface_facets_to_polyhedron(c2t3, p);

	// storing vertices info
	std::map<Polyhedron::Vertex_iterator, int> vertices_index_map;
	std::map<Polyhedron::Vertex_iterator, bool> vertices_used;
	std::vector<int> *vertices_neighbors = new std::vector<int>[tr.number_of_vertices()];

	// getting ready for output
	std::ofstream position_out;
	position_out.open(position_file);
	position_out.precision(17); // double
	std::ofstream neighbors_out;
	neighbors_out.open(neighbors_file);
	std::ofstream triangles_out;
	triangles_out.open(triangles_file);

	// registering vertex positions
	int num = 0;
	for (Polyhedron::Vertex_iterator vit = p.vertices_begin(); vit != p.vertices_end(); vit++) {
		position_out << std::scientific << vit->point().x() << '\t' << std::scientific << vit->point().y() << '\t' << std::scientific << vit->point().z() << std::endl;
		vertices_used[vit] = false;
		vertices_index_map[vit] = num;
		num++;
	} // num now becomes the number of vertices after loop

	// registering neighbours
	Polyhedron::Vertex_iterator vit;
	Polyhedron::Halfedge_iterator hit_new;
	for (Polyhedron::Halfedge_iterator hit = p.halfedges_begin(); hit != p.halfedges_end(); hit++) {
		vit = hit->vertex();
		if (!vertices_used[vit]) {
			vertices_used[vit] = true;
			hit_new = hit;
			do {
				vertices_neighbors[vertices_index_map[vit]].push_back(vertices_index_map[hit_new->opposite()->vertex()]);
				hit_new = hit_new->opposite()->prev();
			} while (hit_new != hit);
		}
	}
	for (int i = 0; i < num; i++) {
		for (unsigned int j = 0; j < vertices_neighbors[i].size(); j++) {
			neighbors_out << vertices_neighbors[i][j] << '\t';
		}
		neighbors_out << std::endl;
	}
	delete[] vertices_neighbors;

	// registering triangles
	for (Polyhedron::Facet_iterator fit = p.facets_begin(); fit != p.facets_end(); fit++) {
		Polyhedron::Halfedge_iterator hit = fit->facet_begin();
		Polyhedron::Halfedge_iterator hit_new = hit;
		do {
			triangles_out << vertices_index_map[hit_new->vertex()] << '\t';
			hit_new = hit_new->next();
		} while (hit_new != hit);
		triangles_out << std::endl;
	}

	std::cout << "Finished writing files.";
	position_out.close();
	neighbors_out.close();
	triangles_out.close();

	return 0;
}

int main() {
	mesh_init();
	return 0;
}