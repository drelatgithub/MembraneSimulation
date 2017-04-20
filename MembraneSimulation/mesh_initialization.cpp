#define _CRT_SECURE_NO_WARNINGS

#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<sys/stat.h>

#include "CGAL/Surface_mesh_default_triangulation_3.h"
#include "CGAL/Complex_2_in_triangulation_3.h"
#include "CGAL/make_surface_mesh.h"
#include "CGAL/Implicit_surface_3.h"

#include "CGAL/Polyhedron_3.h"
#include "CGAL/IO/output_surface_facets_to_polyhedron.h"

#include"surface_mesh.h"
#include"simulation_process.h"

#define USE_SAVED_MESH true

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

int mesh_init(std::vector<MS::vertex*> &vertices) {
	char *data_file = "v_data.txt";

	struct stat buffer;
	if (USE_SAVED_MESH && stat(data_file, &buffer) == 0) { // File exists
		std::cout << "Saved mesh file found. Trying to build from file...\n";

		std::ifstream v_file_in;
		v_file_in.open(data_file);

		std::string line;
		while (std::getline(v_file_in, line)) {
			std::stringstream ss(line);
			double x, y, z;
			if (ss >> x >> y >> z) {
				MS::point_3 *pt = new MS::point_3(x, y, z);
				MS::vertex *new_vertex = new MS::vertex(pt);
				vertices.push_back(new_vertex);
			}
		}

		v_file_in.clear();
		v_file_in.seekg(0, std::ios::beg);

		int num = 0;
		while (std::getline(v_file_in, line)) {
			std::stringstream ss(line);
			double trash;
			if (!(ss >> trash >> trash >> trash))break;
			int n, np, nn;
			while (ss >> n >> np >> nn) {
				vertices[num]->n.push_back(vertices[n]);
				vertices[num]->n_prev.push_back(vertices[np]);
				vertices[num]->n_next.push_back(vertices[nn]);
				vertices[num]->dump_data_vectors();
			}
			num++;
		}

		std::cout << "Number of vertices: " << num << std::endl;
		v_file_in.close();

		return 0;
	}

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
	vertices.reserve(tr.number_of_vertices());
	std::map<Polyhedron::Vertex_iterator, int> vertices_index_map;
	std::map<MS::vertex*, int> vertices_index_map_new;

	int num = 0;
	for (Polyhedron::Vertex_iterator vit = p.vertices_begin(); vit != p.vertices_end(); vit++) {
		MS::point_3 *pt = new MS::point_3(vit->point().x(), vit->point().y(), vit->point().z());
		MS::vertex *new_vertex = new MS::vertex(pt);
		vertices.push_back(new_vertex);
		vertices_index_map[vit] = num;
		vertices_index_map_new[vertices[num]] = num;
		num++;
	} // num now becomes the number of vertices

	// registering neighbours
	std::ofstream v_file_out;
	v_file_out.open(data_file);
	v_file_out.precision(17); // double

	for (Polyhedron::Halfedge_iterator hit = p.halfedges_begin(); hit != p.halfedges_end(); hit++) {
		int this_index = vertices_index_map[hit->vertex()];
		vertices[this_index]->n.push_back(vertices[vertices_index_map[hit->opposite()->vertex()]]);
		vertices[this_index]->n_prev.push_back(vertices[vertices_index_map[hit->opposite()->next()->vertex()]]);
		vertices[this_index]->n_next.push_back(vertices[vertices_index_map[hit->next()->vertex()]]);
		vertices[this_index]->dump_data_vectors();
	}
	for (int i = 0; i < num; i++) {
		v_file_out << std::scientific << vertices[i]->point->x << '\t' << std::scientific << vertices[i]->point->y << '\t' << std::scientific << vertices[i]->point->z << '\t';
		vertices[i]->count_neighbours();
		for (int j = 0; j < vertices[i]->neighbours; j++) {
			v_file_out << vertices_index_map_new[vertices[i]->n[j]] << '\t' << vertices_index_map_new[vertices[i]->n_prev[j]] << '\t' << vertices_index_map_new[vertices[i]->n_next[j]] << '\t';
		}
		v_file_out << std::endl;
	}
	v_file_out.close();

	return 0;
}