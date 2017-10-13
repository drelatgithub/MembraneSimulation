#pragma once

#include"glew/glew.h"

class gl_vertex_holder {
public:
	std::vector<GLfloat> vertex_buffer_data;
	size_t num_tris;

	GLfloat* gen_mesh_vertex_buffer_data(const char* pos_out_path, const char* tri_path, size_t snapshot);

};
