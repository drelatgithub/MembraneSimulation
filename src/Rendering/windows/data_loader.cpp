#include"common.h"
#include"data_loader.h"

const GLfloat* gen_mesh_vertex_buffer_data(const char* pos_out_path, const char* tri_path, size_t snapshot) {
	bool cont = true;

	size_t N = 0;
	std::vector<double> coords;

	std::string vertex_pos_content;
	std::ifstream vertex_pos(pos_out_path, std::ios::in);
	if (vertex_pos.is_open()) {
		// read line by line until the desired snapshot
		// For example, if snapshot == 3, then 3 lines are skipped and try to read the 4th line
		for (int i = 0; i < snapshot && std::getline(vertex_pos, vertex_pos_content); i++);
		if (vertex_pos && std::getline(vertex_pos, vertex_pos_content)) {
			// Now the expected line is stored, let's parse them
			double each_num;
			std::istringstream iss(vertex_pos_content, std::istringstream::in);
			coords.reserve(vertex_pos_content.length() / 12); // Should be slightly bigger than needed.

			while (iss >> each_num) {
				coords.push_back(each_num);
			}

			if (coords.size() % 3) {
				LOG(ERROR) << "Number of coordinates must be multiples of 3.";
				cont = false;
			}
			else {
				N = coords.size() / 3;
			}
		}
		else {
			LOG(ERROR) << "File ended before snapshot " << snapshot << " can be reached.";
			cont = false;
		}
		vertex_pos.close();
	}
	else {
		LOG(ERROR) << "Cannot read position output file";
		return NULL;
	}

	if (!cont) return NULL;

}