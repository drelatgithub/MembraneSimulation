#include"common.h"
#include"data_loader.h"

const GLfloat* gen_mesh_vertex_buffer_data(const char* pos_out_path, const char* tri_path, size_t snapshot) {
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
			std::vector<double> coords;

			while (iss >> each_num) {
				coords.push_back(each_num);
			}

			LOG(INFO) << coords.size();
			return NULL;
		}
		else {
			LOG(ERROR) << "File ended before snapshot " << snapshot << " can be reached.";
			return NULL;
		}
	}
	else {
		LOG(ERROR) << "Cannot read position output file";
		return NULL;
	}

}