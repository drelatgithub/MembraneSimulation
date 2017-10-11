#include"shader_loader.h"

#include"common.h"

GLuint LoadShaders(const char* vertex_file_path, const char* fragment_file_path){
	
	// Create shaders
	GLuint vertex_shader_id = glCreateShader(GL_VERTEX_SHADER);
	GLuint fragment_shader_id = glCreateShader(GL_FRAGMENT_SHADER);

	// Read vertex shader from file
	std::string vertex_shader_code;
	std::ifstream vertex_shader_stream(vertex_file_path, std::ios::in);
	if (vertex_shader_stream.is_open()) {
		std::string line = "";
		while (std::getline(vertex_shader_stream, line))
			vertex_shader_code += "\n" + line;
		vertex_shader_stream.close();
	}
	else {
		LOG(ERROR) << "Cannot open vertex shader file.";
		return 0;
	}

	// Read fragment shader from file
	std::string fragment_shader_code;
	std::ifstream fragment_shader_stream(fragment_file_path, std::ios::in);
	if (fragment_shader_stream.is_open()) {
		std::string line = "";
		while (std::getline(fragment_shader_stream, line))
			fragment_shader_code += "\n" + line;
		fragment_shader_stream.close();
	}
	else {
		LOG(WARNING) << "Cannot open fragment shader file.";
	}

	GLint res = GL_FALSE;
	int info_log_len;

	// Compile vertex shader
	LOG(INFO) << "Compiling shader: " << vertex_file_path;
	char const * vertex_src_ptr = vertex_shader_code.c_str();
	glShaderSource(vertex_shader_id, 1, &vertex_src_ptr, NULL);
	glCompileShader(vertex_shader_id);

	// Check vertex shader
	glGetShaderiv(vertex_shader_id, GL_COMPILE_STATUS, &res);
	glGetShaderiv(vertex_shader_id, GL_INFO_LOG_LENGTH, &info_log_len);
	if (info_log_len > 0) {
		std::vector<char> vertex_shader_err_msg(info_log_len + 1);
		glGetShaderInfoLog(vertex_shader_id, info_log_len, NULL, &vertex_shader_err_msg[0]);
		LOG(ERROR) << &vertex_shader_err_msg[0];
	}

	// Compile fragment shader
	LOG(INFO) << "Compiling shader: " << fragment_file_path;
	char const * fragment_src_ptr = fragment_shader_code.c_str();
	glShaderSource(fragment_shader_id, 1, &fragment_src_ptr, NULL);
	glCompileShader(fragment_shader_id);

	// Check fragment shader
	glGetShaderiv(fragment_shader_id, GL_COMPILE_STATUS, &res);
	glGetShaderiv(fragment_shader_id, GL_INFO_LOG_LENGTH, &info_log_len);
	if (info_log_len > 0) {
		std::vector<char> fragment_shader_err_msg(info_log_len + 1);
		glGetShaderInfoLog(fragment_shader_id, info_log_len, NULL, &fragment_shader_err_msg[0]);
		LOG(ERROR) << &fragment_shader_err_msg[0];
	}

	// Link the program
	LOG(INFO) << "Linking program";
	GLuint program_id = glCreateProgram();
	glAttachShader(program_id, vertex_shader_id);
	glAttachShader(program_id, fragment_shader_id);
	glLinkProgram(program_id);

	// Check the program
	glGetProgramiv(program_id, GL_LINK_STATUS, &res);
	glGetProgramiv(program_id, GL_INFO_LOG_LENGTH, &info_log_len);
	if (info_log_len > 0) {
		std::vector<char> program_err_msg(info_log_len + 1);
		glGetProgramInfoLog(program_id, info_log_len, NULL, &program_err_msg[0]);
		LOG(ERROR) << &program_err_msg[0];
	}

	glDetachShader(program_id, vertex_shader_id);
	glDetachShader(program_id, fragment_shader_id);

	glDeleteShader(vertex_shader_id);
	glDeleteShader(fragment_shader_id);

	return program_id;

}