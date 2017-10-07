#include"glew/glew.h"
#include"GLFW/glfw3.h"

#include"common.h"

int main() {

	// Initialize GLFW
	if (!glfwInit()) {
		LOG(ERROR) << "Failed to initialize GLFW";
		return -1;
	}

	return 0;
}
