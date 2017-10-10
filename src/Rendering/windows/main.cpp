#include"glew/glew.h"
#include"GLFW/glfw3.h"

#include"common.h"

int main() {

	logger::Logger::default_init("Rendering.log");

	// Initialize GLFW
	if (!glfwInit()) {
		LOG(ERROR) << "Failed to initialize GLFW";
		return -1;
	}

	glfwWindowHint(GLFW_SAMPLES, 4); // 4x antialiasing
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3); // We want OpenGL 3.3
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE); // We don't want the old OpenGL

	// Open a window and create its OpenGL context
	GLFWwindow* window;
	window = glfwCreateWindow(1024, 768, "test", NULL, NULL);
	if (window == NULL) {
		LOG(ERROR) << "Failed to open GLFW window.";
		glfwTerminate();
		return -1;
	}

	glfwMakeContextCurrent(window); // Initialize GLEW
	glewExperimental = true; // Needed in core profile
	if (glewInit() != GLEW_OK) {
		LOG(ERROR) << "Failed to initialize GLEW";
		return -1;
	}

	// Ensure we can capture the escape key being pressed below
	glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);

	do {
		// Draw nothing, see you in tutorial 2 !


		// Swap buffers
		glfwSwapBuffers(window);
		glfwPollEvents();

	} // Check if the ESC key was pressed or the window was closed
	while (glfwGetKey(window, GLFW_KEY_ESCAPE) != GLFW_PRESS &&
		glfwWindowShouldClose(window) == 0);

	LOG(INFO) << glGetString(GL_VERSION);

	system("pause");

	return 0;
}
