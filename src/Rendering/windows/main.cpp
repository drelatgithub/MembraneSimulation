#include"glew/glew.h"
#include"GLFW/glfw3.h"
#include"glm/glm.hpp"
#define GLM_ENABLE_EXPERIMENTAL
#include"glm/gtx/transform.hpp"

#include"common.h"
#include"data_loader.h"
#include"shader_loader.h"

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
	GLuint windowWidth = 1024;
	GLuint windowHeight = 768;
	window = glfwCreateWindow(windowWidth, windowHeight, "Plot", NULL, NULL);
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

	// Vertex array object
	GLuint VertexArrayID;
	glGenVertexArrays(1, &VertexArrayID);
	glBindVertexArray(VertexArrayID);

	// Generate vertex array
	gl_vertex_holder glvh;
	glvh.gen_mesh_vertex_buffer_data("p_out.SimOut", "triangles.txt", 0);

	// Drawing
	// This will identify our vertex buffer
	GLuint vertexbuffer;
	// Generate 1 buffer, put the resulting identifier in vertexbuffer
	glGenBuffers(1, &vertexbuffer);
	// The following commands will talk about our 'vertexbuffer' buffer
	glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
	// Give our vertices to OpenGL.
	glBufferData(GL_ARRAY_BUFFER, glvh.num_tris * 9 * sizeof(GLfloat), glvh.vertex_buffer_data.data(), GL_STATIC_DRAW);


	// Ensure we can capture the escape key being pressed below
	glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);


	// Load shaders
	GLuint program_id = LoadShaders("simple_vertex_shader.glsl", "simple_fragment_shader.glsl");

	// Transformations
	glm::mat4 model = glm::mat4(1); // Identity matrix
	glm::mat4 view = glm::lookAt(glm::vec3(4e-6, 3e-6, 4e-6), glm::vec3(0, 0, 0), glm::vec3(0, 1, 0));
	glm::mat4 projection = glm::perspective(
		glm::quarter_pi<double_t>(),
		(double_t)windowWidth / (double_t)windowHeight,
		0.1e-6, 100.0e-6
	);
	glm::mat4 mvp = projection * view * model;

	GLuint matrixId = glGetUniformLocation(program_id, "mvp");

	// Main Loop
	do {

		// Clear screen
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		// Apply transformation
		glUniformMatrix4fv(matrixId, 1, GL_FALSE, &mvp[0][0]);
		// Use shader
		glUseProgram(program_id);

		// 1st attribute buffer : vertices
		glEnableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
		glVertexAttribPointer(
			0,                  // attribute 0. No particular reason for 0, but must match the layout in the shader.
			3,                 // size
			GL_FLOAT,           // type
			GL_FALSE,           // normalized?
			0,                  // stride
			(void*)0            // array buffer offset
		);
		// Draw the triangle !
		glDrawArrays(GL_TRIANGLES, 0, glvh.num_tris*3); // Starting from vertex 0; 6 vertices total -> 2 triangle
		glDisableVertexAttribArray(0);


		// Swap buffers
		glfwSwapBuffers(window);
		glfwPollEvents();

	} // Check if the ESC key was pressed or the window was closed
	while (glfwGetKey(window, GLFW_KEY_ESCAPE) != GLFW_PRESS &&
		glfwWindowShouldClose(window) == 0);

	system("pause");

	return 0;
}
