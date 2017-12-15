#include"glfw/glfw3.h"
#include"glm/glm.hpp"
#include"glm/gtc/constants.hpp"

#include"controls.h"

glm::vec3 position = glm::vec3(0, 0, 5);

GLfloat horizontalAngle = glm::pi<GLfloat>();
GLfloat verticalAngle = 0;
GLfloat initialFoV = glm::quarter_pi<GLfloat>();

GLfloat speed = 3;
GLfloat mouseSpeed = 0.005;