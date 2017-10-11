#version 330 core

// Input vertex data
layout(location = 0) in vec3 vertex_pos_model;

// The model-view-projection matrix
uniform mat4 mvp;

void main(){
    gl_Position = mvp * vec4(vertex_pos_model, 1);
}
