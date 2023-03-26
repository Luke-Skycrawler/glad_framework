#version 330 core

layout(triangles) in;
layout(line_strip, max_vertices = 3) out;

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

void main() {
    for (int i = 0; i < 3; i ++) {
    // gl_Position = projection * (view * model * gl_in[i].gl_Position + vec4(0.0, 0.0, -0.05, 0.0));

    // if no tessllation shader, use this line
    // gl_Position = projection * view * model * gl_in[i].gl_Position + vec4(0.0, 0.0, -0.05, 0.0);

    // else remove projection view model
    gl_Position = gl_in[i].gl_Position + vec4(0.0, 0.0, -0.05, 0.0);
    EmitVertex();
    }
    
}