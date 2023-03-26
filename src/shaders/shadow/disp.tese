#version 410 core

layout(quads, equal_spacing, ccw) in;

in vec2 TexCoords[];
uniform sampler2D disp_map;

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

vec4 interpolate(vec4 v0, vec4 v1, vec4 v2, vec4 v3) {
    vec4 a = mix(v0, v1, gl_TessCoord.x);
    vec4 b = mix(v3, v2, gl_TessCoord.x);
    
    return mix(a, b, gl_TessCoord.y);
}
void main() {
    gl_Position = projection * view * model * 
    interpolate(
        gl_in[0].gl_Position, 
        gl_in[1].gl_Position, 
        gl_in[2].gl_Position, 
        gl_in[3].gl_Position 
    );
    
}