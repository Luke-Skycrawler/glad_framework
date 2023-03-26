#version 410 core

layout(vertices = 4) out;


uniform sampler2D disp_map;
uniform float lod;

// in vec3 FragPos[];
// in vec3 Normal[];
// in vec2 TexCoords[];
// in vec4 FragPosLightSpace[];
// in vec3 CubeMapCoords[];


// out vec3 oFragPos[];
// out vec3 oNormal[];
// out vec2 oTexCoords[];
// out vec4 oFragPosLightSpace[];
// out vec3 oCubeMapCoords[];



void main() {
    gl_TessLevelOuter[0] = lod;
    gl_TessLevelOuter[1] = lod;
    gl_TessLevelOuter[2] = lod;
    gl_TessLevelOuter[3] = lod;
    gl_TessLevelInner[0] = lod;
    gl_TessLevelInner[1] = lod;
    gl_out[gl_InvocationID].gl_Position = gl_in[gl_InvocationID].gl_Position;
}