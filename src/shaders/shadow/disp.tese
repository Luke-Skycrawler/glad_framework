#version 410 core

layout(quads, equal_spacing, ccw) in;

in vec2 tescTexCoords[];
// in vec3 tescFragPos[];
out vec2 teseTexCoords;
out vec3 FragPos;
out vec3 frag_pos_world;
out vec4 FragPosLightSpace;
uniform sampler2D disp_map;

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;
uniform float popup;
uniform mat4 lightView;

vec4 interpolate(vec4 v0, vec4 v1, vec4 v2, vec4 v3) {
    vec4 a = mix(v0, v1, gl_TessCoord.x);
    vec4 b = mix(v3, v2, gl_TessCoord.x);
    
    return mix(a, b, gl_TessCoord.y);
}

vec2 interpolate(vec2 v0, vec2 v1, vec2 v2, vec2 v3) {
    vec2 a = mix(v0, v1, gl_TessCoord.x);
    vec2 b = mix(v3, v2, gl_TessCoord.x);
    
    return mix(a, b, gl_TessCoord.y);
}
void main() {

    vec4 pos = interpolate( 
        gl_in[0].gl_Position, 
        gl_in[1].gl_Position, 
        gl_in[2].gl_Position, 
        gl_in[3].gl_Position 
    );
    gl_Position = projection * view * model * ( pos + 
        vec4(
            0.0, 
            0.0, 
            texture(disp_map, vec2(gl_TessCoord.x, gl_TessCoord.y)).r * popup, 
            0.0
        )
    );

    teseTexCoords = interpolate(
        tescTexCoords[0],
        tescTexCoords[1],
        tescTexCoords[2],
        tescTexCoords[3]
    );
    FragPos = vec3(view * model * pos);
    FragPosLightSpace = lightView * model * pos;
    frag_pos_world = vec3(model * pos);
}