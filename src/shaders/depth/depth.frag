#version 330 core
out vec4 FragColor;
in vec3 frag_pos_world;

uniform vec3 view_pos;
void main()
{
    // vec3 viewDir= normalize(viewPos-FragPos);
    vec3 depth = vec3(length(view_pos-frag_pos_world));
    FragColor = vec4(depth, 1.0);
}
