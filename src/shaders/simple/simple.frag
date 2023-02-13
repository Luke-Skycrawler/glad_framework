#version 330 core

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;
uniform vec3 light_pos;
uniform vec3 view_pos;
uniform vec3 object_color;
uniform vec3 light_color;

out vec4 FragColor;
in vec3 FragPos;
in vec3 Normal;
void main()
{
    mat3 trans = mat3(view);
    trans = transpose(inverse(trans));
    vec3 normal =  trans * Normal;
    vec3 view_dir = normalize(view_pos - FragPos);
    vec3 light_dir = normalize(light_pos - FragPos);
    vec3 reflect_dir = reflect(-light_dir, normal);
    vec3 diffuse = max(dot(normal, light_dir), 0.0) * light_color;

    float spec = pow(max(dot(view_dir,reflect_dir),0.0),2.0);
    vec3 specular = spec * light_color;
    vec3 ambient = vec3 (0.7, 0.0, 0.0); 
    FragColor = vec4(normal, 1.0); 
    FragColor = vec4(diffuse + specular + ambient, 1.0);
}