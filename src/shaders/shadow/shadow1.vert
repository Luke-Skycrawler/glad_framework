#version 330 core
layout (location = 0) in vec3 aPos;
layout (location = 1) in vec3 aNormal;
layout (location = 2) in vec2 aTexCoords;

out vec3 FragPos;
out vec3 Normal;
out vec2 TexCoords;
out vec4 FragPosLightSpace;
out vec3 CubeMapCoords;

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;
uniform mat4 lightView;
void main()
{
    FragPos = vec3(view * model * vec4(aPos, 1.0));
    // FragPos.y*= -1.0;
    Normal = vec3(model * vec4(aNormal,0.0));  
    // Normal.y *= -1.0;
    TexCoords = aTexCoords;
    FragPosLightSpace = lightView * vec4(FragPos,1.0);
    CubeMapCoords = mat3(model) * aPos;
    vec4 a = model * vec4(aPos, 1.0);
    a.y *= -1.0;
    gl_Position = projection * view * a;
}