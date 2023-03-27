#version 330 core
out vec4 FragColor;
struct Material{
    sampler2D diffuse;
    sampler2D specular;
    float shininess;
};
// in vec3 Normal;
in vec3 FragPos;
in vec2 teseTexCoords;
in vec4 FragPosLightSpace;



uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;
uniform vec3 light_pos;
uniform vec3 view_pos;
uniform vec3 object_color;
uniform vec3 light_color;
uniform Material material;
uniform sampler2D shadow_map;
uniform sampler2D normal_map;
uniform samplerCube skybox;
// in vec3 CubeMapCoords;

float ShadowCalc(){
    vec3 projCoords = FragPosLightSpace.xyz / FragPosLightSpace.w;
    projCoords = projCoords * 0.5 + 0.5;
    if(projCoords.x<0.0||projCoords.x>1.0||projCoords.y<0.0||projCoords.y>1.0)
        return 1.0;
    float depth = texture(shadow_map,projCoords.xy).r;
    float currentDepth = projCoords.z;
    float shadow = currentDepth < depth+0.01? 1.0:0.0;
    return shadow;
}
void main()
{
    // vec3 normal = transpose(inverse(mat3(view))) * Normal;
    vec3 normal = texture(normal_map, teseTexCoords).rgb;
    // vec3 normal = transpose(inverse(mat3(view))) * Normal;
    // vec3 normal = Normal;
    normal = normalize(normal);
    // vec3 view_dir = normalize(vec3(view as* vec4(view_pos, 1)) - FragPos);
    vec3 view_dir = normalize(-FragPos);
    vec3 light_dir = light_pos - FragPos;
    float distance = length(light_dir);
    distance = distance * distance;
    light_dir = normalize(light_dir);
    vec3 reflect_dir = reflect(-light_dir, normal);
    vec3 diffuse = max(dot(normal, light_dir), 0.0) * light_color * 0.5 / distance;

    float spec = pow(max(dot(view_dir,reflect_dir),0.0), 20.0) * 0.8;
    vec3 specular = spec * light_color / distance;
    vec3 ambient = vec3(0.6); 
    FragColor = vec4(normal, 1.0); 
    // FragColor = vec4(diffuse + specular + ambient, 1.0);

    reflect_dir = reflect(-view_dir, normal);
    float shadow = ShadowCalc();
    vec3 result = (ambient + diffuse + shadow * specular);

    // vec3 result = (ambient + diffuse) * texture(material.diffuse,TexCoords).rgb+ specular*texture(material.specular,TexCoords).rgb;

    FragColor = vec4(result, 1.0);
    // FragColor = vec4(normal - normalize(Normal) + vec3(1.0) / 2, 1.0);
}
