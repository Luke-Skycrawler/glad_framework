#version 330 core
out vec4 FragColor;
struct Material{
    sampler2D diffuse;
    sampler2D specular;
    float shininess;
};
in vec3 Normal;
in vec3 FragPos;
in vec2 TexCoords;
in vec4 FragPosLightSpace;



uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;
uniform vec3 light_pos;
uniform vec3 view_pos;
uniform vec3 object_color;
uniform vec3 light_color;
uniform Material material;
uniform sampler2D shadowMap;
// float ShadowCalc(){
//     vec3 projCoords = FragPosLightSpace.xyz / FragPosLightSpace.w;
//     projCoords = projCoords * 0.5 + 0.5;
//     if(projCoords.x<0.0||projCoords.x>1.0||projCoords.y<0.0||projCoords.y>1.0)
//         return 1.0;
//     float depth = texture(shadowMap,projCoords.xy).r;
//     float currentDepth = projCoords.z;
//     float shadow = currentDepth < depth+0.01? 1.0:0.0;
//     return shadow;
// }
void main()
{
     mat3 trans = mat3(view);
    trans = transpose(inverse(trans));
    vec3 normal =  normalize(trans * Normal);
    vec3 view_dir = normalize(view_pos - FragPos);
    vec3 light_dir = normalize(light_pos - FragPos);
    vec3 reflect_dir = reflect(-light_dir, normal);
    vec3 diffuse = max(dot(normal, light_dir), 0.0) * light_color;

    float spec = pow(max(dot(view_dir,reflect_dir),0.0),50.0);
    vec3 specular = spec * light_color;
    vec3 ambient = vec3(0.1); 
    FragColor = vec4(normal, 1.0); 
    // FragColor = vec4(diffuse + specular + ambient, 1.0);

    // float spec = pow(max(dot(viewDir,reflectDir),0.0),material.shininess);
    // vec3 specular = specularStrength * spec * lightColor;
    vec3 result = (ambient + diffuse) * texture(material.diffuse,TexCoords).rgb+specular;

    // // float shadow = ShadowCalc();
    // vec3 result = (ambient + diffuse) * texture(material.diffuse,TexCoords).rgb+ specular*texture(material.specular,TexCoords).rgb;

    FragColor = vec4(result, 1.0);
}
