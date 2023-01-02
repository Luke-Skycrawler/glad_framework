#define LIGHT_CPP
#include "light.h"


Light::Light(glm::vec3 LightPositions[], int num):lightCubeShader("shaders/light/light.vert", "shaders/light/light.frag"),lightingShader("shaders/light/multi_light.vert", "shaders/light/multi_light.frag")
{
    pointLightPositions = LightPositions;
    LightNum = num;
    // directional light
    lightingShader.setVec3("dirLight.direction", -0.2f, -1.0f, -0.3f);
    lightingShader.setVec3("dirLight.ambient", 0.05f, 0.05f, 0.05f);
    lightingShader.setVec3("dirLight.diffuse", 0.4f, 0.4f, 0.4f);
    lightingShader.setVec3("dirLight.specular", 0.5f, 0.5f, 0.5f);

    for (int i = 0; i < LightNum; i++)
    {
        lightingShader.setVec3("pointLights[" + to_string(i) + "].position", pointLightPositions[i]);
        lightingShader.setVec3("pointLights[" + to_string(i) + "].ambient", 0.05f, 0.05f, 0.05f);
        lightingShader.setVec3("pointLights[" + to_string(i) + "].diffuse", 0.8f, 0.8f, 0.8f);
        lightingShader.setVec3("pointLights[" + to_string(i) + "].specular", 1.0f, 1.0f, 1.0f);
        lightingShader.setFloat("pointLights[" + to_string(i) + "].constant", 1.0f);
        lightingShader.setFloat("pointLights[" + to_string(i) + "].linear", 0.09);
        lightingShader.setFloat("pointLights[" + to_string(i) + "].quadratic", 0.032);
    }
}

void Light::SetShaderValue(Camera &camera)
{

    // spotLight
    lightingShader.setVec3("spotLight.position", camera.Position);
    lightingShader.setVec3("spotLight.direction", camera.Front);
    lightingShader.setVec3("spotLight.ambient", 0.0f, 0.0f, 0.0f);
    lightingShader.setVec3("spotLight.diffuse", 1.0f, 1.0f, 1.0f);
    lightingShader.setVec3("spotLight.specular", 1.0f, 1.0f, 1.0f);
    lightingShader.setFloat("spotLight.constant", 1.0f);
    lightingShader.setFloat("spotLight.linear", 0.09);
    lightingShader.setFloat("spotLight.quadratic", 0.032);
    lightingShader.setFloat("spotLight.cutOff", glm::cos(glm::radians(12.5f)));
    lightingShader.setFloat("spotLight.outerCutOff", glm::cos(glm::radians(15.0f)));
}
void Light::Draw(Camera &camera){
    // we now draw as many light bulbs as we have point lights.
    // glBindVertexArray(lightCubeVAO);
    // also draw the lamp object(s)
    lightCubeShader.use();
    glm::mat4 projection = glm::perspective(glm::radians(camera.Zoom), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 100.0f);
    glm::mat4 view = camera.GetViewMatrix();
    glm::mat4 model;
    lightCubeShader.setMat4("projection", projection);
    lightCubeShader.setMat4("view", view);
    for (int i = 0; i < LightNum; i++)
    {
        model = glm::translate(glm::mat4(1.0f), pointLightPositions[i]);
        model = glm::scale(model, glm::vec3(0.2f)); // Make it a smaller cube
        lightCubeShader.setMat4("model", model);
        renderCube(1);
    }
}