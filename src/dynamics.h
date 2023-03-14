#pragma once
#include <glm/glm.hpp>
struct Particle {
    glm::vec2 pos, f, v;
    void force(float x, float y);
    void step(float dt, bool centre = false);
    void concentric_explicit();
    void step_implicit(float dt);
    void bounce();
};