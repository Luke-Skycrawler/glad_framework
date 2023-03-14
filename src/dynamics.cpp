#include "dynamics.h"
#include <iostream>
using namespace std;
using namespace glm;
static const float k = 3, m = 1, alpha_damp = 1.0f, bound = 5;
void Particle::force(float xc, float yc)
{
    vec2 pc{xc, yc};
    f = k * (pc - pos);
    //cout << pc[0] << " " << pc[1] << "\n";
}

mat3 skew(vec3 r)
{
    glm::mat3 ret(0.0);
    ret[1][0] = -r[2];
    ret[0][1] = +r[2];
    ret[2][0] = +r[1];
    ret[0][2] = -r[1];
    ret[2][1] = -r[0];
    ret[1][2] = +r[0];
    return ret;
}
void Particle::step_implicit(float dt)
{
    mat3 m = mat3(1.0) - (skew(vec3(0.0, 0.0, 1.0)) * dt * dt);

    vec3 b(pos + dt * v, 0.0f);
    auto J = determinant(m);
    vec3 a1(m[0]), a2(m[1]), a3(m[2]);
    auto x1 = determinant(mat3(b, a2, a3)) / J;
    auto x2 = determinant(mat3(a1, b, a3)) / J;
    auto x3 = determinant(mat3(a1, a2, b)) / J;
    auto opos = pos;
    pos = vec2(x1, x2);
    v = -(opos - pos) / dt;
    bounce();
    cout << dot(pos, pos) << " " << dt << "\n";
}
void Particle::concentric_explicit()
{
    auto tang = cross(vec3(pos, 0.0), vec3(0.0, 0.0, 1.0));
    f = tang * k;
}
void Particle::step(float dt, bool centre)
{
    if (centre)
        concentric_explicit();
    v += f * dt / m;
    v *= alpha_damp;
    pos += v * dt;
    bounce();
}

void Particle::bounce(){
    float& x = pos[0], &y = pos[1];
    const auto reflect = [](float &x, float &v){
        if (x < -bound) {
            x = -2 * bound - x;
            v = -v;
        }
        if (x > bound) {
            x = 2 * bound - x;
            v = -v;
        }
        assert(x < bound && x > -bound);
    };
    reflect(x, v[0]);
    reflect(y, v[1]);
}