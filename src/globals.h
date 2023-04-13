#pragma once
#include <glm/glm.hpp>
#include <nlohmann/json.hpp>
#include "dynamics.h"

using json = nlohmann::json;
struct vRender {
    Vector3f pos, normals;
};
struct shayMesh
{
    vector<Vector3f> vertices;
    vector<vRender> vs;
    vector<unsigned> indices, _indices;
    int nv, nf, nvn;
    shayMesh(const string &filename);
    shayMesh(const vector<vec3> &xcs, const vector<unsigned> &indices);
    unsigned vao, vbo, ebo;
    void setupMesh();
    void draw();
    void update_positions();
    void compute_normals(bool clockwise = true);
};

struct Globals {
    MassSpringMesh *mesh;
    shayMesh* rendered_mesh;
    json config;
    Vector2d xm, dxm;
    int selected;
    Matrix4d P_inv, P;
};