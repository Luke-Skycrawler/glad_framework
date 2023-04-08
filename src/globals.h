#pragma once
#include <glm/glm.hpp>
#include <nlohmann/json.hpp>
#include "dynamics.h"

using json = nlohmann::json;
struct shayMesh
{
    vector<Vector3f> vertices;
    vector<unsigned> indices;
    int nv, nf;
    shayMesh(const string &filename);
    shayMesh(const vector<vec3> &xcs, const vector<unsigned> &indices);
    unsigned vao, vbo, ebo;
    void setupMesh();
    void draw();
    void update_positions();
};

struct Globals {
    MassSpringMesh *mesh;
    shayMesh* rendered_mesh;
    json config;
};