#pragma once
#include <Eigen/Eigen>
#include <vector>
#include <array>

using namespace std;
using namespace Eigen;
using vec3 = Vector3d;
using mat3 = Matrix3d;

static const double dt = 1e-3,
                    bound = 1.0,
                    eps = 1e-1;
inline mat3 skew(const vec3 & r){
    mat3 ret;
    ret.setZero(3,3);
    ret(0, 1) = -r(2);
    ret(1, 0) = +r(2);
    ret(0, 2) = +r(1);
    ret(2, 0) = -r(1);
    ret(1, 2) = -r(0);
    ret(2, 1) = +r(0);
    return ret;
}


struct StateVector {
    vec3 x, p, L;
    mat3 R;
    inline StateVector operator +(const StateVector &S){
        return {x + S.x, p + S.p, L + S.L, R + S.R};
    }
    inline StateVector operator *(const double k){
        return {k * x, k * p, k * L, k * R};
    }
};

struct RigidBody{
    StateVector S[2];
    double M_inv;
    mat3 J0;
    int n_vertices; 
    vector<int> col_set[6];
    vec3 *vertices;
    
    RigidBody(int n_vertices = 0, vec3 * vertices = nullptr): n_vertices(n_vertices), vertices(vertices), M_inv(1.0){
        J0 = mat3::Identity(3,3) / 12.0;
        S[0] = {vec3(0.0, 0.0, 0.0), vec3(1.0, 0.0, 0.0), vec3(0.1, 0.0, 0.1), mat3::Identity(3, 3)};
    }

    mat3 compute_K(const vec3 &r, int t);
    StateVector dSdt(const vec3 &f, int t);
    void step(int ts);
};


static double kd = 1e-2, ks = 1e1;

struct Edge {
    double l0;
    int i, j;
};
struct MassSpringMesh {
    vector<vec3> mass_v;
    vector<vec3> mass_x;
    vector<Edge> edges;
};
// void compute_b(VectorXd &b, const VectorXd &v_plus);
// // M * v_t + dt * f(x_t + v_t+1 * dt)
// mat3 compute_single_spring_K(Edge &e, const vec3 &vi, const vec3 &vj);
// void compute_force(VectorXd &b, const VectorXd &v_plus);
// // b += dt * f(x + v * dt)
// void init_l0();
void implicit_euler();
void init();
void extract_edges(vector<Edge>& edges, const vector<unsigned> indices);
