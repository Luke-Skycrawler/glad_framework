#include "dynamics.h"
#include <spdlog/spdlog.h>
#include <assert.h>
#include <iostream>
mat3 RigidBody::compute_K(const vec3 &r, int t)
{
    auto &R{S[t].R};
    mat3 Rr{skew(R * r)};
    mat3 J_inv{R * J0.inverse() * R.transpose()};
    return M_inv * mat3::Identity(3, 3) - Rr * J_inv * Rr;
}

StateVector RigidBody::dSdt(const vec3 &f, int t)
{
    auto &R{S[t].R};
    auto &pdot = f;
    auto &xdot = M_inv * S[t].p;
    auto &Rdot = skew(R * J0.inverse() * R.transpose() * S[t].L) * R;
    vec3 Ldot(0.0, 0.0, 0.0);
    return {xdot, pdot, Ldot, Rdot};
}

void RigidBody::step(int ts)
{
    auto Sdot = dSdt(vec3(0.0, 0.0, 0.0), 0);

    S[1] = S[0] + Sdot * dt;
    for (int i = 0; i < 6; i++)
        col_set[i].resize(0);
    auto R = S[1].R;

    for (int i = 0; i < n_vertices; i++)
    {
        vec3 &r{vertices[i]}, xi{R * r + S[1].x};
        for (int k = 0; k < 3; k++)
        {
            if (xi(k) > bound)
                col_set[k * 2].push_back(i);
            else if (xi(k) < -bound)
                col_set[k * 2 + 1].push_back(i);
        }
        for (int i = 0; i < 6; i++)
            if (col_set[i].size())
                spdlog::info("collsiion set {} size = {}", i, col_set[i].size());

        vec3 dp(0.0, 0.0, 0.0), dL(0.0, 0.0, 0.0), dx(0.0, 0.0, 0.0);
        int tot = 0;
        for (int i = 0; i < 6; i++)
        {
            tot += col_set[i].size();
            vec3 r(0.0, 0.0, 0.0);
            for (int I : col_set[i])
            {
                vec3 &rr{vertices[I]};
                r += rr;
            }
            if (col_set[i].size())
            {
                r /= col_set[i].size();
                // vec3 r{vertices[I]};
                vec3 Rr{R * r},
                    xi{R * r + S[1].x};
                auto K = compute_K(r, 1);

                vec3 omega = R * J0.inverse() * R.transpose() * S[1].L;
                vec3 vi = M_inv * S[1].p + skew(omega) * Rr;

                vec3 dv(0.0, 0.0, 0.0);

                // for (int k = 0; k < 3; k++) {
                int k = i / 2;
                if (abs(xi(k)) > bound)
                {
                    dv(k) += -vi(k) * (1.0 + eps);
                    if (xi(k) > 0.0)
                        dx(k) += bound - xi(k);
                    else
                        dx(k) += -bound - xi(k);
                }
                // }
                vec3 deltap = K.inverse() * dv;
                dp += deltap;
                dL += skew(Rr) * deltap;

                vec3 omega1 = R * J0.inverse() * R.transpose() * (S[1].L + dL);
                vec3 vi1 = M_inv * (S[1].p + dp) + skew(omega1) * Rr;

                if (dv(k) != 0.0)
                {
                    spdlog::info("dim = {}, vi0 = {:.6f}, vi1 = {:.6f}", k, vi(k), vi1(k));
                }
            }
        }
        if (tot)
        {
            spdlog::info("ts = {}", ts);
        }

        S[1].p += dp;
        S[1].L += dL;
        S[0] = S[1];
        S[0].x += dx;
    }
}

#include "globals.h"
extern Globals globals;

static const double M = 1.0;

void compute_force(VectorXd &b, const VectorXd &v_plus)
{
    // b += dt * f(x + dt * v_plus)
    auto &xs{globals.mesh->mass_x};
    static const vec3 gravity{0.0, -9.8, 0.0};
    for (auto &e : globals.mesh->edges)
    {
        int i = e.i, j = e.j;
        vec3 xji{xs[j] - xs[i]};
        xji += (v_plus.segment<3>(3 * j) - v_plus.segment<3>(3 * i)) * dt;

        assert(xji.squaredNorm() > 0.0);
        vec3 fi = ks * (xji - e.l0 * (xji).normalized());
        if (fi.squaredNorm() > 0.01)
        {
            // cout << fi.transpose() << "\n";
        }
        b.segment<3>(3 * i) += fi * dt;
        b.segment<3>(3 * j) -= fi * dt;
    }
    int n_mass = globals.mesh->mass_x.size();

    for (int i = 0; i < n_mass; i++)
    {
        b.segment<3>(i * 3) += M * gravity * dt;

        // vec3 vi{v_plus.segment<3>(i * 3)};
        // vec3 xi{xs[i] + vi * dt};
        // vec3 dv{0.0, 0.0, 0.0};
        // vec3 dx{0.0, 0.0, 0.0};
        // for (int k = 0; k < 3; k++)
        //     if (abs(xi(k)) > bound && vi(k) * xi(k) > 0.0)
        //     {
        //         dv(k) += -vi(k) * (1.0 + eps);

        //         // if (xi(k) > 0.0)
        //         //     dx(k) += bound - xi(k);
        //         // else
        //         //     dx(k) += -bound - xi(k);
        //     }
        // b.segment<3>(i * 3) += M * dv;
    }
}
void compute_b(VectorXd &b, const VectorXd &v_plus)
/****************************************
b = M * v_t + dt * f(x_t + v_t+1 * dt)
****************************************/
{
    auto &vs{globals.mesh->mass_v}, &xs{globals.mesh->mass_x};

    int n_mass = vs.size();
    b.setZero(3 * n_mass);
    for (int i = 0; i < n_mass; i++)
    {
        b.segment<3>(i * 3) = M * vs[i];
    }
    compute_force(b, v_plus);

    // cout << b.transpose() << "\n";
}

void gen_non_zero_entries(SparseMatrix<double> &sparse_matrix)
{
}

mat3 compute_single_spring_K(Edge &e, const vec3 &vi, const vec3 &vj, bool use_v_t = true)
{
    // partial f partial xi
    auto &xs{globals.mesh->mass_x};
    auto &vs{globals.mesh->mass_v};
    vec3 xji{xs[e.j] - xs[e.i]};
    xji += dt * (vj - vi);

    vec3 vji{vs[e.j] - vs[e.i]};
    if (!use_v_t)
        vji = (vj - vi);

    mat3 K = (xji) * (xji).transpose();
    double l2 = xji.squaredNorm();
    double l = sqrt(l2);
    mat3 term1 = mat3::Identity(3, 3) - K / l2;
    mat3 K_spring = ks * (-mat3::Identity(3, 3) + (e.l0 / l) * term1);
    // mat3 K_damp = -(kd / l) * vji.transpose() * term1;
    return K_spring;
    // return K_spring + K_damp;
}

void compute_A(SparseMatrix<double> &sparse_matrix, const VectorXd &v_n, bool use_v_t = true)
{
    sparse_matrix.setZero();
    auto &xs{globals.mesh->mass_x};
    int n_mass = globals.mesh->mass_x.size();
    auto &vs{globals.mesh->mass_v};
    for (int i = 0; i < n_mass * 3; i++)
    {
        sparse_matrix.coeffRef(i, i) = M;
    }
    static const double h2_neg = -dt * dt;
    for (auto &e : globals.mesh->edges)
    {
        auto I = e.i * 3, J = e.j * 3;
        vec3 vi, vj;
        vi = v_n.segment<3>(I);
        vj = v_n.segment<3>(J);
        mat3 K = compute_single_spring_K(e, vi, vj, use_v_t);
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            {
                sparse_matrix.coeffRef(i + I, j + I) += K(i, j) * h2_neg;
                sparse_matrix.coeffRef(i + I, j + J) += -K(i, j) * h2_neg;
                sparse_matrix.coeffRef(i + J, j + I) += -K(i, j) * h2_neg;
                sparse_matrix.coeffRef(i + J, j + J) += K(i, j) * h2_neg;
            }
    }
    // for (int i = 0; i < n_mass; i++)
    // {
    //     vec3 vi{v_n.segment<3>(i * 3)};
    //     vec3 xi{xs[i] + vi * dt};
    //     vec3 dv{0.0, 0.0, 0.0};
    //     vec3 dx{0.0, 0.0, 0.0};
    //     for (int k = 0; k < 3; k++)
    //         if (abs(xi(k)) > bound && vi(k) * xi(k) > 0.0)
    //         {
    //             dv(k) += -vi(k) * (1.0 + eps);
    //             sparse_matrix.coeffRef(k + i * 3, k + i * 3) += - M * (1.0 + eps) / dt;
    //             // if (xi(k) > 0.0)
    //             //     dx(k) += bound - xi(k);
    //             // else
    //             //     dx(k) += -bound - xi(k);
    //         }
    // }
}

void init_l0()
{
    auto &xs{globals.mesh->mass_x};
    for (auto &e : globals.mesh->edges)
    {
        int i = e.i, j = e.j;
        vec3 xji{xs[j] - xs[i]};
        e.l0 = xji.norm();
    }
}

void init()
{
    init_l0();
    // globals.mesh->mass_x[0] += vec3{0.0, 0.0, 0.3};
}

static const double tol = 1e-4;
static const int max_iters = 100;

VectorXd cat(const vector<vec3> &v)
{
    VectorXd ret;
    int n_mass = v.size();
    ret.setZero(n_mass * 3);
    for (int i = 0; i < n_mass; i++)
    {
        ret.segment<3>(i * 3) = v[i];
    }
    return ret;
}

void implicit_euler()
{

    int n_mass = globals.mesh->mass_x.size();
    auto &vs{globals.mesh->mass_v};
    auto &xs{globals.mesh->mass_x};

    int n_unknowns = n_mass * 3;
    VectorXd b, v_plus;
    SparseMatrix<double> sparse_matrix(n_unknowns, n_unknowns);
    b.setZero(n_unknowns);
    v_plus.setZero(n_unknowns);
    int iter = 0;
    do
    {
        // Newton iteration

        compute_b(b, v_plus);
        // 1st iter v_plus = 0, b = M * v_t + dt * f(x_t)

        compute_A(sparse_matrix, v_plus, iter == 0);
        SimplicialLDLT<SparseMatrix<double>> ldlt_solver;
        ldlt_solver.compute(sparse_matrix);
        v_plus = ldlt_solver.solve(b);
        double residue;
        {

            VectorXd f = M * (cat(vs) - v_plus);
            compute_force(f, v_plus);
            // f += dt * f(x + dt * v_plus)
            residue = f.norm();
        }
        bool term_cond = residue < tol;
        if (term_cond)
            break;
        if (++iter >= max_iters)
            break;
    } while (true);
    cout << "iter = " << iter << ", norm v = " << v_plus.norm() << "\n";
    for (int i = 0; i < n_mass; i++)
    {
        auto vi = v_plus.segment<3>(i * 3);
        vs[i] = vi;
        xs[i] += vi * dt;
        vec3 xi{xs[i]};
        vec3 dv{0.0, 0.0, 0.0};
        vec3 dx{0.0, 0.0, 0.0};
        for (int k = 0; k < 3; k++)
            if (abs(xi(k)) > bound && vi(k) * xi(k) > 0.0)
            {
                dv(k) += -vi(k) * (1.0 + eps);

                if (xi(k) > 0.0)
                    dx(k) += bound - xi(k);
                else
                    dx(k) += -bound - xi(k);
            }
        // b.segment<3>(i * 3) += M * dv;
        vs[i] += dv;
        xs[i] += dx;
    }
}

#include <set>
void extract_edges(vector<Edge> &edges, const vector<unsigned> indices)
{
    set<array<unsigned, 2>> e;
    edges.resize(0);
    static const auto insert = [&](unsigned a, unsigned b)
    {
        e.insert({min(a, b), max(a, b)});
    };
    int n_faces = indices.size() / 3;
    for (int i = 0; i < n_faces; i++)
    {
        auto t0 = indices[3 * i], t1 = indices[3 * i + 1], t2 = indices[3 * i + 2];
        insert(t0, t1);
        insert(t1, t2);
        insert(t0, t2);
    }
    edges.reserve(e.size());
    for (auto &ei : e)
    {
        edges.push_back({0.0, int(ei[0]), int(ei[1])});
    }
}
