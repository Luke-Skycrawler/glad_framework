#include "dynamics.h"
#include <spdlog/spdlog.h>
#include <assert.h>
#include <iostream>
#include <set>
#include "globals.h"

extern Globals globals;
static const double M = 1.0;

void compute_force(VectorXd &b, const VectorXd &v_plus)
{
    // b += dt * f(x + dt * v_plus)
    auto &xs{globals.mesh->mass_x};
    static const vec3 gravity{0.0, -9.8, 0.0};
    auto &is_static{globals.mesh->is_static};
    int n_mass = globals.mesh->mass_x.size();

    for (auto &e : globals.mesh->edges)
    {
        int i = e.i, j = e.j;
        vec3 xji{xs[j] - xs[i]};
        xji += (v_plus.segment<3>(3 * j) - v_plus.segment<3>(3 * i)) * dt;

        assert(xji.squaredNorm() > 0.0);
        vec3 fi = globals.config["ks"] * (xji - e.l0 * (xji).normalized());
        if (fi.squaredNorm() > 0.01)
        {
            // cout << fi.transpose() << "\n";
        }
        if (!is_static[i])
            b.segment<3>(3 * i) += fi * dt;
        if (!is_static[j])
            b.segment<3>(3 * j) -= fi * dt;
    }
    if (globals.config["gravity"])
    for (int i = 0; i < n_mass; i++)
        if (!is_static[i])
        {
            b.segment<3>(i * 3) += M * gravity * dt;

            // vec3 vi{v_plus.segment<3>(i * 3)};
            // vec3 xi{xs[i] + vi * dt};
            // vec3 dv{0.0, 0.0, 0.0};
            // vec3 dx{0.0, 0.0, 0.0};
            // for (int k = 0; k < 3; k++)
            //     if (abs(xi(k)) > globals.config["bound"] && vi(k) * xi(k) > 0.0)
            //     {
            //         dv(k) += -vi(k) * (1.0 + eps);

            //         // if (xi(k) > 0.0)
            //         //     dx(k) += globals.config["bound"] - xi(k);
            //         // else
            //         //     dx(k) += -globals.config["bound"] - xi(k);
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
    auto &is_static{globals.mesh->is_static};
    int n_mass = vs.size();
    b.setZero(3 * n_mass);
    for (int i = 0; i < n_mass; i++)
        if (!is_static[i])
        {
            b.segment<3>(i * 3) = M * (vs[i] - v_plus.segment<3>(i * 3));
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
    mat3 K_spring = globals.config["ks"] * (-mat3::Identity(3, 3) + (e.l0 / l) * term1);
    // mat3 K_damp = -(kd / l) * vji.transpose() * term1;
    return K_spring;
    // return K_spring + K_damp;
}

inline int starting_offset(int i, int j, const std::map<std::array<int, 2>, int> &lut, int *outers)
{
    auto it = lut.find({i, j});
    int k = it->second;
    return k * 3 + outers[j * 3];
}

inline int stride(int j, int *outers)
{
    return outers[j * 3 + 1] - outers[j * 3];
    // full 3x3 matrix, no overflow issue
}

void put(double *values, int offset, int _stride, const mat3 &block)
{
    for (int j = 0; j < 3; j++)
        for (int i = 0; i < 3; i++)
        {
            values[offset + _stride * j + i] += block(i, j);
        }
}

void compute_A(SparseMatrix<double> &sparse_matrix, const VectorXd &v_n, bool use_v_t = true)
{
    sparse_matrix.setZero();
    auto &xs{globals.mesh->mass_x};
    int n_mass = globals.mesh->mass_x.size();
    auto &vs{globals.mesh->mass_v};
    auto &is_static{globals.mesh->is_static};
#define FANCY

#ifdef FANCY

    // TODO:  move this to initialization
    map<array<int, 2>, int> lut;
    gen_empty_sm(n_mass, globals.mesh->edges, sparse_matrix, lut);

    auto outers = sparse_matrix.outerIndexPtr();
    auto values = sparse_matrix.valuePtr();
    for (int k = 0; k < n_mass; k++)
    {
        int offset = starting_offset(k, k, lut, outers);
        int _stride = stride(k, outers);
        for (int i = 0; i < 3; i++)
        {
            values[offset + _stride * i + i] = is_static[k] ? 1.0 : M;
        }
    }
#endif
#ifdef NO_FANCY
    for (int i = 0; i < n_mass * 3; i++)
    {
        sparse_matrix.coeffRef(i, i) = M;
    }
#endif
    static const double h2_neg = -dt * dt;
    for (auto &e : globals.mesh->edges)
    {
        auto I = e.i * 3, J = e.j * 3;
        vec3 vi, vj;
        vi = v_n.segment<3>(I);
        vj = v_n.segment<3>(J);
        mat3 K = compute_single_spring_K(e, vi, vj, use_v_t);

#ifdef FANCY
        int ii = e.i, jj = e.j;
        auto stride_j = stride(jj, outers), stride_i = stride(ii, outers);
        auto oii = starting_offset(ii, ii, lut, outers), ojj = starting_offset(jj, jj, lut, outers), oij = starting_offset(ii, jj, lut, outers), oji = starting_offset(jj, ii, lut, outers);
        if (!is_static[ii])
        {
            put(values, oii, stride_i, K * h2_neg);
            if (!is_static[jj])
                put(values, oji, stride_i, -K * h2_neg);
        }
        if (!is_static[jj])
        {
            if (!is_static[ii])
                put(values, oij, stride_j, -K * h2_neg);
            put(values, ojj, stride_j, K * h2_neg);
        }
#endif
#ifdef NO_FANCY
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            {
                sparse_matrix.coeffRef(i + I, j + I) += K(i, j) * h2_neg;
                sparse_matrix.coeffRef(i + I, j + J) += -K(i, j) * h2_neg;
                sparse_matrix.coeffRef(i + J, j + I) += -K(i, j) * h2_neg;
                sparse_matrix.coeffRef(i + J, j + J) += K(i, j) * h2_neg;
            }
#endif
    }
    // for (int i = 0; i < n_mass; i++)
    // {
    //     vec3 vi{v_n.segment<3>(i * 3)};
    //     vec3 xi{xs[i] + vi * dt};
    //     vec3 dv{0.0, 0.0, 0.0};
    //     vec3 dx{0.0, 0.0, 0.0};
    //     for (int k = 0; k < 3; k++)
    //         if (abs(xi(k)) > globals.config["bound"] && vi(k) * xi(k) > 0.0)
    //         {
    //             dv(k) += -vi(k) * (1.0 + eps);
    //             sparse_matrix.coeffRef(k + i * 3, k + i * 3) += - M * (1.0 + eps) / dt;
    //             // if (xi(k) > 0.0)
    //             //     dx(k) += globals.config["bound"] - xi(k);
    //             // else
    //             //     dx(k) += -globals.config["bound"] - xi(k);
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
        auto dv = ldlt_solver.solve(b);
        double residue;
        {
            v_plus += dv;
            VectorXd f = M * (cat(vs) - v_plus);
            compute_force(f, v_plus);
            // f += dt * f(x + dt * v_plus)
            residue = f.norm();
        }
        bool term_cond = residue < globals.config["tol"];
        if (term_cond)
            break;
        if (++iter >= globals.config["max_iters"])
            break;
    } while (true);
    cout << "iter = " << iter << ", norm v = " << v_plus.norm() << "\n";
    for (int i = 0; i < n_mass; i++)
    {
        auto vi = v_plus.segment<3>(i * 3);
        vs[i] = vi;
        xs[i] += vi * dt;

#ifdef BOUND_ENABLED
        vec3 xi{xs[i]};
        vec3 dv{0.0, 0.0, 0.0};
        vec3 dx{0.0, 0.0, 0.0};
        for (int k = 0; k < 3; k++)
            if (abs(xi(k)) > globals.config["bound"] && vi(k) * xi(k) > 0.0)
            {
                dv(k) += -vi(k) * (1.0 + eps);

                if (xi(k) > 0.0)
                    dx(k) += globals.config["bound"] - xi(k);
                else
                    dx(k) += -globals.config["bound"] - xi(k);
            }
        // b.segment<3>(i * 3) += M * dv;
        vs[i] += dv;
        xs[i] += dx;
#endif
    }
}

void extract_edges(vector<Edge> &edges, const vector<unsigned> indices)
{
    static set<array<unsigned, 2>> e;
    e.clear();
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

void gen_empty_sm(
    int n_mass,
    // vector<array<int, 4>>& idx,
    // vector<array<int, 4>>& eidx,
    vector<Edge> &edges,
    SparseMatrix<double> &sparse_hess,
    map<array<int, 2>, int> &lut)
{
    static const int n_submat = 3;
    static set<array<int, 2>> cset;
    cset.clear();
    const auto insert = [&](int a, int b)
    {
        cset.insert({a, b});
        cset.insert({b, a});
    };
    for (auto &e : edges)
        insert(e.i, e.j);
    for (int i = 0; i < n_mass; i++)
        cset.insert({i, i});

    auto old = cset.begin();
    auto old_col = (*old)[0];
    for (auto it = cset.begin();; it++)
    {
        if (it == cset.end() || (*it)[0] != old_col)
        {
            for (int c = 0; c < n_submat; c++)
            {

                auto cc = c + old_col * n_submat;
                sparse_hess.startVec(cc);
                int k = 0;
                for (auto kt = old; kt != it; ++kt)
                {
                    lut[{(*kt)[1], (*kt)[0]}] = k++;
                    auto row = (*kt)[1];
                    for (int r = 0; r < n_submat; r++)
                    {
                        auto rr = row * n_submat + r;
                        sparse_hess.insertBack(rr, cc) = 0.0;
                    }
                }
            }
            if (it == cset.end())
                break;
            old = it;
            old_col = (*it)[0];
        }
    }
    sparse_hess.makeCompressed();
    // spdlog::info("\nsparse matrix : #non-zero blocks = {}", cset.size());
    // cout << sparse_hess;
}
