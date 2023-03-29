#include "dynamics.h"
#include <spdlog/spdlog.h>
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
void compute_b(VectorXd &b)
{
    auto &vs{globals.mesh->mass_v}, &xs{globals.mesh->mass_x};

    int n_mass = vs.size();
    b.setZero(3 * n_mass);
    for (int i = 0; i < n_mass; i++)
    {
        b.segment<3>(i * 3) = M * vs[i];
    }
    compute_force(b);
}

void compute_force(VectorXd &b)
{
    auto &xs{globals.mesh->mass_x};
    for (auto &e : globals.mesh->edges)
    {
        int i = e.i, j = e.j;
        vec3 xji{xs[j] - xs[i]};
        vec3 fi = ks * (xji - e.l0 * (xji).normalized());
        b.segment<3>(3 * i) += fi;
        b.segment<3>(3 * j) -= fi;
    }
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

void compute_single_spring_K(Edge &e)
{
    auto &xs{globals.mesh->mass_x};

    vec3 xji{xs[e.j] - xs[e.i]};

    mat3 K = (xji) * (xji).transpose();
    double l = xji.norm();
    // K = K * (kd + ks * l0)
    K = K * ks * e.l0 / (l * l * l) + ks * mat3::Identity(3, 3) * (e.l0 - l) / l;
}