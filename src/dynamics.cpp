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

            int k = i / 2;
            dv(k) += -vi(k) * (1.0 + eps);
            if (xi(k) > 0.0)
                dx(k) += bound - xi(k);
            else
                dx(k) += -bound - xi(k);
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
