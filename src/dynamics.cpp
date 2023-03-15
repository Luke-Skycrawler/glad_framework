#include "dynamics.h"
#include <spdlog/spdlog.h>
mat3 RigidBody::compute_K(const vec3 &r, int t){
    auto & R {S[t].R};
    mat3 Rr {skew(R * r)};
    mat3 J_inv {R * J0.inverse() * R.transpose()};
    return M_inv * mat3::Identity(3, 3) - Rr * J_inv * Rr;
}

StateVector RigidBody::dSdt(const vec3 &f, int t){
    auto & R{S[t].R};
    auto &pdot = f;
    auto & xdot = M_inv * S[t].p;
    auto & Rdot = skew(R * J0.inverse() * R.transpose() * S[t].L) * R;
    vec3 Ldot(0.0, 0.0, 0.0);
    return {xdot, pdot, Ldot, Rdot};
}

void RigidBody::step(int ts) {
    auto Sdot = dSdt(vec3(0.0, 0.0, 0.0), 0);
    
    S[1] = S[0] + Sdot * dt;
    col_set.resize(0);
    auto R = S[1].R;

    for (int i = 0; i < n_vertices;i++) {
        vec3 &r {vertices[i]}, xi {R * r + S[1].x};
        if ((xi.array().abs() > bound).any()) {
            col_set.push_back(i);
        }
    }
    if (col_set.size()) spdlog::info("collsiion set size = {}", col_set.size());

    vec3 dp(0.0, 0.0, 0.0), dL(0.0, 0.0, 0.0), dx(0.0, 0.0, 0.0);
    vec3 r(0.0, 0.0, 0.0);
    
    for (int I:col_set) {
    //     vec3 &rr {vertices[I]};
    //     r +=rr;
    // }
    // if (col_set.size()) {
    //     r /= col_set.size();
    //     {
        vec3 r{vertices[I]};
            vec3 Rr{ R * r },
                xi{ R * r + S[1].x };
            auto K = compute_K(r, 1);

            vec3 omega = R * J0.inverse() * R.transpose() * S[1].L;
            vec3 vi = M_inv * S[1].p + skew(omega) * Rr;

            vec3 dv(0.0, 0.0, 0.0);

            for (int k = 0; k < 3; k++) {
                if (abs(xi(k)) > bound) {
                    dv(k) += -vi(k) * (1.0 + eps);
                    if (xi(k) > 0.0) dx(k) += bound - xi(k);
                    else dx(k) += -bound - xi(k);
                }
            }
            vec3 deltap = K.inverse() * dv;
            dp += deltap;
            dL += skew(Rr) * deltap;

            vec3 omega1 = R * J0.inverse() * R.transpose() * (S[1].L + dL);
            vec3 vi1 = M_inv * (S[1].p + dp) + skew(omega1) * Rr;

            for (int k = 0; k < 3; k++) if (dv(k) != 0.0) {
                spdlog::info("dim = {}, vi0 = {:.6f}, vi1 = {:.6f}", k, vi(k), vi1(k));
            }

        }
    //}
    if (col_set.size()) {
        spdlog::info("ts = {}", ts);
        col_set.resize(0);
    }


    S[1].p += dp;
    S[1].L += dL;
    S[0] = S[1];
    S[0].x += dx;
}



