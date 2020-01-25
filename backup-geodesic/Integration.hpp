/**
 *  @2019 Kazuki IWANAGA, the University of Tsukuba.
 *  All rights reserved.
 */

#ifndef ___INTEGRATION___
#define ___INTEGRATION___

// std::array<double, 8> tmp, k1, k2, k3, k4;

struct bh_system {
    using state = std::array<double, 8>;

};

struct BH_System{
    using state = std::array<double, 8>;

    double t;
    double r;
    double theta;
    double phi;
    double p_t;
    double p_r;
    double p_theta;
    double p_phi;

    double dot_t;
    double dot_r;
    double dot_theta;
    double dot_phi;
    double dot_p_t;
    double dot_p_r;
    double dot_p_theta;
    double dot_p_phi;

    double Sigma;
    double Delta;
    double E;
    double L;
    double Q;
    double kappa;

    BH_System(const std::array<double, 8>& ray)
        : E(-ray[4]),
          L(ray[7]),
          Q(ray[6]*ray[6] + cos(ray[2])*cos(ray[2])
            * (L*L / (sin(ray[2])*sin(ray[2])) - BH_A*BH_A * E*E)),
          kappa(Q + L*L + BH_A*BH_A * E*E) {}
    
    void operator()(const state& x, state& dx, double lamda) {
        t       = x[0];
        r       = x[1];
        theta   = x[2];
        phi     = x[3];
        p_t     = x[4];
        p_r     = x[5];
        p_theta = x[6];
        p_phi   = x[7];

        Sigma = calcSigma(r, theta);
        Delta = calcDelta(r);

        dot_t = E + (2.0*r*E*(r*r + BH_A*BH_A) - 2.0*BH_A*r*L) / (Sigma*Delta);
        dot_r = p_r * Delta / Sigma;
        dot_theta = p_theta / Sigma;
        dot_phi = (2.0*BH_A*r*E + (Sigma - 2.0*r)*L
                            / (sin(theta)*sin(theta))) / (Sigma * Delta);
        dot_p_t = 0.0;
        dot_p_r = ((1.0 - r)*kappa + 2.0*r*E*E*(r*r 
                            + BH_A*BH_A) - 2.0*BH_A*E*L) / (Sigma * Delta)
                            - 2.0*p_r*p_r*(r - 1.0) / Sigma;
        dot_p_theta = (L*L/(sin(theta)*sin(theta)*sin(theta)*sin(theta))
                                - BH_A*BH_A*E*E) * sin(theta)*cos(theta) / Sigma;
        dot_p_phi = 0.0;

        dx[0] = dot_t;
        dx[1] = dot_r;
        dx[2] = dot_theta;
        dx[3] = dot_phi;
        dx[4] = dot_p_t;
        dx[5] = dot_p_r;
        dx[6] = dot_p_theta;
        dx[7] = dot_p_phi;
    }
};

// void ODE(std::array<double, 8>& k, const std::array<double, 8>& ray) {
//     double t       = ray[0];
//     double r       = ray[1];
//     double theta   = ray[2];
//     double phi     = ray[3];
//     double p_t     = ray[4];
//     double p_r     = ray[5];
//     double p_theta = ray[6];
//     double p_phi   = ray[7];

//     double E = -p_t;
//     double L = p_phi;
//     double Q = p_theta*p_theta + cos(theta)*cos(theta)
//                 * (L*L/(sin(theta)*sin(theta)) - BH_A*BH_A*E*E);

//     double Sigma = calcSigma(r, theta);
//     double Delta = calcDelta(r);
//     double kappa = Q + L*L + BH_A*BH_A*E*E;

//     double dot_t = E + (2.0*r*E*(r*r + BH_A*BH_A) - 2.0*BH_A*r*L) / (Sigma*Delta);
//     double dot_r = p_r * Delta / Sigma;
//     double dot_theta = p_theta / Sigma;
//     double dot_phi = (2.0*BH_A*r*E + (Sigma - 2.0*r)*L
//                         / (sin(theta)*sin(theta))) / (Sigma * Delta);
//     double dot_p_t = 0.0;
//     double dot_p_r = ((1.0 - r)*kappa + 2.0*r*E*E*(r*r 
//                         + BH_A*BH_A) - 2.0*BH_A*E*L) / (Sigma * Delta)
//                         - 2.0*p_r*p_r*(r - 1.0) / Sigma;
//     double dot_p_theta = (L*L/(sin(theta)*sin(theta)*sin(theta)*sin(theta))
//                             - BH_A*BH_A*E*E) * sin(theta)*cos(theta) / Sigma;
//     double dot_p_phi = 0.0;

//     k[0] = dot_t;
//     k[1] = dot_r;
//     k[2] = dot_theta;
//     k[3] = dot_phi;
//     k[4] = dot_p_t;
//     k[5] = dot_p_r;
//     k[6] = dot_p_theta;
//     k[7] = dot_p_phi;
// }

// void integrate(std::array<double, 8>& ray) {
//     tmp = ray;
//     ODE(k1, tmp);
//     for (int i = 0; i < 8; ++i)
//         k1[i] *= STEP;
    
//     tmp = ray;
//     for (int i = 0; i < 8; ++i)
//         tmp[i] += k1[i] / 2.0;
//     ODE(k2, tmp);
//     for (int i = 0; i < 8; ++i)
//         k2[i] *= STEP/2.0;

//     tmp = ray;
//     for (int i = 0; i < 8; ++i)
//         tmp[i] += k2[i] / 2.0;
//     ODE(k3, tmp);
//     for (int i = 0; i < 8; ++i)
//         k3[i] *= STEP/2.0;

//     tmp = ray;
//     for (int i = 0; i < 8; ++i)
//         tmp[i] += k3[i];
//     ODE(k4, tmp);
//     for (int i = 0; i < 8; ++i)
//         k4[i] *= STEP;

//     for (int i = 0; i < 8; ++i)
//         ray[i] += (k1[i] + k2[i]*2.0 + k3[i]*2.0 + k4[i]) / 6.0;
// }

#endif // !___INTEGRATION___
