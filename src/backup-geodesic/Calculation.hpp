/**
 *  @2019 Kazuki IWANAGA, the University of Tsukuba.
 *  All rights reserved.
 */

#ifndef ___CALCULATION___
#define ___CALCULATION___

double calcSigma(double r, double theta) {
    return r * r + BH_A * BH_A * cos(theta) * cos(theta);
}

double calcDelta(double r) {
    return r * r - 2.0 * BH_M * r + BH_A * BH_A;
}

double calcE(double r, double theta, double dot_r, double dot_theta, double dot_phi) {
    double Sigma = calcSigma(r, theta);
    double Delta = calcDelta(r);
    double E_sq = (Sigma - 2.0 * r) * (Sigma * dot_r * dot_r 
                    + Sigma * Delta * dot_theta * dot_theta) / (Sigma * Delta)
                    + Delta * dot_phi * dot_phi * sin(theta) * sin(theta);
    return sqrt(E_sq);
}

double calcL(double r, double theta, double dot_r, double dot_theta, double dot_phi) {
    double Sigma = calcSigma(r, theta);
    double Delta = calcDelta(r);
    double E = calcE(r, theta, dot_r, dot_theta, dot_phi);
    return (sin(theta) * sin(theta) * (Sigma * Delta * dot_phi 
                        - 2.0 * BH_A * r * E)) / (Sigma - 2.0 * r);
}

double calcQ(double r, double theta, double dot_r, double dot_theta, double dot_phi) {
    double Sigma = calcSigma(r, theta);
    double Delta = calcDelta(r);
    double E = calcE(r, theta, dot_r, dot_theta, dot_phi);
    double L = calcL(r, theta, dot_r, dot_theta, dot_phi);
    return Sigma * Sigma * dot_theta * dot_theta + cos(theta) * cos(theta)
            * (L * L / (sin(theta) * sin(theta)) - BH_A * BH_A * E * E);
}

#endif // !___CALCULATION___
