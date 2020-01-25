/**
 *  @2019 Kazuki IWANAGA, the University of Tsukuba.
 *  All rights reserved.
 */

#ifndef ___INITIALIZATION___
#define ___INITIALIZATION___

void printRay(std::array<double, 8>& ray) {
    std::cout << "Ray  --->";
    for (int i = 0; i < 8; ++i)
        std::cout << "  [" << i << "] = " << ray[i];
    std::cout << std::endl;
}

void Initialize(std::array<double, 8>& ray, std::array<double, 3>& pixel) {
    double r     = pixel[0];
    double theta = pixel[1];
    double phi   = pixel[2];
    double theta_screen = PI / 2.0 - degree2radian(ELEVATION);
    double phi_screen   = degree2radian(AZIMUTH);

    double RR  = sqrt(r*r + BH_A*BH_A);
    double PHI = phi - phi_screen;

    double Sigma = calcSigma(r, theta);
    double Delta = calcDelta(r);

    double dot_r = -(r * RR * sin(theta) * sin(theta_screen) * cos(PHI)
                        + RR*RR * cos(theta) * cos(theta_screen)) / Sigma;
    double dot_theta = (r * sin(theta) * cos(theta_screen)
                        - RR * cos(theta) * sin(theta_screen) * cos(PHI))/ Sigma;
    double dot_phi = sin(theta_screen) * sin(PHI) / (RR * sin(theta));

    double E = calcE(r, theta, dot_r, dot_theta, dot_phi);
    double L = calcL(r, theta, dot_r, dot_theta, dot_phi);
    double Q = calcQ(r, theta, dot_r, dot_theta, dot_phi);

    ray[0] = 0.0;
    ray[1] = r;
    ray[2] = theta;
    ray[3] = phi;
    ray[4] = -E;
    ray[5] = dot_r * Sigma / Delta;
    ray[6] = dot_theta * Sigma;
    ray[7] = L;
}

#endif // !___INITIALIZATION___
