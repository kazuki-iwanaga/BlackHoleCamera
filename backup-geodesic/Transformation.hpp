/**
 *  @2019 Kazuki IWANAGA, the University of Tsukuba.
 *  All rights reserved.
 */

#ifndef ___TRANSFORMATION___
#define ___TRANSFORMATION___

double degree2radian(const double& deg) {
    return deg * PI / 180.0;
}

double radian2degree(const double& rad) {
    return rad * 180.0 / PI;
}

void printPixel(std::array<double, 3>& p) {
    std::cout << "pixel  ---->";
    for (int i = 0; i < 3; ++i)
        std::cout << "  [" << i << "] = " << p[i];
    std::cout << std::endl;
}

void Screen2Cartesian(std::array<double, 3>& p) {
    double i = p[0];
    double j = p[1];
    double distance = DISTANCE;
    double theta    = PI / 2.0 - degree2radian(ELEVATION);
    double phi      = degree2radian(AZIMUTH);

    double x = - j * cos(theta) * cos(phi) - i * sin(phi)
                        + distance * sin(theta) * cos(phi);
    double y = - j * cos(theta) * sin(phi) + i * cos(phi)
                        + distance * sin(theta) * sin(phi);
    double z = distance * cos(theta) + j * sin(theta);

    p[0] = x;
    p[1] = y;
    p[2] = z;
}

void Cartesian2BoyerLindquist(std::array<double, 3>& p) {
    double x = p[0];
    double y = p[1];
    double z = p[2];
    double useful = x*x + y*y + z*z - BH_A*BH_A;

    double r = sqrt((useful + sqrt(useful*useful
                    + 4.0 * BH_A*BH_A * z*z)) / 2.0);
    double theta = acos(z / r);
    double phi = atan2(y, x);

    p[0] = r;
    p[1] = theta;
    p[2] = phi;
}

void BoyerLindquist2Cartesian(std::array<double, 3>& p) {
    double x = sqrt(p[0] * p[0] + BH_A * BH_A) * sin(p[1]) * cos(p[2]);
    double y = sqrt(p[0] * p[0] + BH_A * BH_A) * sin(p[1]) * sin(p[2]);
    double z = p[0] * cos(p[1]);

    p[0] = x;
    p[1] = y;
    p[2] = z;
}

#endif // !___TRANSFORMATION___
