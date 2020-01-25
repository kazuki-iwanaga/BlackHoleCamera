/**
 *  @2019 Kazuki IWANAGA, the University of Tsukuba.
 *  All rights reserved.
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <array>
#include <string>
#include <boost/numeric/odeint.hpp>

constexpr double PI = 4.0 * atan(1.0);

constexpr double BH_M = 1.0;
constexpr double BH_R = 1.0;
constexpr double BH_A = 0.0;

constexpr double DISTANCE = 100.0;
constexpr double ELEVATION = 0.0;
constexpr double AZIMUTH = 0.0;

constexpr double SCREEN_X = -5.0;
constexpr double SCREEN_Y = 0.0;

constexpr double STEP = 0.01;

#include "Calculation.hpp"
#include "Transformation.hpp"
#include "Initialization.hpp"
#include "Integration.hpp"
#include "Observation.hpp"

int main() {
    std::array<double, 8> ray;
    std::array<double, 3> pixel = {SCREEN_X, SCREEN_Y, 0.0};

    // printPixel(pixel);
    Screen2Cartesian(pixel);
    // printPixel(pixel);

    Cartesian2BoyerLindquist(pixel);
    // printPixel(pixel);

    // BoyerLindquist2Cartesian(pixel);
    // printPixel(pixel);

    // printRay(ray);
    Initialize(ray, pixel);
    // printRay(ray);
    
    BH_System System(ray);
    BH_System::state State = ray;
    // printRay(State);

    CSV_Observer Observer("output.csv");

    boost::numeric::odeint::runge_kutta_dopri5<BH_System::state> Stepper;

    // auto Stepper = boost::numeric::odeint::make_controlled<
    //     boost::numeric::odeint::runge_kutta_dopri5<BH_System::state>
    //     >(0.01, 0.01);

    boost::numeric::odeint::integrate_const(
        Stepper, System, State, 0.0, 1000.0, 0.01, std::ref(Observer)
    );
}
