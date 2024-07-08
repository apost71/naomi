//
// Created by alex_ on 6/10/2024.
//

#ifndef CONSTANTS_H
#define CONSTANTS_H

#ifndef GRAVITATIONAL_CONSTANT
#define GRAVITATIONAL_CONSTANT 6.67428e-11L
#endif

// #ifndef EARTH_MU
// #define EARTH_MU (3.986004418*10e14)
// // #define EARTH_MU (3.9800000*10e14)

namespace naomi::constants
{
constexpr double EARTH_MU = 3.986004418*1e14;
constexpr double EARTH_MU_KM = 3.986004418*1e5;
const arma::vec3 PLUS_I({1, 0, 0});
const arma::vec3 PLUS_J({0, 1, 0});
const arma::vec3 PLUS_K({0, 0, 1});
const arma::vec3 MINUS_I({-1, 0, 0});
const arma::vec3 MINUS_J({0, -1, 0});
const arma::vec3 MINUS_K({0, 0, -1});

}



// #define EARTH_MU 398600.4418
#endif
