//
// Created by alex_ on 6/2/2024.
//

#ifndef EARTH_H
#define EARTH_H
#include "celestial_body.h"
#include "constants.h"
#include <symengine/expression.h>
#include <symengine/lambda_double.h>

using SymEngine::Expression;

#ifndef EARTH_SOI
#define EARTH_SOI 9.24 * 10e5 * 1000.0 // meters
#endif

#ifndef EARTH_J2
#define EARTH_J2  1082.63 * 10e-6
#endif

#ifndef EARTH_EQ_RADIUS
#define EARTH_EQ_RADIUS 6378.1 * 1000.0 // meters
#endif


class earth : public celestial_body
{
  double earth_j2 = 1082.63 * 10e-6;
  double earth_radius = 6378.1 * 1000.0; // meters
  double a_j2 = 0.5 * earth_j2 * pow(earth_radius, 2);

public:
  earth()
    : celestial_body(constants::EARTH_MU, EARTH_SOI, EARTH_EQ_RADIUS, {EARTH_J2}){}

  auto get_potential_partial_derivative(
      Eigen::Vector3d position)
    -> Eigen::Vector3d override
  {
    double x = position[0];
    double y = position[1];
    double z = position[2];
    Eigen::Vector3d result;

    const double r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
    double u_x = m_mu
        * (-x / pow(r, 3)
          + a_j2 * (15 * x * pow(z, 2) / pow(r, 7) - 3 * x / pow(r, 5)));
    double u_y = m_mu
        * (-y / pow(r, 3)
          + a_j2 * (15 * y * pow(z, 2) / pow(r, 7) - 3 * y / pow(r, 5)));
    double u_z = m_mu
        * (-z / pow(r, 3)
          + a_j2 * (15 * z * pow(z, 2) / pow(r, 7) - 9 * z / pow(r, 5)));

    result << u_x, u_y, u_z;
    return result;
  }

  auto get_potential_partial_derivative(arma::vec position)
    -> arma::vec override
  {
    double x = position[0];
    double y = position[1];
    double z = position[2];

    const double r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
    double u_x = m_mu
        * (-x / pow(r, 3)
          + a_j2 * (15 * x * pow(z, 2) / pow(r, 7) - 3 * x / pow(r, 5)));
    double u_y = m_mu
        * (-y / pow(r, 3)
          + a_j2 * (15 * y * pow(z, 2) / pow(r, 7) - 3 * y / pow(r, 5)));
    double u_z = m_mu
        * (-z / pow(r, 3)
          + a_j2 * (15 * z * pow(z, 2) / pow(r, 7) - 9 * z / pow(r, 5)));

    arma::vec result = {u_x, u_y, u_z};
    return result;
  }
};

#endif //EARTH_H