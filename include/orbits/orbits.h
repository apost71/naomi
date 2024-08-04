//
// Created by alex_ on 6/11/2024.
//

#ifndef ORBITS_H
#define ORBITS_H
#include <armadillo>

#include "constants.h"
#include "naomi.h"
namespace naomi::orbits
{
inline vector_type get_circular_orbit(const arma::vec3& initial_position, const double mu = constants::EARTH_MU)
{
  const double vn = sqrt(mu / norm(initial_position));
  const arma::vec k_hat {0.0, 0.0, 1.0};
  const arma::vec h_dir = arma::normalise(
      cross(initial_position, cross(k_hat, initial_position)));
  const arma::vec v_dir = arma::normalise(arma::cross(h_dir, initial_position));
  const arma::vec v = v_dir*vn;
  return join_cols(initial_position, v);
}

/**
 * Compute the velocity of an object given the radius, semi-major axis, and mu
 *
 * @param radius The radius at the current location in orbit in meters
 * @param sma The semi-major axis of the orbit in meters
 * @param mu Optional gravitational parameter for the central body, defaults to
 * the gravitational parameter for Earth
 * @return The velocity of the object
 */
[[nodiscard]] inline double vis_viva(const double radius,
                              const double sma,
                              const double mu = constants::EARTH_MU)
{
  return sqrt(mu * (2/radius - 1/sma));
}
}
#endif //ORBITS_H
