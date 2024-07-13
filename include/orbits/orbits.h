//
// Created by alex_ on 6/11/2024.
//

#ifndef ORBITS_H
#define ORBITS_H
#include <armadillo>

#include "constants.h"
#include "spacecraft/state_vector.h"

namespace naomi::orbits
{
inline state_type get_circular_orbit(const arma::vec3& initial_position, const double mu = constants::EARTH_MU)
{
  const double vn = sqrt(mu / norm(initial_position));
  const arma::vec k_hat {0.0, 0.0, 1.0};
  const arma::vec h_dir = arma::normalise(
      arma::cross(initial_position, cross(k_hat, initial_position)));
  const arma::vec v_dir = arma::normalise(arma::cross(h_dir, initial_position));
  const arma::vec v = v_dir*vn;
  return arma::join_cols(initial_position, v);
}
}
#endif //ORBITS_H
