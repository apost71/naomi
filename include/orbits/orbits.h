//
// Created by alex_ on 6/11/2024.
//

#ifndef ORBITS_H
#define ORBITS_H
#include <armadillo>

#include "constants.h"
#include "spacecraft/state_vector.h"

inline state_vector get_circular_orbit(const arma::vec3& initial_position, const double mu = constants::EARTH_MU)
{
  double vn = sqrt(mu/norm(initial_position));
  arma::vec k_hat {0.0, 0.0, 1.0};
  arma::vec h_dir = arma::normalise(arma::cross(initial_position, cross(k_hat, initial_position)));
  arma::vec v_dir = arma::normalise(arma::cross(h_dir, initial_position));
  arma::vec v = v_dir*vn;
  state_vector sv(initial_position, v);
  return sv;
}

#endif //ORBITS_H
