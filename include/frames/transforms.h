//
// Created by alex_ on 6/23/2024.
//

#ifndef TRANSFORMS_H
#define TRANSFORMS_H

#include <armadillo>

#include "spacecraft/state_vector.h"

arma::mat33 eci2ric(arma::vec3 r, arma::vec3 v)
{
  auto rn = normalise(r);
  auto vn = normalise(v);
  auto h = cross(r, v);
  auto c = normalise(h);
  auto i = cross(c, rn);
  return arma::join_rows(rn, i, c);
}

arma::mat33 eci2ric(state_vector& sv)
{
  return eci2ric(sv.get_position(), sv.get_velocity());
}

#endif //TRANSFORMS_H
