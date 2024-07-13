//
// Created by alex_ on 6/23/2024.
//

#ifndef TRANSFORMS_H
#define TRANSFORMS_H

#include <armadillo>

#include "spacecraft/state_vector.h"

inline arma::mat33 eci2ric(arma::vec3 r, arma::vec3 v)
{
  auto rn = normalise(r);
  auto vn = normalise(v);
  auto h = cross(r, v);
  auto c = normalise(h);
  auto i = normalise(cross(c, rn));
  return arma::join_rows(rn, i, c);
}

inline arma::mat33 eci2ric(const state_type& sv)
{
  return eci2ric(sv(arma::span(0, 2)), sv(arma::span(3, 5)));
}

#endif //TRANSFORMS_H
