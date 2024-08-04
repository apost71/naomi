//
// Created by alex on 7/10/2024.
//

#ifndef TORQUE_FREE_H
#define TORQUE_FREE_H

#include "frames/transforms.h"
#include "math/vector_utils.h"
#include "spacecraft/state_provider.h"

namespace naomi::attitude
{
using namespace math;

class torque_free_eoms final :
  public forces::equations_of_motion
{

  arma::mat33 _inertia_matrix;

public:
  explicit torque_free_eoms(const arma::mat33& inertia_matrix): _inertia_matrix(inertia_matrix){}
  ~torque_free_eoms() override = default;
  [[nodiscard]] state_type get_derivative(const state_type& state, double t) const override
  {
    const auto q = state(arma::span(0, 3));
    const auto w = state(arma::span(4, 6));
    const state_type w4 = {0, w[0], w[1], w[2]};
    const state_type q_dot = 0.5 * q_skew(q) * w4;
    const arma::vec3 w_dot = {
      -(_inertia_matrix(2, 2) - _inertia_matrix(1, 1))*w[1]*w[2] / _inertia_matrix(0, 0),
      -(_inertia_matrix(0, 0) - _inertia_matrix(2, 2))*w[2]*w[0] / _inertia_matrix(1, 1),
      -(_inertia_matrix(1, 1) - _inertia_matrix(0, 0))*w[0]*w[1] / _inertia_matrix(2, 2),
    };
    state_type res(10);
    res(arma::span(0, 3)) = q_dot;
    res(arma::span(4, 6)) = w_dot;
    return res;
  }

};
}
#endif //TORQUE_FREE_H
