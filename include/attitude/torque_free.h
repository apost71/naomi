//
// Created by alex on 7/10/2024.
//

#ifndef TORQUE_FREE_H
#define TORQUE_FREE_H

#include "attitude_law.h"
#include "frames/transforms.h"
#include "math/vector_utils.h"

namespace naomi::attitude
{
using namespace math;

class torque_free_attitude final : public attitude_law
{
public:

  [[nodiscard]] state_type get_derivative(
      const state_type& state) const override
  {
    const auto t = eci2ric(state(arma::span(0, 5)));
    const state_type r = state(arma::span(0, 2));
    const state_type v = state(arma::span(3, 5));
    const state_type w = t * cross(r, v) / dot(r, r);
    const state_type w4 = {0, w[0], w[1], w[2]};
    const quaternion_type q_state = state(arma::span(6, 9));
    state_type res = (0.5 * q_skew(q_state) * w4);
    return res;
  }

  [[nodiscard]] std::size_t get_size() const override
  {
    return 4;
  }
};
}
#endif //TORQUE_FREE_H
