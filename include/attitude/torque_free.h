//
// Created by alex on 7/10/2024.
//

#ifndef TORQUE_FREE_H
#define TORQUE_FREE_H

#include "attitude_provider.h"
#include "frames/transforms.h"
#include "math/vector_utils.h"

namespace naomi::attitude
{
using namespace math;

class torque_free_attitude final : public attitude_provider
{
  static arma::vec3 compute_angular_velocity(const state_type& state)
  {
    const auto t = eci2ric(state(arma::span(0, 5)));
    const state_type r = state(arma::span(0, 2));
    const state_type v = state(arma::span(3, 5));
    const state_type w = t * cross(r, v) / dot(r, r);
    return w;
  }
public:
  explicit torque_free_attitude(const state_type& attitude, const arma::mat33& inertia_matrix): attitude_provider(attitude, inertia_matrix){}

  [[nodiscard]] state_type get_derivative(
      const state_type& state) const override
  {
    const auto w = compute_angular_velocity(state);
    const state_type w4 = {0, w[0], w[1], w[2]};
    const quaternion_type q_state = state(arma::span(6, 9));
    state_type res = 0.5 * q_skew(q_state) * w4;
    return res;
    // return {
    //   -(m_inertia_matrix(2, 2) - m_inertia_matrix(1, 1))*w[1]*w[2] / m_inertia_matrix(0, 0),
    //   -(m_inertia_matrix(0, 0) - m_inertia_matrix(2, 2))*w[2]*w[0] / m_inertia_matrix(1, 1),
    //   -(m_inertia_matrix(1, 1) - m_inertia_matrix(0, 0))*w[0]*w[1] / m_inertia_matrix(2, 2),
    // };
  }

  quaternion_type get_rotation() override
  {
    quaternion_type res;
    return res;
  }

  state_type get_angular_momentum() override
  {
    return {1, 2, 3};
  }
  state_type get_angular_velocity() override
  {
    return {1, 2, 3};
  }
};
}
#endif //TORQUE_FREE_H
