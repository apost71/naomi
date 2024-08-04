//
// Created by alex_ on 7/31/2024.
//

#ifndef TORQUE_FREE_PROVIDER_H
#define TORQUE_FREE_PROVIDER_H

#include <utility>

#include "attitude/attitude_provider.h"
#include "spacecraft/state_provider.h"
#include "torque_free.h"

class torque_free_attitude_provider:
  public attitude::attitude_provider,
  public integrated_provider
{

  quaternion_type _q = {1, 0, 0, 0};
  state_type _w = {0, 0, 0};
  state_type _wdot = {0, 0, 0};
  arma::mat33 _inertia_matrix;
  std::shared_ptr<forces::equations_of_motion> _eoms;

  static state_type compute_angular_velocity(const state_type& state)
  {
    const auto t = eci2ric(state(arma::span(0, 5)));
    const state_type r = state(arma::span(0, 2));
    const state_type v = state(arma::span(3, 5));
    const state_type w = t * cross(r, v) / dot(r, r);
    return w;
  }

  static state_type compute_angular_velocity(const pv_coordinates& state)
  {
    const state_type r = state.get_position();
    const state_type v = state.get_velocity();
    const state_type pv = join_cols(r, v);
    const auto t = eci2ric(pv);
    const state_type w = t * cross(r, v) / dot(r, r);
    return w;
  }

public:
  torque_free_attitude_provider(
      const arma::mat33& inertia_matrix,
      const quaternion_type& q,
      const pv_coordinates& pv):
      _q(q)
      , _w(compute_angular_velocity(pv))
      , _inertia_matrix(inertia_matrix)
      , _eoms(std::make_shared<attitude::torque_free_eoms>(inertia_matrix))
  {
  }

  torque_free_attitude_provider(
    const arma::mat33& inertia_matrix,
    const quaternion_type& q,
    const state_type& pv):
    _q(q)
    , _w(compute_angular_velocity(pv))
    , _inertia_matrix(inertia_matrix)
    , _eoms(std::make_shared<attitude::torque_free_eoms>(inertia_matrix))
  {
  }

  quaternion_type get_rotation() override
  {
    return _q;
  }
  state_type get_angular_momentum() override
  {
    return {0, 0, 1};
  }
  state_type get_angular_velocity() override
  {
    return _w;
  }
  std::shared_ptr<forces::equations_of_motion> get_eoms() override
  {
    return _eoms;
  }

  std::size_t get_size() override
  {
    return 10;
  }

  [[nodiscard]] state_type get_integrated_state() override
  {
    state_type state(10);
    state(arma::span(0, 3)) = _q;
    state(arma::span(4, 6)) = _w;
    state(arma::span(7, 9)) = _wdot;
    return state;
  }

  void set_integrated_state(const state_type& state) override
  {
    _q = state(arma::span(0, 3));
    _w = state(arma::span(4, 6));
    _wdot = state(arma::span(7, 9));
  }

  void apply_force(const arma::vec3& forces) override {}
};

#endif //TORQUE_FREE_PROVIDER_H
