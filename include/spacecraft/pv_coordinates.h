//
// Created by alex on 7/23/2024.
//

#ifndef PV_COORDINATES_H
#define PV_COORDINATES_H

#include <armadillo>

#include <fmt/format.h>

class pv_coordinates
{
  arma::vec3 _pos;
  arma::vec3 _vel;
  arma::vec3 _acc;

public:
  explicit pv_coordinates(const naomi::state_type& state)
  {
    if (state.size() == 6) {
      _pos = state(arma::span(0, 2));
      _vel = state(arma::span(3, 5));
    } else if (state.size() == 9) {
      _pos = state(arma::span(0, 2));
      _vel = state(arma::span(3, 5));
      _acc = state(arma::span(6, 8));
    } else {
      throw std::runtime_error(
        fmt::format("PV state vector must have size 6 or 9 but was {}",
          state.size()));
    }
  }

  explicit pv_coordinates(const arma::vec9& state):
    _pos(state(arma::span(0, 2)))
  , _vel(state(arma::span(3, 5)))
  , _acc(state(arma::span(6, 8))){}


  auto get_position() -> arma::vec3&
  {
    return _pos;
  }

  auto get_velocity() -> arma::vec3&
  {
    return _vel;
  }

  auto get_acceleration() -> arma::vec3
  {
    return _acc;
  }

  auto to_vec() const -> arma::vec9
  {
    const auto pv = join_cols(_pos, _vel);
    auto pva = join_cols(pv, _acc);
    return pva;
  }
};

#endif //PV_COORDINATES_H