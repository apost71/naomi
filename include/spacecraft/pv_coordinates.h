//
// Created by alex on 7/23/2024.
//

#ifndef PV_COORDINATES_H
#define PV_COORDINATES_H

#include <armadillo>

class pv_coordinates
{
  arma::vec3 m_pos;
  arma::vec3 m_vel;
  arma::vec3 m_acc;

public:
  explicit pv_coordinates(const arma::vec6& state):
    m_pos(state(arma::span(0, 2))), m_vel(state(arma::span(3, 5))){}

  auto get_position() -> arma::vec3&
  {
    return m_pos;
  }

  auto get_velocity() -> arma::vec3&
  {
    return m_vel;
  }

  auto get_acceleration() -> arma::vec3
  {
    return m_acc;
  }
};

#endif //PV_COORDINATES_H
