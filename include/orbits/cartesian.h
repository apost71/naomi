//
// Created by alex_ on 6/15/2024.
//

#ifndef CARTESIAN_ORBIT_H
#define CARTESIAN_ORBIT_H

#include <armadillo>

class cartesian_orbit
{
  arma::vec3 m_pos;
  arma::vec3 m_vel;

public:
  cartesian_orbit(const arma::vec3& pos, const arma::vec3& vel): m_pos(pos), m_vel(vel){}

  auto get_position() -> arma::vec3&
  {
    return m_pos;
  }

  auto get_velocity()  -> arma::vec3&
  {
    return m_vel;
  }
};

#endif //CARTESIAN_ORBIT_H
