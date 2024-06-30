//
// Created by alex_ on 5/25/2024.
//

#ifndef SPACECRAFT_H
#define SPACECRAFT_H
#include "state_vector.h"

class spacecraft
{
  state_vector m_state;  // TODO: probably should be a pointer
  double m_mass;

public:
  spacecraft(state_vector& state, const double& mass): m_state(state), m_mass(mass) {}

  auto get_state() -> state_vector&
  {
    return m_state;
  }

  auto get_mass() -> double
  {
    return m_mass;
  }

};

#endif //SPACECRAFT_H
