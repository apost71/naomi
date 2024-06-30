//
// Created by alex_ on 5/24/2024.
//

#ifndef COWELL_PROPAGATOR_H
#define COWELL_PROPAGATOR_H

#include "integrators/integrator.h"
#include "spacecraft/spacecraft.h"
#include "systems/two_body.h"

template <typename TStepper>
class numerical_propagator
{

  integrator<TStepper> m_integrator;
  two_body_system& m_system;
  double m_t = 0.0;
public:
  ~numerical_propagator() = default;
  explicit numerical_propagator(two_body_system& system)
      : m_system(system){}

  numerical_propagator(const numerical_propagator&) = delete;
  numerical_propagator(numerical_propagator&&) = delete;
  auto operator=(const numerical_propagator&) -> numerical_propagator& = delete;
  auto operator=(const numerical_propagator&&) -> numerical_propagator& = delete;

  void propagate_by(double dt)
  {
    m_integrator.integrate( m_system, m_system.get_state() , m_t , m_t + dt , 0.1 );
    m_t += dt;
  };
};

#endif //COWELL_PROPAGATOR_H
