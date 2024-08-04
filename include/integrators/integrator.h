//
// Created by alex_ on 6/10/2024.
//

#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <naomi.h>
#include <systems/two_body.h>

#include "propagators/event_detector.h"
#include "forces/force_model.h"

namespace naomi::numeric
{
using namespace events;
using namespace forces;

typedef std::function<void(const vector_type&, vector_type&, double)> system_t;

template< class Stepper>
class integrator
{
public:
  integrator(const integrator& other)
      : m_stepper(other.m_stepper)
  {
  }
  integrator(integrator&& other) noexcept
      : m_stepper(std::move(other.m_stepper))
  {
  }
  integrator& operator=(const integrator& other)
  {
    if (this == &other)
      return *this;
    m_stepper = other.m_stepper;
    return *this;
  }
  integrator& operator=(integrator&& other) noexcept
  {
    if (this == &other)
      return *this;
    m_stepper = std::move(other.m_stepper);
    return *this;
  }

private:
  Stepper m_stepper;

public:
  ~ integrator() = default;
  integrator() = default;

  std::pair<double, vector_type> find_event_time(const system_t& system, double start_time, double end_time, std::shared_ptr<event_detector> e, state_and_time_type state, double step_size)
  {
    auto stepper = make_dense_output(1.0e-6, 1.0e-6, Stepper());
    vector_type s = state.first;
    stepper.initialize(s, start_time, end_time-start_time);
    stepper.do_step(system);
    double mid_time;
    vector_type next_state = s;
    while(std::abs(end_time - start_time) > 1e-6) {
      mid_time = 0.5 * (start_time + end_time);  // get the mid point time
      stepper.calc_state(mid_time, next_state); // obtain the corresponding state
      state_and_time_type prev = {s, start_time};
      state_and_time_type curr = {next_state, mid_time};
      if ((*e)(prev, curr))
        end_time = mid_time;  // condition changer lies before midpoint
      else {
        start_time = mid_time;  // condition changer lies after midpoint
        s = next_state;
      }
    }
    // we found the interval of size eps, take it's midpoint as final guess
    mid_time = 0.5 * (start_time + end_time);
    stepper.calc_state(mid_time, s);
    return {mid_time, s};
  }

  double integrate(const system_t& system, vector_type& state, double start_time, double end_time, double step_size)
  {
    integrate_adaptive(make_controlled( 1.0e-6 , 1.0e-6, m_stepper), system, state, start_time, end_time, step_size);
    // integrate_const(m_stepper, system, state, start_time, end_time, step_size);
    return end_time;
  }
};
}


#endif //INTEGRATOR_H
