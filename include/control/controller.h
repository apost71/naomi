//
// Created by alex on 7/4/2024.
//

#ifndef SPACECRAFT_CONTROLLER_H
#define SPACECRAFT_CONTROLLER_H
#include "naomi.h"
#include "spacecraft/spacecraft_state.h"

namespace naomi::control
{
using namespace naomi;

class control_input
{
  arma::vec3 m_position_control_input;
  quaternion_type m_attitude_control_input;

public:
  control_input(const arma::vec3& position_control_input, const quaternion_type& attitude_control_input):
    m_position_control_input(position_control_input), m_attitude_control_input(attitude_control_input){}

};

class controller
{
public:
  virtual ~controller() = default;
  virtual control_input get_control_input(const vector_type& state, const vector_type& attitude,
                                             double t) = 0;
  virtual void initialize(const vector_type& state, const vector_type& attitude, double t) = 0;

  /**
   * @brief
   * @param state Spacecraft state in a stacked vector with state in the
   *  inertial frame and attitude represented as a quaternion
   * @param attitude The spacecraft attitude in quaternion form
   * @param t current time
   * @return control input
   */
  virtual vector_type get_desired_state(const vector_type& state, const vector_type& attitude, double t) = 0;

};

// class spacecraft_subsystem: simulation_component<spacecraft_state>
// {
//   std::shared_ptr<controller> m_controller = nullptr;
//   std::shared_ptr<attitude::additional_state_provider> m_state_provider = nullptr;
//
// public:
//   ~spacecraft_subsystem() override = default;
// };

}
#endif //SPACECRAFT_CONTROLLER_H
