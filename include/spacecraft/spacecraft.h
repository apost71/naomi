//
// Created by alex_ on 5/25/2024.
//

#ifndef SPACECRAFT_H
#define SPACECRAFT_H
#include <utility>

#include "attitude/attitude_provider.h"
#include "body_shape.h"
#include "control/controller.h"
#include "frames/frame.h"
#include "frames/transforms.h"
#include "maneuvers/maneuver_plan.h"
#include "pv_coordinates.h"

namespace naomi {

using namespace naomi::control;
using namespace naomi::attitude;
using namespace naomi::geometry;

class spacecraft: public simulation_component<spacecraft_state>
{
  std::string m_identifier;
  pv_coordinates m_pv_coordinates;
  spacecraft_state m_state;
  quaternion_type m_attitude = {1, 0, 0, 0};
  std::vector<std::shared_ptr<controller>> m_controllers;
  std::shared_ptr<maneuvers::maneuver_plan> m_maneuver_plan;
  std::vector<std::shared_ptr<spacecraft_subsystem>> _subsystems;
  state_type m_forces;
  body_shape m_body_shape = body_shape::make_rectangle(1, 1, 1, 100);

public:
  spacecraft(std::string identifier, const state_type& state, const double& mass):
    m_identifier(std::move(identifier)),
    m_pv_coordinates(state),
    m_state(state, mass) {}

  spacecraft(std::string identifier, const state_type& state, const double& mass, const std::shared_ptr<attitude_provider>& attitude_law):
    m_identifier(std::move(identifier)),
    m_pv_coordinates(state),
    m_state(state, mass) {}

  spacecraft(std::string identifier, const state_type& state, const quaternion_type& attitude, const double& mass, const std::shared_ptr<attitude::attitude_provider>& attitude_law):
    m_identifier(std::move(identifier)),
    m_pv_coordinates(state),
    m_state(state, mass),
    m_attitude(attitude) {}


  spacecraft(std::string identifier, const state_type& state, const double& mass, const std::shared_ptr<maneuvers::maneuver_plan>& mp):
    m_identifier(std::move(identifier)), m_pv_coordinates(state), m_state(state, mass), m_maneuver_plan(mp) {}

  spacecraft(std::string identifier, const state_type& state, const quaternion_type& attitude, const double& mass, const std::shared_ptr<maneuvers::maneuver_plan>& mp):
    m_identifier(std::move(identifier)), m_pv_coordinates(state), m_state(state, mass), m_attitude(attitude), m_maneuver_plan(mp) {}

  auto get_state() const -> state_type
  {
    auto pva = m_pv_coordinates.to_vec();
    return pva;
  }

  auto get_propagation_state() const -> state_type
  {

  }

  /**
   * The current positional coordinates of the spacecraft in a given frame.
   *  TODO: Actually implement frame transformation logic.
   *
   * @param frame The desired frame of the pv coordinates
   * @return The current pv coordinates copied by value
   */
  auto get_pv_coordinates(const std::shared_ptr<frames::frame>& frame) -> pv_coordinates
  {
    return m_state.get_pv_coordinates();
  }

  auto get_attitude() -> quaternion_type&
  {
    return m_attitude;
  }

  /**
   * Get the 3x3 intertia matrix of the spacecraft.
   *
   * @return Inertia matrix for the spacecraft copied by value.
   */
  auto get_inertia_matrix() -> arma::mat33
  {
    return m_body_shape.get_inertia_tensor();
  }

  void update(const spacecraft_state& state) override
  {
    m_state = state;
  }

  // void tick(const state_type& state, const double t)
  // {
  //   m_state.set_state(state(arma::span(0, 5)));
  //   m_attitude = state(arma::span(6, 9));
  //   for (const auto& controller: m_controllers) {
  //     control_input inp = controller->get_control_input(m_state.get_state(), m_attitude, t);
  //   }
  // }

  auto get_additional_state_providers() -> std::vector<std::shared_ptr<additional_state_provider>>
  {
    return { };
  }

  [[nodiscard]] auto get_identifier() const -> std::string
  {
    return m_identifier;
  }

  [[nodiscard]] auto get_controllers() const -> std::vector<std::shared_ptr<controller>>
  {
    return m_controllers;
  }

  auto get_mass() -> double
  {
    return m_state.get_mass();
  }

  /**
   * Returns the string identifier of this spacecraft
   * @return Copy of the spacecraft identifier
   */
  auto get_identifier() -> std::string
  {
    return m_identifier;
  }

  void set_state(const state_type& state)
  {
    m_state.set_pv_coordinates(pv_coordinates(state));
  }

  void set_attitude(const quaternion_type& attitude)
  {
    m_attitude = attitude;
  }

  auto get_maneuver_plan() -> std::shared_ptr<maneuvers::maneuver_plan>
  {
    return m_maneuver_plan;
  }

  void apply_force(const arma::vec3& forces)
  {
    auto state = m_state.get_pv_coordinates();
    // const auto eci_forces = (eci2ric(state) * forces).eval();
    // state(arma::span(3, 5)) += eci_forces;
  }
};
}

#endif //SPACECRAFT_H
