//
// Created by alex_ on 5/25/2024.
//

#ifndef SPACECRAFT_H
#define SPACECRAFT_H
#include <utility>

#include "frames/transforms.h"
#include "maneuvers/maneuver_plan.h"

namespace naomi {

class spacecraft_controller;
class spacecraft
{
  std::string m_identifier;
  pv_state_type m_state;  // TODO: probably should be a pointer
  double m_mass;
  std::vector<std::shared_ptr<spacecraft_controller>> m_controllers;
  std::shared_ptr<maneuvers::maneuver_plan> m_maneuver_plan;

public:
  spacecraft(std::string identifier, const pv_state_type& state, const double& mass):
    m_identifier(std::move(identifier)),
    m_state(state),
    m_mass(mass) {}

  spacecraft(std::string identifier, const pv_state_type& state, const double& mass, const std::shared_ptr<maneuvers::maneuver_plan>& mp):
    m_identifier(std::move(identifier)), m_state(state), m_mass(mass), m_maneuver_plan(mp) {}

  auto get_state() -> pv_state_type&
  {
    return m_state;
  }

  [[nodiscard]] auto get_identifier() const -> std::string
  {
    return m_identifier;
  }

  [[nodiscard]] auto get_controllers() const -> std::vector<std::shared_ptr<spacecraft_controller>>
  {
    return m_controllers;
  }

  auto get_mass() -> double
  {
    return m_mass;
  }

  auto get_identifier() -> std::string
  {
    return m_identifier;
  }

  void set_state(const pv_state_type& state)
  {
    m_state = state;
  }

  auto get_maneuver_plan() -> std::shared_ptr<maneuvers::maneuver_plan>
  {
    return m_maneuver_plan;
  }

  void apply_force(const arma::vec3& forces)
  {
    const auto eci_forces = (eci2ric(m_state) * forces).eval();
    m_state(arma::span(3, 5)) += eci_forces;
  }
};
}

#endif //SPACECRAFT_H
