//
// Created by alex_ on 5/25/2024.
//

#ifndef SPACECRAFT_H
#define SPACECRAFT_H
#include <utility>

#include "attitude/attitude_provider.h"
#include "attitude/constant_attitude_provider.h"
#include "body_shape.h"
#include "control/controller.h"
#include "frames/frame.h"
#include "maneuvers/maneuver_plan.h"
#include "pv_coordinates.h"
#include "pv_coordinates_provider.h"

namespace naomi {

using namespace naomi::control;
using namespace naomi::attitude;
using namespace naomi::geometry;


class spacecraft
{
  std::string m_identifier;
  pv_coordinates m_pv_coordinates;
  quaternion_type m_attitude = {1, 0, 0, 0};
  std::shared_ptr<attitude_provider> _attitude_provider;
  std::shared_ptr<maneuvers::maneuver_plan> m_maneuver_plan;
  body_shape m_body_shape = body_shape::make_rectangle(1, 1, 1, 100);
  spacecraft_state _state;

public:
  spacecraft(std::string identifier, const vector_type& state, const double& mass):
    m_identifier(std::move(identifier)),
    m_pv_coordinates(state),
    _attitude_provider(std::make_shared<constant_attitude_provider>()),
    _state(
      std::make_shared<pv_coordinates_provider>(pv_coordinates(state)),
      std::make_shared<constant_attitude_provider>(),
      mass
    ){}

  spacecraft(std::string identifier, const vector_type& state, const double& mass, const std::shared_ptr<attitude_provider>& attitude_provider):
    m_identifier(std::move(identifier)),
    m_pv_coordinates(state),
    _state(
      std::make_shared<pv_coordinates_provider>(pv_coordinates(state)),
      attitude_provider,
      mass ),
    _attitude_provider(attitude_provider){}

  spacecraft(std::string identifier, const vector_type& state, const double& mass, const std::shared_ptr<maneuvers::maneuver_plan>& mp):
    m_identifier(std::move(identifier)),
    m_pv_coordinates(state),
    _state(
      std::make_shared<pv_coordinates_provider>(pv_coordinates(state)),
      std::make_shared<constant_attitude_provider>(),
      mass),
    _attitude_provider(std::make_shared<constant_attitude_provider>()),
    m_maneuver_plan(mp) {}

  auto get_state() -> spacecraft_state
  {
    return _state;
  }

  /**
   * The current positional coordinates of the spacecraft in a given frame.
   *  TODO: Actually implement frame transformation logic.
   *
   * @param frame The desired frame of the pv coordinates
   * @return The current pv coordinates copied by value
   */
  [[nodiscard]] auto get_pv_coordinates(const std::shared_ptr<frames::frame>& frame = nullptr)
      -> pv_coordinates
  {
    return _state.get_state_provider()->get_pv_coordinates();
  }

  [[nodiscard]] auto get_attitude() const -> quaternion_type
  {
    return _attitude_provider->get_rotation();
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


  [[nodiscard]] auto get_identifier() const -> std::string
  {
    return m_identifier;
  }

  /**
   * Returns the string identifier of this spacecraft
   * @return Copy of the spacecraft identifier
   */
  auto get_identifier() -> std::string
  {
    return m_identifier;
  }

  void set_attitude(const quaternion_type& attitude)
  {
    m_attitude = attitude;
  }

  auto get_maneuver_plan() -> std::shared_ptr<maneuvers::maneuver_plan>
  {
    return m_maneuver_plan;
  }

  void update(const double dt)
  {
    const auto control_inp = m_maneuver_plan->get_control_input(dt, _state);
    _state.get_state_provider()->apply_control(control_inp);
  }
};
}

#endif //SPACECRAFT_H
