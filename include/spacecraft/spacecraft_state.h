//
// Created by alex on 7/21/2024.
//

#ifndef SPACECRAFT_STATE_H
#define SPACECRAFT_STATE_H
#include <utility>

#include "attitude/attitude_provider.h"
#include "simulation/simulation.h"
#include "state_vector.h"

namespace naomi
{

using namespace naomi::attitude;


class state_provider
{

public:
  virtual ~state_provider() = default;
  virtual state_type get_state() = 0;
  virtual void update_state(const state_type& state) = 0;
  virtual std::size_t get_size() = 0;
};


class mass_provider
{

};

class spacecraft_state
{
  state_type m_state;
  std::shared_ptr<attitude_provider> m_attitude;
  double m_mass;

public:
  spacecraft_state(state_type m_state,
                   const std::shared_ptr<attitude_provider>& m_attitude,
                   const double m_mass)
      : m_state(std::move(m_state))
      , m_attitude(m_attitude)
      , m_mass(m_mass)
  {
  }

  auto get_state() -> state_type&
  {
    return m_state;
  }

  void set_state(const state_type& state)
  {
    m_state = state;
  }

  auto get_attitude_provider() -> std::shared_ptr<attitude_provider>&
  {
    return m_attitude;
  }

  [[nodiscard]] auto get_mass() const -> double
  {
    return m_mass;
  }
};

class spacecraft_component: simulation_component<spacecraft_state>{};

}
#endif //SPACECRAFT_STATE_H
