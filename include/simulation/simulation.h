//
// Created by alex on 7/6/2024.
//

#ifndef SIMULATION_H
#define SIMULATION_H
#include "observers/simulation_observer.h"

namespace naomi
{
using namespace observers;

template<typename state_t>
class simulation_component
{
public:
  virtual ~simulation_component() = default;
  virtual void initialize(state_t& state){}
  virtual void update(state_t& state){}
  virtual void terminate(state_t& state){}
};

template<typename system_t>
class simulation
{
  std::shared_ptr<system_t> m_system;
  std::vector<std::shared_ptr<simulation_observer<system_t>>> m_observers;
  double m_t = 0;

public:
  explicit simulation(std::shared_ptr<system_t> system): m_system(system){}
  simulation(std::shared_ptr<system_t> system, const std::initializer_list<std::shared_ptr<simulation_observer<system_t>>>& observers):
    m_system(system), m_observers(observers){}

  std::pair<double, std::shared_ptr<simulation_observer<system_t>>> get_next_update()
  {
    if (m_observers.empty()) {
      return {0, nullptr};
    }
    auto obs = m_observers.at(0);
    return {obs->get_next_update(), obs};
  }

  void simulate(const double duration)
  {
    for(const auto& observer: m_observers) observer->initialize(m_system);

    while (m_t < duration)
    {
      auto next_update = get_next_update();
      double interval = 0;
      if (next_update.second == nullptr) {
        interval = m_t + duration;
      } else {
        interval = next_update.first;
      }
      const double curr_time = m_system->simulate_to(interval);
      if (next_update.second != nullptr) {
        next_update.second->observe_state(m_system);
      }
      m_t = curr_time;
    }
    for(const auto& observer: m_observers) observer->terminate(m_system);
  }
};
}
#endif //SIMULATION_H
