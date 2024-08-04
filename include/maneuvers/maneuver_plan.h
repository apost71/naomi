//
// Created by alex on 7/4/2024.
//

#ifndef MANEUVER_PLAN_H
#define MANEUVER_PLAN_H
#include <vector>

#include "maneuver.h"

namespace naomi
{
class spacecraft;
namespace maneuvers
{
class no_such_maneuver_exception : public std::exception {
public:
  const char * what() {
    return std::string("No such maneuver exists.").c_str();
  }
};

class maneuver_plan : public event_detector
{
  std::vector<maneuver> m_maneuvers;
  std::vector<maneuver> _active_maneuvers;
  std::size_t stage = 0;
  double m_start_time = 0.0;

public:
  explicit maneuver_plan(const std::vector<maneuver>& maneuvers)
      : event_detector(ALL)
      , m_maneuvers(maneuvers)
  {
  }

  maneuver_plan(const maneuver_plan& other)
      : event_detector(ALL)
      , m_maneuvers(other.m_maneuvers)
      , m_start_time(other.m_start_time)
  {
    stage = other.stage;
  }

  maneuver_plan(maneuver_plan&& other) noexcept
      : event_detector(ALL)
      , m_maneuvers(std::move(other.m_maneuvers))
      , stage(other.stage)
      , m_start_time(other.m_start_time)
  {
  }

  maneuver_plan& operator=(const maneuver_plan& other)
  {
    if (this == &other)
      return *this;
    m_maneuvers = other.m_maneuvers;
    stage = other.stage;
    m_start_time = other.m_start_time;
    return *this;
  }

  maneuver_plan& operator=(maneuver_plan&& other) noexcept
  {
    if (this == &other)
      return *this;
    m_maneuvers = std::move(other.m_maneuvers);
    stage = other.stage;
    m_start_time = other.m_start_time;
    return *this;
  }

  auto get_maneuvers() -> std::vector<maneuver> { return m_maneuvers; }

  [[nodiscard]] auto get_total_delta_v() const -> double
  {
    double total = 0;
    for (const auto& man : m_maneuvers) {
      total += man.get_delta_v_mag();
    }
    return total;
  }

  auto pop_stage() -> maneuver& { return m_maneuvers.at(stage++); }

  void execute_maneuver(const std::shared_ptr<spacecraft>& sc);

  friend auto operator<<(std::ostream& output, const maneuver_plan& mp)
      -> std::ostream&
  {
    for (std::size_t i = 0; i < mp.m_maneuvers.size(); i++) {
      output << "Manuever " << i << ": \n";
      output << mp.m_maneuvers[i] << "\n\n";
    }
    return output;
  }

  [[nodiscard]] double g(const state_and_time_type& sv) const override
  {
    return m_maneuvers.at(stage).get_trigger()->g(sv);
  }

  vector_type get_control_input(double dt, spacecraft_state& state)
  {
    vector_type control_inp(9);
    for (const auto& maneuver: _active_maneuvers) {
      control_inp += maneuver.get_control_input(dt, state);
    }
    _active_maneuvers.clear();
    return control_inp;
  }

  void handle_event(const std::shared_ptr<spacecraft>& sc, double t) override
  {
    _active_maneuvers.push_back(m_maneuvers.at(stage++));
    if (stage >= m_maneuvers.size()) m_is_active = false;
  }
};
}
}

#endif //MANEUVER_PLAN_H
