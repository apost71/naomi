//
// Created by alex_ on 6/19/2024.
//

#ifndef MANEUVER_H
#define MANEUVER_H
#include <armadillo>

#include "frames/transforms.h"
#include "propagators/event_detector.h"
#include "spacecraft/spacecraft_state.h"

namespace naomi
{
class spacecraft_state;
}
namespace naomi::maneuvers
{
using namespace events;

class maneuver
{
  double m_delta_v;
  arma::vec3 m_direction;  // RIC Frame
  std::shared_ptr<event_detector> m_trigger;

public:

  maneuver(const double dv, const arma::vec3& dir, std::shared_ptr<event_detector>& rt): m_delta_v(dv), m_direction(dir), m_trigger(rt){}

  friend auto operator<<(std::ostream &output, const maneuver &maneuver ) -> std::ostream & {
    output << "delta-v: " << maneuver.m_delta_v << "\n";
    output << "direction: [" << maneuver.m_direction[0] << ", " << maneuver.m_direction[1] << ", " << maneuver.m_direction[2] << "]\n";
    return output;
  }

  [[nodiscard]] auto get_delta_v() const -> arma::vec3
  {
    return m_delta_v * m_direction;
  }

  [[nodiscard]] auto get_delta_v_mag() const -> double
  {
    return m_delta_v;
  }

  [[nodiscard]] auto get_direction() const -> arma::vec3
  {
    return m_direction;
  }

  [[nodiscard]] state_type get_control_input(double dt, spacecraft_state& state) const
  {
    const auto dv = get_delta_v();
    auto curr_pva = state.get_state_provider()->get_pv_coordinates().to_vec();
    state_type curr_pv = curr_pva(arma::span(0, 5));
    const auto eci_control = eci2ric(curr_pv(arma::span(0, 5))) * dv;
    state_type control_inp(9);
    control_inp(arma::span(3, 5)) = eci_control;

    return control_inp;
  }

  auto get_trigger() const -> std::shared_ptr<event_detector>
  {
    return m_trigger;
  }
};
}



#endif //MANEUVER_H
