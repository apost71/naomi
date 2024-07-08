//
// Created by alex_ on 6/19/2024.
//

#ifndef MANEUVER_H
#define MANEUVER_H
#include <armadillo>

#include "propagators/event_detector.h"

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

  auto get_trigger() const -> std::shared_ptr<event_detector>
  {
    return m_trigger;
  }
};
}



#endif //MANEUVER_H
