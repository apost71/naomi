//
// Created by alex_ on 6/19/2024.
//

#ifndef MANEUVER_PLAN_H
#define MANEUVER_PLAN_H
#include <vector>
#include <armadillo>
#include <optional>

struct maneuver
{
  double delta_v;
  arma::vec3 direction;  // RIC Frame
  double relative_time;

  maneuver(const double dv, const arma::vec3& dir, const double rt): delta_v(dv), direction(dir), relative_time(rt){}

  friend auto operator<<(std::ostream &output, const maneuver &maneuver ) -> std::ostream & {
    output << "delta-v: " << maneuver.delta_v << "\n";
    output << "direction: [" << maneuver.direction[0] << ", " << maneuver.direction[1] << ", " << maneuver.direction[2] << "]\n";
    output << "relative time: " << maneuver.relative_time;
    return output;
  }

};

class no_such_maneuver_exception : public std::exception {
public:
  const char * what() {
    return std::string("No such maneuver exists.").c_str();
  }
};

class maneuver_plan
{
  std::vector<maneuver> m_maneuvers;
  double m_start_time = 0.0;
public:
  explicit maneuver_plan(const std::vector<maneuver>& maneuvers): m_maneuvers(maneuvers){}

  auto get_maneuvers() -> std::vector<maneuver>
  {
    return m_maneuvers;
  }

  auto get_total_delta_v() -> double
  {
    double total = 0;
    for (const auto& man: m_maneuvers) {
      total += man.delta_v;
    }
    return total;
  }

  auto has_next_maneuver(double t) -> bool
  {
    if (t <= m_start_time && ! m_maneuvers.empty()) {
      return true;
    }
    for (maneuver& m: m_maneuvers) {
      if (t <= m_start_time + m.relative_time) {
        return true;
      }
    }
    return false;
  }

  auto get_next_maneuver(double t) -> maneuver
  {
    if (t <= m_start_time) {
      return m_maneuvers[0];
    }
    for (maneuver& m: m_maneuvers) {
      if (t <= m_start_time + m.relative_time) {
        return m;
      }
    }
    throw no_such_maneuver_exception();
  }

  friend auto operator<<(std::ostream &output, const maneuver_plan &mp ) -> std::ostream & {
    for (int i = 0; i < mp.m_maneuvers.size(); i ++) {
      output << "Manuever " << i << ": \n";
      output << mp.m_maneuvers[i] << "\n\n";
    }
    return output;
  }
};

#endif //MANEUVER_PLAN_H
