//
// Created by alex_ on 6/19/2024.
//

#ifndef HOHMANN_TRANSFER_H
#define HOHMANN_TRANSFER_H
#include <fmt/core.h>

#include "constants.h"
#include "maneuver.h"
#include "maneuver_plan.h"
#include "orbits/keplerian.h"
#include "spacecraft/state_vector.h"

namespace naomi::maneuvers
{
using namespace orbits;
class hohmann_transfer
{
  pv_state_type m_initial_state;
  double m_target_radius;
  std::vector<double> m_dvs;
  double m_transit_time;

public:
  hohmann_transfer(const pv_state_type& initial_state, double target_radius ):
    m_initial_state(initial_state), m_target_radius(target_radius)
  {
    auto r = m_initial_state(arma::span(0, 2));
    auto v = m_initial_state(arma::span(3, 5));
    arma::vec3 h_vec = cross(r, v);
    arma::vec3 ih = normalise(h_vec);
    arma::vec3 rn = normalise(r);
    arma::vec3 vn = normalise(v);
    arma::vec3 e_vec = cross(v, h_vec) / constants::EARTH_MU - rn;

    if (const double e = norm(e_vec); e > 1e-6) {
      throw std::runtime_error("WARNING: Initial orbit is not near circular");
    }

    auto initial_radius = norm(r);
    double h_t = sqrt(2*constants::EARTH_MU*((initial_radius * m_target_radius)/(initial_radius + m_target_radius)));
    double v_tp = h_t / initial_radius;
    double v_ta = h_t / m_target_radius;
    double v_i = sqrt(constants::EARTH_MU / initial_radius);
    double v_f = sqrt(constants::EARTH_MU / m_target_radius);
    double sma_t = (initial_radius + m_target_radius) / 2;
    m_transit_time = boost::math::double_constants::pi * sqrt(pow(sma_t, 3)/ constants::EARTH_MU);
    m_dvs = {
      v_tp - v_i,
      v_f - v_ta
    };
  }

  hohmann_transfer(const double initial_radius, const double target_radius):
    hohmann_transfer(keplerian_orbit(initial_radius).to_cartesian(), target_radius)
  {}

  auto get_maneuver_plan(double start_time = 0) -> std::shared_ptr<maneuver_plan>
  {
    double initial_radius = norm(m_initial_state(arma::span(0, 2)));
    std::shared_ptr<event_detector> start_detector = std::make_shared<time_detector>(start_time);
    if (m_target_radius > initial_radius) {
      std::shared_ptr<event_detector> apoapsis_detector = std::make_shared<apside_detector>(apside_detector(DECREASING));
      maneuver first(m_dvs.at(0), constants::PLUS_J, start_detector);
      maneuver second(m_dvs.at(1),  constants::PLUS_J, apoapsis_detector);
      return std::make_shared<maneuver_plan>(maneuver_plan({first, second}));
    }
    std::shared_ptr<event_detector> periapsis_detector = std::make_shared<apside_detector>(apside_detector(INCREASING));
    maneuver first(m_dvs.at(0), constants::PLUS_J, start_detector);
    maneuver second(m_dvs.at(1),  constants::PLUS_J, periapsis_detector);
    return std::make_shared<maneuver_plan>(maneuver_plan({first, second}));
  }

  [[nodiscard]] auto get_total_dv() const -> double
  {
    return std::accumulate(m_dvs.begin(), m_dvs.end(), 0.0,
      [](auto x, auto y) {return x + std::abs(y);});
  }

  [[nodiscard]] auto get_dvs() const -> std::vector<double>
  {
    return m_dvs;
  }

  [[nodiscard]] auto get_transit_time() const -> double
  {
    return m_transit_time;
  }
};

class bielliptic_hohmann_transfer
{
  double m_initial_sma;
  double m_target_sma;
  double m_transfer_magnitude;
  double m_transit_time = 0;
  std::vector<double> m_dvs;

public:
  bielliptic_hohmann_transfer(pv_state_type& initial_state, double target_sma, double transfer_magnitude)
  {
    auto r = initial_state(arma::span(0, 2));
    auto v = initial_state(arma::span(3, 5));
    const arma::vec3 h_vec = cross(r, v);
    const arma::vec3 ih = normalise(h_vec);
    const arma::vec3 rn = normalise(r);
    const arma::vec3 vn = normalise(v);
    const arma::vec3 e_vec = cross(v, h_vec) / constants::EARTH_MU - rn;

    if (const double e = norm(e_vec); e > 1e-6) {
      throw std::runtime_error("Initial orbit is not near circular");
    }
    m_initial_sma = norm(initial_state(arma::span(0, 2)));
    m_target_sma = target_sma;
    m_transfer_magnitude = transfer_magnitude;
    double transfer_sma1 = (m_target_sma * m_transfer_magnitude + m_initial_sma) / 2;
    double v_pt = sqrt(constants::EARTH_MU * ((2/m_initial_sma) - (1/transfer_sma1)));
    double r_a1 = transfer_sma1 * 2 - m_initial_sma;
    double transfer_sma2 = (2 * transfer_sma1 + (m_target_sma - m_initial_sma))/2;
    double v_a1 = sqrt(constants::EARTH_MU *  ((2/r_a1) - (1/transfer_sma1)));
    double v_a2 = sqrt(constants::EARTH_MU * ((2/r_a1) - (1/transfer_sma2)));
    double v_p2 = sqrt(constants::EARTH_MU * (2/m_target_sma - 1/transfer_sma2));
    double v_t = sqrt(constants::EARTH_MU * (2/m_target_sma - 1/m_target_sma));
    m_dvs = {
      v_pt - sqrt(constants::EARTH_MU / m_initial_sma),
      v_a2 - v_a1,
      v_t - v_p2
    };

  }

  bielliptic_hohmann_transfer(const double initial_sma,
                              const double target_sma,
                              const double transfer_radius):
    m_initial_sma(initial_sma), m_target_sma(target_sma), m_transfer_magnitude(transfer_radius)
  {
    double transfer_sma1 = (transfer_radius + m_initial_sma) / 2.0;
    double v_initial = vis_viva(m_initial_sma, m_initial_sma);
    double v_pt = vis_viva(m_initial_sma, transfer_sma1);
    double r_a1 = transfer_sma1 * 2.0 - m_initial_sma;
    double transfer_sma2 = (r_a1 + m_target_sma)/2.0;
    double v_a1 = vis_viva(r_a1, transfer_sma1);
    double v_a2 = vis_viva(r_a1, transfer_sma2);

    double v_p2 = vis_viva(m_target_sma, transfer_sma2);
    double v_t = vis_viva(m_target_sma, m_target_sma);

    m_dvs = {
      v_pt - v_initial,
      v_a2 - v_a1,
      v_t - v_p2
    };
  }

  [[nodiscard]] static double vis_viva(const double radius,
                                const double sma,
                                const double mu = constants::EARTH_MU)
  {
    return sqrt(mu * (2/radius - 1/sma));
  }

  [[nodiscard]] auto get_maneuver_plan(double start_time = 0.0) const -> std::shared_ptr<maneuver_plan>
  {
    std::shared_ptr<event_detector> start_detector = std::make_shared<time_detector>(start_time);
    std::shared_ptr<event_detector> periapsis_detector = std::make_shared<apside_detector>(apside_detector(DECREASING));
    std::shared_ptr<event_detector> apoapsis_detector = std::make_shared<apside_detector>(apside_detector(DECREASING));
    maneuver first(m_dvs.at(0), constants::PLUS_J, start_detector);
    maneuver second(m_dvs.at(1),  constants::PLUS_J, apoapsis_detector);
    maneuver third(m_dvs.at(2), constants::PLUS_J, periapsis_detector);

    return std::make_shared<maneuver_plan>(maneuver_plan({first, second, third}));
  }

  [[nodiscard]] auto get_total_dv() const -> double
  {
    return std::accumulate(m_dvs.begin(), m_dvs.end(), 0.0,
      [](auto x, auto y) {return x + std::abs(y);});
  }

  [[nodiscard]] auto get_dvs() const -> std::vector<double>
  {
    return m_dvs;
  }

  [[nodiscard]] auto get_transit_time() const -> double
  {
    return m_transit_time;
  }
};
}

#endif //HOHMANN_TRANSFER_H
