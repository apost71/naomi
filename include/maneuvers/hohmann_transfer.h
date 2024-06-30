//
// Created by alex_ on 6/19/2024.
//

#ifndef HOHMANN_TRANSFER_H
#define HOHMANN_TRANSFER_H
#include <fmt/core.h>

#include "constants.h"
#include "maneuver.h"
#include "spacecraft/state_vector.h"

class hohmann_transfer
{
  state_vector m_initial_state;
  double m_target_radius;
public:
  hohmann_transfer(const state_vector& initial_state, double target_radius ): m_initial_state(initial_state), m_target_radius(target_radius){}

  auto get_maneuver_plan() -> maneuver_plan
  {
    auto r = m_initial_state.get_position();
    auto v = m_initial_state.get_velocity();
    arma::vec3 h_vec = cross(r, v);
    h_vec.print("h_vec = ");
    arma::vec3 ih = normalise(h_vec);
    ih.print("ih = ");
    arma::vec3 rn = normalise(r);
    rn.print("rn = ");
    arma::vec3 vn = normalise(v);
    arma::vec3 e_vec = cross(v, h_vec) / constants::EARTH_MU - rn;
    double e = norm(e_vec);



    if (e > 1e-6) {
      fmt::print("WARNING: Initial orbit is not circular");
    }

    auto initial_radius = norm(r);
    double h_t = sqrt(2*constants::EARTH_MU*((initial_radius * m_target_radius)/(initial_radius + m_target_radius)));
    double v_tp = h_t / initial_radius;
    double v_ta = h_t / m_target_radius;
    double v_i = sqrt(constants::EARTH_MU / initial_radius);
    double v_f = sqrt(constants::EARTH_MU / m_target_radius);
    double sma_t = (initial_radius + m_target_radius) / 2;
    double transit_time = boost::math::double_constants::pi * sqrt(pow(sma_t, 3)/ constants::EARTH_MU);
    double dv1 = fabs(v_i - v_tp);
    double dv2 = fabs(v_f - v_ta);

    if (m_target_radius > initial_radius) {
      maneuver first(dv1, constants::PLUS_J, 0.0);
      maneuver second(dv2,  constants::MINUS_J, transit_time);
      return maneuver_plan({first, second});
    } else {
      maneuver first(dv1, constants::MINUS_J, 0.0);
      maneuver second(dv2,  constants::PLUS_J, transit_time);
      return maneuver_plan({first, second});
    }
  }
};

#endif //HOHMANN_TRANSFER_H
