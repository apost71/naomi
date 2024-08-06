//
// Created by alex_ on 6/23/2024.
//

#include <armadillo>

#include <gtest/gtest.h>

#include "bodies/earth.h"
#include "forces/two_body_force_model.h"
#include "maneuvers/hohmann_transfer.h"
#include "observers/simulation_observer.h"
#include "orbits/keplerian.h"
#include "orbits/orbits.h"
#include "propagators/numerical_propagator.h"
#include "simulation/simulation.h"
#include "systems/system.h"

using namespace naomi;
using namespace naomi::maneuvers;
using namespace naomi::orbits;
using namespace naomi::numeric;
using namespace naomi::bodies;
using namespace naomi::forces;

TEST(TestHohmann, TestDeltaV)
{
  constexpr double initial_r = 6378000 + 250000;
  constexpr double target_r = 42164154.0;
  const arma::vec3 r = {initial_r, 0, 0};
  const vector_type sv = get_circular_orbit(r);
  hohmann_transfer ht(sv, target_r);
  const auto mp = ht.get_maneuver_plan();

  EXPECT_NEAR(mp->get_total_delta_v(), 3912.172254, 1e-6);
  EXPECT_NEAR(mp->get_maneuvers().at(0).get_delta_v_mag(), 2440.123785, 1e-6);
  EXPECT_NEAR(mp->get_maneuvers().at(1).get_delta_v_mag(), 1472.048468, 1e-6);
}

TEST(TestHohmann, TestDeltaV2)
{
  constexpr double initial_r = 385000000;
  constexpr double target_r = 6378000 + 500000;

  const hohmann_transfer ht(initial_r, target_r);

  EXPECT_NEAR(ht.get_total_dv(), 3885.251730, 1e-6);
  for (const auto & dv: ht.get_dvs()) {
    EXPECT_LT(dv, 0);
  }

  const bielliptic_hohmann_transfer bht(initial_r, target_r, initial_r*2);
  const auto dvs = bht.get_dvs();

  EXPECT_NEAR(bht.get_total_dv(), 3755.0, 3754.642191);
  EXPECT_GT(dvs.at(0), 0);
  EXPECT_LT(dvs.at(1), 0);
  EXPECT_LT(dvs.at(2), 0);
}