//
// Created by alex_ on 6/23/2024.
//

#include <gtest/gtest.h>
#include <armadillo>

#include "maneuvers/hohmann_transfer.h"
#include "orbits/keplerian.h"
#include "orbits/orbits.h"


TEST(TestHohmann, TestDeltaV)
{
  double initial_r = 6378000 + 250000;
  double target_r = 42164154.0;
  arma::vec3 r = {initial_r, 0, 0};
  state_vector sv = get_circular_orbit(r);
  hohmann_transfer ht(sv, target_r);
  auto mp = ht.get_maneuver_plan();
  std::cout << mp << "\n";

  EXPECT_NEAR(mp.get_total_delta_v(), 3912.17, 1e-2);
}

