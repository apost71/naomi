//
// Created by alex_ on 6/23/2024.
//

#include <armadillo>

#include <gtest/gtest.h>

#include "bodies/earth.h"
#include "forces/two_body_force_model.h"
#include "maneuvers/hohmann_transfer.h"
#include "orbits/keplerian.h"
#include "orbits/orbits.h"
#include "propagators/numerical_propagator.h"
#include "systems/system.h"

using namespace naomi;
using namespace naomi::maneuvers;
using namespace naomi::orbits;
using namespace naomi::numeric;
using namespace naomi::bodies;
using namespace naomi::forces;

TEST(TestHohmann, TestDeltaV)
{
  const double initial_r = 6378000 + 250000;
  const double target_r = 42164154.0;
  const arma::vec3 r = {initial_r, 0, 0};
  const pv_state_type sv = get_circular_orbit(r);
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

TEST(TestHohmann, TestPropagation)
{
  double initial_r = 6378000 + 250000;
  double target_r = 42164154.0;
  arma::vec3 r = {initial_r, 0, 0};
  pv_state_type state_vec = get_circular_orbit(r);
  hohmann_transfer ht(state_vec, target_r);
  auto mp = ht.get_maneuver_plan(keplerian_orbit::get_orbital_period(initial_r));
  std::shared_ptr<spacecraft> sc = std::make_shared<spacecraft>("test", state_vec, 100.0, mp);
  std::shared_ptr<celestial_body> earth_body = std::make_shared<earth>();
  std::shared_ptr<force_model> two_body_forces = std::make_shared<two_body_force_model>(earth_body);
  typedef numerical_propagator<rk_dopri5_stepper> propagator;
  physical_system<propagator> system(sc, two_body_forces);
  std::fstream fout;

  // opens an existing csv file or creates a new file.
  fout.open("/home/alexpost/code/naomi/test_hohmann_transfer.csv", std::ios::out);

  fout << "x,y,z,vx,vy,vz\n";

  auto curr_state = sc->get_state();
  fout << curr_state[0]  << "," << curr_state[1] << "," << curr_state[2] << ",";
  fout << curr_state[3]  << "," << curr_state[4] << "," << curr_state[5] << "\n";

  double t = 0.0;
  std::cout << "Starting propagation..\n";
  while (t < 60.0*60.0*32) {
    system.simulate_by(10.0);
    curr_state = sc->get_state();
    // std::cout << t << "\n";
    fout << curr_state[0]  << "," << curr_state[1] << "," << curr_state[2] << ",";
    fout << curr_state[3]  << "," << curr_state[4] << "," << curr_state[5] << "\n";
    t += 10.0;
  }
  fout.close();
}

TEST(TestHohmann, TestPropagationInclined)
{
  // double initial_r = 6378000 + 250000;
  double target_r = 42164154.0;
  arma::vec3 r = {3900000.0, 3900000.0, 3900000.0};
  pv_state_type state_vec = get_circular_orbit(r);
  hohmann_transfer ht(state_vec, target_r);
  auto mp = ht.get_maneuver_plan(keplerian_orbit::get_orbital_period(norm(r)));
  std::shared_ptr<spacecraft> sc = std::make_shared<spacecraft>("test", state_vec, 100.0, mp);
  std::shared_ptr<celestial_body> earth_body = std::make_shared<earth>();
  std::shared_ptr<force_model> two_body_forces = std::make_shared<two_body_force_model>(earth_body);
  typedef numerical_propagator<rk_dopri5_stepper> propagator;
  physical_system<propagator> system(sc, two_body_forces);
  std::fstream fout;

  // opens an existing csv file or creates a new file.
  fout.open("/home/alexpost/code/naomi/test_hohmann_transfer.csv", std::ios::out);

  fout << "x,y,z,vx,vy,vz\n";

  auto curr_state = sc->get_state();
  fout << curr_state[0]  << "," << curr_state[1] << "," << curr_state[2] << ",";
  fout << curr_state[3]  << "," << curr_state[4] << "," << curr_state[5] << "\n";

  double t = 0.0;
  std::cout << "Starting propagation..\n";
  while (t < 60.0*60.0*32) {
    system.simulate_by(10.0);
    curr_state = sc->get_state();
    // std::cout << t << "\n";
    fout << curr_state[0]  << "," << curr_state[1] << "," << curr_state[2] << ",";
    fout << curr_state[3]  << "," << curr_state[4] << "," << curr_state[5] << "\n";
    t += 10.0;
  }
  fout.close();
}

TEST(TestBiEllipticHohmann, TestSomethingElse)
{
  double initial_r = 7000000;
  double target_r = 105000000;
  arma::vec3 r = {initial_r, 0, 0};
  pv_state_type state_vec = get_circular_orbit(r);
  bielliptic_hohmann_transfer ht(initial_r, target_r, target_r*2);
  auto mp = ht.get_maneuver_plan(0);
  std::shared_ptr<spacecraft> sc = std::make_shared<spacecraft>("test", state_vec, 100.0, mp);
  std::shared_ptr<celestial_body> earth_body = std::make_shared<earth>();
  std::shared_ptr<force_model> two_body_forces = std::make_shared<two_body_force_model>(earth_body);
  typedef numerical_propagator<rk_dopri5_stepper> propagator;
  physical_system<propagator> system(sc, two_body_forces);
  std::fstream fout;

  // opens an existing csv file or creates a new file.
  fout.open("/home/alexpost/code/naomi/test_hohmann_transfer_bielliptic.csv", std::ios::out);

  fout << "x,y,z,vx,vy,vz\n";

  auto curr_state = sc->get_state();
  fout << curr_state[0]  << "," << curr_state[1] << "," << curr_state[2] << ",";
  fout << curr_state[3]  << "," << curr_state[4] << "," << curr_state[5] << "\n";

  double t = 0.0;
  std::cout << "Starting propagation..\n";
  while (t < 60.0*60.0*190) {
    system.simulate_by(10.0);
    curr_state = sc->get_state();
    // std::cout << t << "\n";
    fout << curr_state[0]  << "," << curr_state[1] << "," << curr_state[2] << ",";
    fout << curr_state[3]  << "," << curr_state[4] << "," << curr_state[5] << "\n";
    t += 10.0;
  }
  fout.close();
}