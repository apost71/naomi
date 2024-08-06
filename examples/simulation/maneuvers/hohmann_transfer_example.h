//
// Created by alex_ on 8/5/2024.
//

#ifndef HOHMANN_TRANSFER_EXAMPLE_H
#define HOHMANN_TRANSFER_EXAMPLE_H

#include <armadillo>

#include "bodies/celestial_body.h"
#include "bodies/earth.h"
#include "forces/two_body_force_model.h"
#include "maneuvers/hohmann_transfer.h"
#include "naomi.h"
#include "observers/simulation_observer.h"
#include "orbits/orbits.h"
#include "propagators/numerical_propagator.h"
#include "simulation/simulation.h"
#include "spacecraft/spacecraft.h"
#include "systems/system.h"

using namespace naomi;
using namespace naomi::orbits;
using namespace naomi::bodies;
using namespace naomi::forces;
using namespace naomi::numeric;
using namespace naomi::observers;
using namespace naomi::maneuvers;

inline void simple_hohmann_transfer_maneuver()
{
  constexpr double initial_r = 6378000 + 250000;
  constexpr double target_r = 42164154.0;
  const arma::vec3 r = {initial_r, 0, 0};
  vector_type state_vec = get_circular_orbit(r);
  hohmann_transfer ht(state_vec, target_r);
  auto mp =
      ht.get_maneuver_plan();
  auto sc = std::make_shared<spacecraft>("test", state_vec, 100.0, mp);
  std::shared_ptr<celestial_body> earth_body = std::make_shared<earth>();
  std::shared_ptr<equations_of_motion> two_body_forces = std::make_shared<two_body_force_model_eoms>(earth_body);

  typedef numerical_propagator<rk_dopri5_stepper> propagator;
  typedef physical_system<propagator> system_t;
  const auto system = std::make_shared<system_t>(sc, two_body_forces);
  auto file_observer = std::make_shared<results_csv_writer_observer<system_t>>(
        results_csv_writer_observer<system_t>(10.0, "./test_hohmann_transfer.csv"));
  simulation<system_t> sim(system, {file_observer});
  sim.simulate(60*60*32);
}

inline void simple_hohmann_transfer_maneuver_inclined_delayed_start()
{
  constexpr double target_r = 42164154.0;
  const arma::vec3 r = {3900000.0, 3900000.0, 3900000.0};
  vector_type state_vec = get_circular_orbit(r);
  hohmann_transfer ht(state_vec, target_r);
  auto mp = ht.get_maneuver_plan(keplerian_orbit::get_orbital_period(norm(r)));
  auto sc = std::make_shared<spacecraft>("test", state_vec, 100.0, mp);
  std::shared_ptr<celestial_body> earth_body = std::make_shared<earth>();
  std::shared_ptr<equations_of_motion> two_body_forces = std::make_shared<two_body_force_model_eoms>(earth_body);
  typedef numerical_propagator<rk_dopri5_stepper> propagator;
  typedef physical_system<propagator> system_t;
  auto file_observer = std::make_shared<results_csv_writer_observer<system_t>>(
      observers::results_csv_writer_observer<system_t>(
          10.0,
          "./test_hohmann_transfer_inclined.csv"));

  const auto system = std::make_shared<system_t>(sc, two_body_forces);
  simulation<system_t> sim(system, {file_observer});
  sim.simulate(60.0*60.0*32);
}

inline void bielliptic_hohmann_transfer_maneuver()
{
  constexpr double initial_r = 7000000;
  constexpr double target_r = 105000000;
  const arma::vec3 r = {initial_r, 0, 0};
  vector_type state_vec = get_circular_orbit(r);
  const bielliptic_hohmann_transfer ht(initial_r, target_r, target_r*2);
  auto mp = ht.get_maneuver_plan(0);
  auto sc = std::make_shared<spacecraft>("test", state_vec, 100.0, mp);
  std::shared_ptr<celestial_body> earth_body = std::make_shared<earth>();
  std::shared_ptr<equations_of_motion> two_body_forces = std::make_shared<two_body_force_model_eoms>(earth_body);
  typedef numerical_propagator<rk_dopri5_stepper> propagator;
  typedef physical_system<propagator> system_t;
  const auto system = std::make_shared<system_t>(sc, two_body_forces);
  auto file_observer = std::make_shared<results_csv_writer_observer<system_t>>(
        observers::results_csv_writer_observer<system_t>(10.0, "./test_hohmann_transfer_bielliptic.csv"));
  simulation<system_t> sim(system, {file_observer});
  sim.simulate(60.0*60.0*240);
}

#endif //HOHMANN_TRANSFER_EXAMPLE_H
