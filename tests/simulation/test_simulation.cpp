//
// Created by alex on 7/6/2024.
//
#include <gtest/gtest.h>

#include "attitude/torque_free.h"
#include "attitude/torque_free_provider.h"
#include "bodies/celestial_body.h"
#include "bodies/earth.h"
#include "forces/two_body_force_model.h"
#include "orbits/orbits.h"
#include "propagators/numerical_propagator.h"
#include "simulation/simulation.h"
#include "systems/system.h"

using namespace naomi;
using namespace naomi::numeric;
using namespace naomi::observers;
using namespace naomi::orbits;
using namespace naomi::bodies;
using namespace naomi::forces;

typedef std::shared_ptr<attitude_provider> att_provider_ptr;

TEST(TestSimulation, TestSimulationInitialization)
{
  arma::vec r {3900000.0, 3900000.0, 3900000.0};
  vector_type state_vec = get_circular_orbit(r);
  std::shared_ptr<celestial_body> earth_body = std::make_shared<earth>();

  std::shared_ptr<spacecraft> sc = std::make_shared<spacecraft>("test", state_vec, 100.0);
  std::shared_ptr<equations_of_motion> two_body_forces = std::make_shared<two_body_force_model_eoms>(earth_body);


  typedef numerical_propagator<rk_dopri5_stepper> propagator;
  typedef physical_system<propagator> system_t;


  auto system = std::make_shared<system_t>(sc, two_body_forces);
  auto file_observer = std::make_shared<results_csv_writer_observer<system_t>>(results_csv_writer_observer<system_t>(10.0, "/home/alexpost/code/naomi/tmp/test_simulation.csv"));
  simulation<system_t> sim(system, {file_observer});
  sim.simulate(60*60*24);
}

TEST(TestSimulation, TestAttitude)
{
  arma::vec r {3900000.0, 3900000.0, 3900000.0};
  vector_type state_vec = get_circular_orbit(r);
  std::shared_ptr<celestial_body> earth_body = std::make_shared<earth>();
  body_shape geom = body_shape::make_rectangle(1, 1, 1, 100);
  vector_type initial_attitude = {1, 0, 0, 0};
  std::shared_ptr<attitude_provider> torque_free_attitude =
    std::make_shared<torque_free_attitude_provider>(
      torque_free_attitude_provider(geom.get_inertia_tensor(), {1, 0, 0, 0}, pv_coordinates(state_vec))
    );
  auto sc = std::make_shared<spacecraft>("test", state_vec, 100.0, torque_free_attitude);

  std::shared_ptr<equations_of_motion> two_body_forces = std::make_shared<two_body_force_model_eoms>(earth_body);
  typedef numerical_propagator<rk_dopri5_stepper> propagator;
  typedef physical_system<propagator> system_t;
  auto system = std::make_shared<system_t>(system_t(sc, two_body_forces));
  auto file_observer = std::make_shared<results_csv_writer_observer<system_t>>(results_csv_writer_observer<system_t>(10.0, "/home/alexpost/code/naomi/tmp/test_simulation_attitude.csv"));
  simulation<system_t> sim(system, {file_observer});
  sim.simulate(60*60*24);
}


// TEST(TestSimulation, TestStructure)
// {
//   arma::vec r {3900000.0, 3900000.0, 3900000.0};
//   vector_type state_vec = get_circular_orbit(r);
//   std::shared_ptr<celestial_body> earth_body = std::make_shared<earth>();
//   body_shape geom = body_shape::make_rectangle(1, 1, 1, 100);
//   vector_type initial_attitude = {1, 0, 0, 0};
//
//   std::shared_ptr<attitude_provider> torque_free_attitude =
//     std::make_shared<attitude::torque_free_attitude>(
//       geom.get_inertia_tensor()
//       );
//
//   std::shared_ptr<spacecraft> sc = std::make_shared<spacecraft>(
//     "test",
//     state_vec,
//     100.0,
//     torque_free_attitude
//     );
//
//   std::shared_ptr<force_model> two_body_forces = std::make_shared<two_body_force_model>(earth_body);
//   typedef numerical_propagator<rk_dopri5_stepper> propagator;
//   typedef physical_system<propagator> system_t;
//   auto system = std::make_shared<system_t>(system_t(sc, two_body_forces));
//   auto file_observer = std::make_shared<results_csv_writer_observer<system_t>>(results_csv_writer_observer<system_t>(10.0, "/home/alexpost/code/naomi/test_simulation_attitude.csv"));
//   simulation<system_t> sim(system, {file_observer});
//   sim.simulate(60*60*24);
// }