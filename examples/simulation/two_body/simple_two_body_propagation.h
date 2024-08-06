//
// Created by alex on 7/6/2024.
//

#ifndef SIMPLE_TWO_BODY_PROPAGATION_H
#define SIMPLE_TWO_BODY_PROPAGATION_H

#include <armadillo>

#include "attitude/torque_free_provider.h"
#include "bodies/celestial_body.h"
#include "bodies/earth.h"
#include "forces/two_body_force_model.h"
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

template <>
class fmt::formatter<pv_coordinates> {
public:
  constexpr auto parse (format_parse_context& ctx) { return ctx.begin(); }
  template <typename Context>
  constexpr auto format (const pv_coordinates& pv, Context& ctx) const {
    const auto pos = pv.get_position();
    const auto vel = pv.get_velocity();
    return format_to(ctx.out(), "\n\tposition: [{}, {}, {}]\n\tvelocity: [{}, {}, {}]", pos.at(0), pos.at(1), pos.at(2), vel.at(0), vel.at(1), vel.at(2));
  }
};

inline void simple_simulation() {
  const arma::vec r {3900000.0, 3900000.0, 3900000.0};
  vector_type state_vec = get_circular_orbit(r);
  std::shared_ptr<celestial_body> earth_body = std::make_shared<earth>();

  auto sc = std::make_shared<spacecraft>("test", state_vec, 100.0);
  std::shared_ptr<equations_of_motion> two_body_forces = std::make_shared<two_body_force_model_eoms>(earth_body);

  typedef numerical_propagator<rk_dopri5_stepper> propagator;
  typedef physical_system<propagator> system_t;

  const auto system = std::make_shared<system_t>(sc, two_body_forces);
  fmt::print("Start Spacecraft state: {}\n", sc->get_pv_coordinates());
  simulation sim(system);
  sim.simulate(60*60*24); // simulate 24 hours
  fmt::print("End Spacecraft state: {}", sc->get_pv_coordinates());
}

inline void simple_simulation_with_file_observer() {
  const arma::vec r {3900000.0, 3900000.0, 3900000.0};
  vector_type state_vec = get_circular_orbit(r);
  std::shared_ptr<celestial_body> earth_body = std::make_shared<earth>();

  auto sc = std::make_shared<spacecraft>("test", state_vec, 100.0);
  std::shared_ptr<equations_of_motion> two_body_forces = std::make_shared<two_body_force_model_eoms>(earth_body);

  typedef numerical_propagator<rk_dopri5_stepper> propagator;
  typedef physical_system<propagator> system_t;

  const auto system = std::make_shared<system_t>(sc, two_body_forces);
  auto file_observer = std::make_shared<results_csv_writer_observer<system_t>>(results_csv_writer_observer<system_t>(10.0, "/home/alexpost/code/naomi/tmp/test_simulation.csv"));

  fmt::print("Start Spacecraft state: {}\n", sc->get_pv_coordinates());
  simulation<system_t> sim(system, {file_observer});
  sim.simulate(60*60*24); // simulate 24 hours
  fmt::print("End Spacecraft state: {}", sc->get_pv_coordinates());
}

inline void simple_simulation_with_attitude_provider()
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

#endif //SIMPLE_TWO_BODY_PROPAGATION_H
