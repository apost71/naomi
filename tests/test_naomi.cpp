//
// Created by alex_ on 6/10/2024.
//

#include <algorithm>
#include <armadillo>
#include <array>
#include <cassert>
#include <iostream>
#include <iterator>
#include <utility>

#include <boost/numeric/odeint/iterator/adaptive_iterator.hpp>
#include <boost/numeric/odeint/iterator/adaptive_time_iterator.hpp>
#include <boost/numeric/odeint/stepper/generation.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta_cash_karp54.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta_dopri5.hpp>
#include <boost/range/adaptor/filtered.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range/numeric.hpp>
#include <gtest/gtest.h>

#include "bodies/celestial_body.h"
#include "bodies/earth.h"
#include "forces/two_body_force_model.h"
#include "naomi.h"
#include "orbits/orbits.h"
#include "propagators/numerical_propagator.h"
#include "spacecraft/state_vector.h"
#include "systems/system.h"

using namespace naomi;
using namespace naomi::numeric;
using namespace naomi::orbits;
using namespace naomi::bodies;
using namespace naomi::forces;
TEST(ArmadilloTest, BasicAssertions)
{
    arma::vec test = {1, 2, 3};
    arma::vec test2 = {2, 3, 4};
    std::cout << test << "\n";
    std::cout << cross(test, test2) << "\n";

    EXPECT_EQ(test[0], 1.0);
}

TEST(HigherOrderGravityTest, BasicAssertions)
{
    Expression x("x");
    Expression y("y");
    Expression z("z");
    std::shared_ptr<celestial_body> earth_body = std::make_shared<earth>();
    auto potential = earth_body->get_potential();
    auto diff_x = potential.diff(x);
    auto diff_y = potential.diff(y);
    auto diff_z = potential.diff(z);
    // simplify(diff_x);
    std::cout << potential << "\n";
    std::cout << diff_x << "\n";
    std::cout << diff_y << "\n";
    std::cout << diff_z << "\n";
}

TEST(HigherOrderGravityTest, PartialDerivative)
{
    arma::vec r {3900000.0, 3900000.0, 3900000.0};

    std::shared_ptr<celestial_body> earth_body = std::make_shared<earth>();
    double vn = sqrt(constants::EARTH_MU/norm(r));
    arma::vec k_hat {0.0, 0.0, 1.0};
    arma::vec h_dir = normalise(cross(r, cross(k_hat, r)));
    arma::vec v_dir = normalise(cross(h_dir, r));
    arma::vec v = v_dir*vn;

    auto sym_partial = earth_body->get_potential_partial(r);
    auto num_partial = earth_body->get_potential_partial_derivative(r);

    std::cout << sym_partial << "\n\n";
    std::cout << num_partial;
}

TEST(NumericalPropagatorTest, CircularOrbitPropagation)
{
  arma::vec r {3900000.0, 3900000.0, 3900000.0};
  state_type state_vec = get_circular_orbit(r);
  std::shared_ptr<celestial_body> earth_body = std::make_shared<earth>();
  std::shared_ptr<spacecraft> sc = std::make_shared<spacecraft>("test", state_vec, 100.0);

  std::shared_ptr<equations_of_motion> two_body_forces = std::make_shared<two_body_force_model_eoms>(earth_body);
  typedef numerical_propagator<rk_dopri5_stepper> propagator;
  physical_system<propagator> system(sc, two_body_forces);
  std::fstream fout;

  // opens an existing csv file or creates a new file.
  fout.open("/home/alexpost/code/naomi/test_sym_prop_j2.csv", std::ios::out);

  fout << "x,y,z,vx,vy,vz\n";


  auto curr_state = sc->get_state().get_integrated_state();
  fout << curr_state[0]  << "," << curr_state[1] << "," << curr_state[2] << ",";
  fout << curr_state[3]  << "," << curr_state[4] << "," << curr_state[5] << "\n";

  double t = 0.0;
  std::cout << "Starting propagation..\n";
  while (t < 60.0*60.0*12) {
    system.simulate_by(10.0);
    curr_state = sc->get_state().get_integrated_state();
    std::cout << t << "\n";
    fout << curr_state[0]  << "," << curr_state[1] << "," << curr_state[2] << ",";
    fout << curr_state[3]  << "," << curr_state[4] << "," << curr_state[5] << "\n";
    t += 10.0;
  }
}

TEST(NumericalPropagatorTest, EllipticalOrbitPropagation)
{
  arma::vec3 r = {1000000, 5000000, 7000000};
  arma::vec3 v = {3000, 4000, 5000};
  state_type state = join_cols(r, v);
  // state_type state = get_circular_orbit(r);
  std::shared_ptr<celestial_body> earth_body = std::make_shared<earth>();
  std::vector<spacecraft> spacecrafts;
  std::shared_ptr<spacecraft> sc = std::make_shared<spacecraft>("test", state, 100.0);

  std::shared_ptr<equations_of_motion> two_body_forces = std::make_shared<two_body_force_model_eoms>(earth_body);
  typedef numerical_propagator<rk_dopri5_stepper> propagator;
  physical_system<propagator> system(sc, two_body_forces);
  std::fstream fout;

  // opens an existing csv file or creates a new file.
  fout.open("/home/alexpost/code/naomi/test_sym_prop_j2_elliptical.csv", std::ios::out);

  fout << "x,y,z,vx,vy,vz\n";

  auto curr_state = sc->get_state().get_integrated_state();
  fout << curr_state[0]  << "," << curr_state[1] << "," << curr_state[2] << ",";
  fout << curr_state[3]  << "," << curr_state[4] << "," << curr_state[5] << "\n";

  double t = 0.0;
  std::cout << "Starting propagation..\n";
  while (t < 60.0*60.0*24) {
    system.simulate_by(10.0);
    curr_state = sc->get_state().get_integrated_state();
    // std::cout << t << "\n";
    fout << curr_state[0]  << "," << curr_state[1] << "," << curr_state[2] << ",";
    fout << curr_state[3]  << "," << curr_state[4] << "," << curr_state[5] << "\n";
    t += 10.0;
  }
  fout.close();
}

