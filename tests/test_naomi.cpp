//
// Created by alex_ on 6/10/2024.
//

#include <gtest/gtest.h>
#include "naomi.h"
#include "bodies/celestial_body.h"
#include <armadillo>
#include <symengine/simplify.h>

#include "orbits/orbits.h"
#include "bodies/earth.h"
#include "propagators/numerical_propagator.h"
#include "spacecraft/state_vector.h"

// Demonstrate some basic assertions.
TEST(HelloTest, BasicAssertions) {
    // Expect two strings not to be equal.
    EXPECT_STRNE("hello", "world");
    // Expect equality.
    EXPECT_EQ(7 * 6, 42);
}

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

TEST(NumericalPropagatorTest, Propagation)
{
  arma::vec r {3900000.0, 3900000.0, 3900000.0};
  state_vector state = get_circular_orbit(r);
  std::shared_ptr<celestial_body> earth_body = std::make_shared<earth>();
  std::vector<spacecraft> spacecrafts;
  spacecraft sc(state, 100.0);
  spacecrafts.push_back(sc);
  // spacecrafts.push_back(sc2);

  two_body_system system(spacecrafts, constants::EARTH_MU, earth_body);
  typedef
    boost::numeric::odeint::runge_kutta_dopri5<
      state_vector,
      double,
      state_vector,
      double,
      boost::numeric::odeint::vector_space_algebra> rk4;
  numerical_propagator< rk4> propagator(system);
  std::fstream fout;

  // opens an existing csv file or creates a new file.
  fout.open("C:\\Users\\alex_\\code\\naomi\\test_sym_prop_j2.csv", std::ios::out);

  fout << "x,y,z,vx,vy,vz\n";

  auto curr_state = system.get_spacecrafts()[0].get_state();
  fout << curr_state.get_position()[0]  << "," << curr_state.get_position()[1] << "," << curr_state.get_position()[2] << ",";
  fout << curr_state.get_velocity()[0]  << "," << curr_state.get_velocity()[1] << "," << curr_state.get_velocity()[2] << "\n";

  double t = 0.0;
  std::cout << "Starting propagation..\n";
  while (t < 60.0*60.0*12) {
    propagator.propagate_by(10.0);
    auto curr_state = system.get_spacecrafts()[0].get_state();
    std::cout << t << "\n";
    fout << curr_state.get_position()[0]  << "," << curr_state.get_position()[1] << "," << curr_state.get_position()[2] << ",";
    fout << curr_state.get_velocity()[0]  << "," << curr_state.get_velocity()[1] << "," << curr_state.get_velocity()[2] << "\n";
    t += 10.0;
  }
}