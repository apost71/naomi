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
#include "systems/system.h"

using namespace naomi;
using namespace naomi::numeric;
using namespace naomi::orbits;
using namespace naomi::bodies;
using namespace naomi::forces;


TEST(HigherOrderGravityTest, BasicAssertions)
{
  const Expression x("x");
  const Expression y("y");
  const Expression z("z");
  const std::shared_ptr<celestial_body> earth_body = std::make_shared<earth>();
  const auto potential = earth_body->get_potential();
  const auto diff_x = potential.diff(x);
  const auto diff_y = potential.diff(y);
  const auto diff_z = potential.diff(z);
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

