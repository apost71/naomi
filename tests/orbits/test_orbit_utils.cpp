//
// Created by alex_ on 6/11/2024.
//
#include <gtest/gtest.h>
#include <armadillo>

// #include <symengine/add.h>
#include <symengine/expression.h>
// #include <symengine/functions.h>
// #include <symengine/mul.h>
// #include <symengine/symbol.h>
// #include <symengine/symengine_rcp.h>

#include "orbits/cartesian.h"
#include "orbits/keplerian.h"
#include "orbits/orbits.h"

TEST(TestOrbitUtils, CircularOrbitVelocity)
{
  arma::vec3 r {6878000.0, 0.0, 0.0};
  arma::vec3 r2 {3900000.0, 3900000.0, 3900000.0};

  state_vector s = get_circular_orbit(r);
  std::cout << s << "\n";
  state_vector s2 = get_circular_orbit(r2);
  std::cout << s2 << "\n";
}

TEST(TestKeplerian, KepFromCart)
{
  arma::vec3 r {5492000,3984000,2955.0};
  arma::vec3 v {-3931.0,5498.0,3665.0};

  cartesian_orbit c = cartesian_orbit(r, v);
  keplerian_orbit k = keplerian_orbit::from_cartesian(c);
  EXPECT_NEAR(k.get_a(), 6827216.613733, 1e-6);
  EXPECT_NEAR(k.get_e(), 0.00880335, 1e-6);
  EXPECT_NEAR(k.get_i(), 28.469632, 1e-6);
  EXPECT_NEAR(k.get_raan(), 35.911819, 1e-6);
  EXPECT_NEAR(k.get_aop(), 314.502161, 1e-6);
  EXPECT_NEAR(k.get_anomaly(), 44.833374, 1e-6);
}

TEST(TestKeplerian, KepFromCart2)
{
  arma::vec3 r {5492000,3984000,2955.0};
  arma::vec3 v {-3000.931,5000.498,3000.665};

  cartesian_orbit c = cartesian_orbit(r, v);
  keplerian_orbit k = keplerian_orbit::from_cartesian(c);
  arma::vec6 res = k.to_vec();
  fmt::print("a = {}\ne = {}\ni = {}\nraan = {}\naop = {}\nma = {}",
    k.get_a(), k.get_e(), k.get_i(true), k.get_raan(true),
    k.get_aop(true), k.get_anomaly(true));

}

TEST(TestKeplerian, CartFromKep)
{
  arma::vec3 r {5492000,3984000,2000.955};
  arma::vec3 v {-3000.931,5000.498,3000.665};

  cartesian_orbit c = cartesian_orbit(r, v);
  keplerian_orbit k = keplerian_orbit::from_cartesian(c);
  arma::vec6 sv = k.to_cartesian();
  EXPECT_NEAR(sv[0], r[0], 1e-6);
  EXPECT_NEAR(sv[1], r[1], 1e-6);
  EXPECT_NEAR(sv[2], r[2], 1e-6);
  EXPECT_NEAR(sv[3], v[0], 1e-6);
  EXPECT_NEAR(sv[4], v[1], 1e-6);
  EXPECT_NEAR(sv[5], v[2], 1e-6);
}


TEST(TestKeplerian, NewtonRaphson)
{
  SymEngine::RCP<const SymEngine::Basic> psi, e, m;
  psi = SymEngine::symbol("psi");
  e = SymEngine::symbol("e");
  m = SymEngine::symbol("m");
  SymEngine::Expression f;

  f = SymEngine::sub(sub(psi, mul(e, sin(psi))),  m);
  auto df = f.diff(psi);
  std::cout << f << "\n";
  std::cout << df << "\n";
}

TEST(TestKeplerian, EccentricAnomalyFromMeanAnomaly)
{
  arma::vec3 r {5492,3984,2.955};
  arma::vec3 v {-3.931,5.498,3.665};

  cartesian_orbit c = cartesian_orbit(r, v);
  keplerian_orbit k = keplerian_orbit::from_cartesian(c);

  auto ecc_anomaly = k.get_eccentric_anomaly();
  fmt::print("Eccentric anomaly: {}", ecc_anomaly);
}

TEST(TestKeplerian, CircularOrbitFromKepOE)
{
  double initial_r = 6378000 + 250000;
  keplerian_orbit k(initial_r, 0, 0, 0, 0, 0, AnomalyType::MEAN);
  arma::vec6 c = k.to_cartesian();
  c.print("state = ");

  arma::vec3 r = {initial_r, 0, 0};
  state_vector sv = get_circular_orbit(r);
  sv.get_position().print("pos = ");
  sv.get_velocity().print("vel = ");
}

TEST(TestKeplerian, TestEccentricity)
{
  arma::vec3 r = {1000000, 5000000, 7000000};
  arma::vec3 v = {3000, 4000, 5000};
  double e_m = keplerian_orbit::compute_eccentricity(r, v, constants::EARTH_MU);
  EXPECT_NEAR(e_m, 0.94754095, 1e-6);
}