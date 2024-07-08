//
// Created by alex_ on 6/23/2024.
//

#include <gtest/gtest.h>
#include <armadillo>

#include "frames/transforms.h"
#include "orbits/orbits.h"
#include "spacecraft/state_vector.h"

using namespace naomi::orbits;

TEST(TestFrameTransforms, TestECI2RIC)
{
  double initial_r = 6378000 + 250000;
  arma::vec3 r = {initial_r, 0, 0};
  pv_state_type sv = get_circular_orbit(r);
  auto t = eci2ric(sv);
  auto r_eci = sv(arma::span(0, 2));
  auto v_eci = sv(arma::span(3, 5));

  arma::vec3 r_ric = (r_eci.t() * t).t();
  arma::vec3 v_ric =( v_eci.t() * t).t();

  arma::vec3 r_back = (r_ric.t() * t.t()).t();
  arma::vec3 v_back = (v_ric.t() * t.t()).t();

  r_eci.print("r_eci = ");
  v_eci.print("v_eci = ");
  r_ric.print("r_ric = ");
  v_ric.print("v_ric = ");
  r_back.print("r_back = ");
  v_back.print("v_back = ");


}