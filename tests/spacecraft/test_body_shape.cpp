//
// Created by alex on 7/15/2024.
//
#include <armadillo>

#include <gtest/gtest.h>

#include "math/volume_integrals.h"
#include "spacecraft/body_shape.h"

using namespace naomi::geometry::test;

TEST(TestBodyShape, TestPolygon)
{
  POLYHEDRON p;
  double density, mass;
  double r[3];            /* center of mass */
  double J[3][3];         /* inertia tensor */

  std::string filename = R"(/home/alexpost/code/naomi/test_cube.txt)";
  const int length = filename.length();
  const char* s = filename.c_str();


  // declaring character array (+1 for null terminator)
  readPolyhedron(s, &p);

  compVolumeIntegrals(&p);


  printf("\nT1 =   %+20.6f\n\n", T0);

  printf("Tx =   %+20.6f\n", T1[X]);
  printf("Ty =   %+20.6f\n", T1[Y]);
  printf("Tz =   %+20.6f\n\n", T1[Z]);

  printf("Txx =  %+20.6f\n", T2[X]);
  printf("Tyy =  %+20.6f\n", T2[Y]);
  printf("Tzz =  %+20.6f\n\n", T2[Z]);

  printf("Txy =  %+20.6f\n", TP[X]);
  printf("Tyz =  %+20.6f\n", TP[Y]);
  printf("Tzx =  %+20.6f\n\n", TP[Z]);

  density = 1.0;  /* assume unit density */

  mass = density * T0;

  /* compute center of mass */
  r[X] = T1[X] / T0;
  r[Y] = T1[Y] / T0;
  r[Z] = T1[Z] / T0;

  /* compute inertia tensor */
  J[X][X] = density * (T2[Y] + T2[Z]);
  J[Y][Y] = density * (T2[Z] + T2[X]);
  J[Z][Z] = density * (T2[X] + T2[Y]);
  J[X][Y] = J[Y][X] = - density * TP[X];
  J[Y][Z] = J[Z][Y] = - density * TP[Y];
  J[Z][X] = J[X][Z] = - density * TP[Z];

  /* translate inertia tensor to center of mass */
  J[X][X] -= mass * (r[Y]*r[Y] + r[Z]*r[Z]);
  J[Y][Y] -= mass * (r[Z]*r[Z] + r[X]*r[X]);
  J[Z][Z] -= mass * (r[X]*r[X] + r[Y]*r[Y]);
  J[X][Y] = J[Y][X] += mass * r[X] * r[Y];
  J[Y][Z] = J[Z][Y] += mass * r[Y] * r[Z];
  J[Z][X] = J[X][Z] += mass * r[Z] * r[X];

  printf("center of mass:  (%+12.6f,%+12.6f,%+12.6f)\n\n", r[X], r[Y], r[Z]);

  printf("inertia tensor with origin at c.o.m. :\n");
  printf("%+15.6f  %+15.6f  %+15.6f\n", J[X][X], J[X][Y], J[X][Z]);
  printf("%+15.6f  %+15.6f  %+15.6f\n", J[Y][X], J[Y][Y], J[Y][Z]);
  printf("%+15.6f  %+15.6f  %+15.6f\n\n", J[Z][X], J[Z][Y], J[Z][Z]);
  /*
  T1 =           +8000.000000

  Tx =              +0.000000
  Ty =              +0.000000
  Tz =              +0.000000

  Txx =        +266666.666667
  Tyy =        +266666.666667
  Tzz =        +266666.666667

  Txy =             +0.000000
  Tyz =             +0.000000
  Tzx =             +0.000000

  center of mass:  (   +0.000000,   +0.000000,   +0.000000)

  inertia tensor with origin at c.o.m. :
  +533333.333333        +0.000000        +0.000000
  +0.000000   +533333.333333        +0.000000
  +0.000000        +0.000000   +533333.333333
  */
}

TEST(TestBodyShape, TestBodyShapeAlgorithm)
{
  auto body = naomi::geometry::body_shape::make_rectangle(20, 20, 20, 100);
  auto Im = body.get_inertia_tensor();
  const arma::mat33 expected = arma::eye(3, 3) * 533333.333333;
  const arma::mat33 res = expected - Im;
  ASSERT_TRUE(res.is_zero(1e-6));
}