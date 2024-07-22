//
// Created by alex on 7/10/2024.
//

#include <armadillo>

#include <gtest/gtest.h>

#include "attitude/torque_free.h"
#include "naomi.h"

using namespace naomi;
using namespace naomi::attitude;
TEST(TestTorqueFree, TestInitialize)
{
  torque_free_attitude tf({1, 2, 3, 4, 5, 6}, arma::mat33().eye());
}