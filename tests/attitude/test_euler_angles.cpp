//
// Created by alex_ on 7/13/2024.
//

#include <armadillo>

#include <fmt/core.h>
#include <gtest/gtest.h>

#include "attitude/euler_angles.h"
#include "naomi.h"

using namespace naomi;
using namespace naomi::attitude;
TEST(TestEulerAngles, TestRotation)
{
  auto order = rotation_order::XYZ;
  auto arr = dynamic_cast<std::array<int, 3>>(order);
  fmt::print("array: {}", arr);

}