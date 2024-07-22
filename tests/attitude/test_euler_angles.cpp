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
  euler_angles ea(0.5, 0.5, 0.5, rotation_order::XYZ);
  auto dcm = ea.get_dcm();
  dcm.print("DCM: ");

}