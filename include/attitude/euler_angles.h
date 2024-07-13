//
// Created by alex on 7/12/2024.
//

#ifndef EULER_ANGLES_H
#define EULER_ANGLES_H
#include "rotation.h"

namespace naomi::attitude
{

enum RotationOrder
{
  XYZ = std::array<int, 3>{1, 2, 3},
  XYX = std::array<int, 3>{1, 2, 1},
  XZY = std::array<int, 3>{1, 3, 2},
  XZX = std::array<int, 3>{1, 3, 1},
  YXZ = std::array<int, 3>{2, 1, 3},
  YXY = std::array<int, 3>{2, 1, 2},
  YZX = std::array<int, 3>{2, 3, 1},
  YZY = std::array<int, 3>{2, 3, 2},
  ZXY = std::array<int, 3>{3, 1, 2},
  ZXZ = std::array<int, 3>{3, 1, 3},
  ZYX = std::array<int, 3>{3, 2, 1},
  ZYZ = std::array<int, 3>{3, 2, 3}
};

class euler_angles: public rotation
{
  int m_first;
  int m_second;
  int m_third;

public:
  euler_angles(const int first, const int second, const int third):
    m_first(first), m_second(second), m_third(third){}
  std::shared_ptr<rotation> apply_to(
      const std::shared_ptr<rotation>& r) override;
  arma::vec3 apply_to(const arma::vec3& r) override;
  arma::mat33 get_dcm() override
  {

  }

  arma::mat33 get_axis_rotation(int axis, double th)
  {
    if (axis == 1) {
      return {
        {1, 0, 0},
        {0, cos(th), sin(th)},
        {0, -sin(th), cos(th)}
      };
    }
    if (axis == 2) {
      return {
          {cos(th), 0, -sin(th)},
          {0, 1, 0},
          {sin(th), 0, cos(th)}
      };
    }
    if (axis == 3) {
      return {
          {cos(th), sin(th), 0},
          {-sin(th), cos(th), 0},
          {0, 0, 1}
      };
    }
  }
};
}

#endif //EULER_ANGLES_H
