//
// Created by alex on 7/12/2024.
//

#ifndef EULER_ANGLES_H
#define EULER_ANGLES_H
#include "rotation.h"

namespace naomi::attitude
{

enum rotation_order
{
  XYZ,
  XYX,
  XZY,
  XZX,
  YXZ,
  YXY,
  YZX,
  YZY,
  ZXY,
  ZXZ,
  ZYX,
  ZYZ
};

class euler_angles: public rotation
{
  double m_alpha;
  double m_beta;
  double m_gamma;
  rotation_order m_rotation_order;

  std::map<rotation_order, std::array<int, 3>> m_rotation_order_map = {
    {XYZ, {1, 2, 3}},
    {XYX, {1, 2, 1}},
    {XZY, {1, 3, 2}},
    {XZX, {1, 3, 1}},
    {YXZ, {2, 1, 3}},
    {YXY, {2, 1, 2}},
    {YZX, {2, 3, 1}},
    {YZY, {2, 3, 2}},
    {ZXY, {3, 1, 2}},
    {ZXZ, {3, 1, 3}},
    {ZYX, {3, 2, 1}},
    {ZYZ, {3, 2, 3}}
  };

public:
  /**
   * Radians
   * @brief
   * @param alpha
   * @param beta
   * @param gamma
   * @param order
   */
  euler_angles(const double alpha, const double beta, const double gamma, const rotation_order order):
    m_alpha(alpha), m_beta(beta), m_gamma(gamma), m_rotation_order(order){}
  ~euler_angles() override = default;
  std::shared_ptr<rotation> apply_to(
      const std::shared_ptr<rotation>& r) override
  {
    return r;
  }

  arma::vec3 apply_to(const arma::vec3& r) override
  {
    return r;
  }

  arma::mat33 get_dcm() override
  {
    std::array<int, 3> rotation_order_arr = m_rotation_order_map[m_rotation_order];
    const auto rot1 = get_axis_rotation(rotation_order_arr[0], m_alpha);
    const auto rot2 = get_axis_rotation(rotation_order_arr[1], m_beta);
    const auto rot3 = get_axis_rotation(rotation_order_arr[2], m_gamma);
    return rot1*rot2*rot3;

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
    throw std::runtime_error(fmt::format("Invalid axis: {}", axis));
  }
};
}

#endif //EULER_ANGLES_H
