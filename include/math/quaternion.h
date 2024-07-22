//
// Created by alex on 7/21/2024.
//

#ifndef QUATERNION_H
#define QUATERNION_H

#include <armadillo>

namespace naomi::math::quaternion
{

typedef arma::vec4 quaternion_type;

/**
 * @brief
 *
 * Taken from https://public.websites.umich.edu/~jbreeden/Quaternions.pdf
 * @param q1
 * @param q2
 * @return
 */
inline quaternion_type product(const quaternion_type& q1, const quaternion_type& q2)
{
  const auto p0 = q1[0]*q2[0] - q1[1]*q2[1] - q1[2]*q2[2] - q1[3]*q2[3];
  const auto p1 = q1[0]*q2[1] + q1[1]*q2[0] + q1[2]*q2[3] - q1[3]*q2[2];
  const auto p2 = q1[0]*q2[2] - q1[1]*q2[3] + q1[2]*q2[0] + q1[3]*q2[1];
  const auto p3 = q1[0]*q2[3] + q1[1]*q2[2] - q1[2]*q2[1] + q1[3]*q2[0];

  return {p0, p1, p2, p3};
}

inline quaternion_type conjugate(const quaternion_type& q)
{
  return {q[0], -q[1], -q[2], -q[3]};
}


}
#endif //QUATERNION_H
