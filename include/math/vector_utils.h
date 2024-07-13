//
// Created by alex on 7/10/2024.
//

#ifndef VECTOR_UTILS_H
#define VECTOR_UTILS_H
#include <naomi.h>

namespace naomi::math
{
inline arma::mat44 q_skew(const quaternion_type& q)
{
  return {
      {q[0], -q[1], -q[2], -q[3]},
      {q[1], q[0], -q[3], q[2]},
      {q[2], q[3], q[0], -q[1]},
      {q[3], -q[2], q[1], q[0]}
  };
}

inline arma::mat44 skew(const arma::vec3& x)
{
  return {
      {0, -x[2], x[1]},
      {x[2], 0, -x[0]},
      {-x[1], x[0], 0}
  };
}
}
#endif //VECTOR_UTILS_H
