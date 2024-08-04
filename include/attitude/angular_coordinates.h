//
// Created by alex_ on 7/24/2024.
//

#ifndef ANGULAR_COORDINATES_H
#define ANGULAR_COORDINATES_H
#include "math/quaternion.h"
#include <naomi.h>

using namespace naomi::math::quaternion;
using namespace naomi;

class angular_coordinates
{
  quaternion_type _rotation;

public:
  explicit angular_coordinates(const quaternion_type& rotation): _rotation(rotation){}

  auto get_rotation() -> quaternion_type&
  {
    return _rotation;
  }

  static angular_coordinates identity()
  {
    return angular_coordinates({1, 0, 0, 0});
  }

  auto to_vec() -> state_type&
  {
   return _rotation;
  }
};

#endif //ANGULAR_COORDINATES_H
