//
// Created by alex_ on 6/24/2024.
//

#ifndef EVENT_DETECTOR_H
#define EVENT_DETECTOR_H
#include <armadillo>
#include "spacecraft/state_vector.h"

struct apside_detector
{

  bool operator()(const state_vector& sv) const
  {
    return dot(sv.get_position(), sv.get_velocity()) == norm(cross(sv.get_position(), sv.get_velocity()));
  }
};

#endif //EVENT_DETECTOR_H
