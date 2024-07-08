//
// Created by alex on 7/4/2024.
//

#ifndef SPACECRAFT_CONTROLLER_H
#define SPACECRAFT_CONTROLLER_H
#include <memory>

#include "propagators/event_detector.h"
#include "spacecraft/spacecraft.h"

class spacecraft_controller: event_detector
{
public:
  virtual void apply_control(std::shared_ptr<spacecraft>& sc) = 0;

};
#endif //SPACECRAFT_CONTROLLER_H
