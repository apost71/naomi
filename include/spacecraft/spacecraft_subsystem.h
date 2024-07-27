//
// Created by alex_ on 7/24/2024.
//

#ifndef SPACECRAFT_SUBSYSTEM_H
#define SPACECRAFT_SUBSYSTEM_H
#include <memory>
#include <string>
#include <utility>

#include "attitude/attitude_provider.h"
#include "spacecraft_state.h"

class spacecraft_subsystem_state
{

};

class spacecraft_subsystem: simulation_component<spacecraft_state>
{
  std::string _name;

public:
  ~spacecraft_subsystem() override = default;
  explicit spacecraft_subsystem(std::string name): _name(std::move(name)){}
  auto get_name() const -> std::string
  {
    return _name;
  }
  virtual std::vector<std::shared_ptr<attitude::additional_state_provider>> get_additional_state_providers() = 0;
};


#endif //SPACECRAFT_SUBSYSTEM_H
