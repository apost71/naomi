//
// Created by alex_ on 6/16/2024.
//

#ifndef SYSTEM_H
#define SYSTEM_H

#include <spacecraft/spacecraft.h>

#include <utility>

#include "bodies/celestial_body.h"
#include "bodies/earth.h"

class physical_system
{
  std::vector<spacecraft> m_spacecrafts;
  std::shared_ptr<celestial_body> m_central_body;
  // A system is made up of a central body and some number of additional perturbing bodies
  // It can have one or many spacecraft
  // fidelity defined at the system level
public:
  explicit physical_system(const spacecraft& spacecraft)
  {
    m_spacecrafts.push_back(spacecraft);
    m_central_body = std::make_shared<celestial_body>(earth());
  }

  physical_system(std::initializer_list<spacecraft>& spacecrafts)
  {
    this(spacecrafts, std::make_shared<celestial_body>(earth()));
  }
  physical_system(const std::initializer_list<spacecraft>& spacecrafts, std::shared_ptr<celestial_body> central_body): m_spacecrafts(spacecrafts), m_central_body(std::move(central_body)){}
};

#endif //SYSTEM_H
