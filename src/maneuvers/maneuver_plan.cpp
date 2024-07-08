//
// Created by alex on 7/4/2024.
//

#include "maneuvers/maneuver_plan.h"
#include "spacecraft/spacecraft.h"

namespace naomi::maneuvers
{
void maneuver_plan::execute_maneuver(const std::shared_ptr<spacecraft>& sc)
{
  sc->apply_force(m_maneuvers.at(stage++).get_delta_v());
}
}