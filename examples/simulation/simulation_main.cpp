//
// Created by alex on 7/6/2024.
//

#include "maneuvers/hohmann_transfer_example.h"
#include "two_body/simple_two_body_propagation.h"

using namespace naomi;

int main(int, char*[])
{
  // Simple simulation with one object, earth J2 eoms, and file output
  simple_simulation();

  // Same simple simulation but it outputs states to a file
  simple_simulation_with_file_observer();

  // Same simple simulation with torque-free attitude propagation
  simple_simulation_with_attitude_provider();

  //////////////////////////////////////////////////////
  ////                 Maneuvers                    ////
  //////////////////////////////////////////////////////

  // Simple hohmann transfer in equatorial orbit
  simple_hohmann_transfer_maneuver();

  simple_hohmann_transfer_maneuver_inclined_delayed_start();

  bielliptic_hohmann_transfer_maneuver();
}
