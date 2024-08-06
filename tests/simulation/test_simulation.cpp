//
// Created by alex on 7/6/2024.
//
#include <gtest/gtest.h>

#include "bodies/celestial_body.h"
#include "orbits/orbits.h"
#include "propagators/numerical_propagator.h"
#include "simulation/simulation.h"
#include "systems/system.h"

using namespace naomi;
using namespace naomi::numeric;
using namespace naomi::observers;
using namespace naomi::orbits;
using namespace naomi::bodies;
using namespace naomi::forces;

typedef std::shared_ptr<attitude_provider> att_provider_ptr;
