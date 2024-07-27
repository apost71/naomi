//
// Created by alex on 7/21/2024.
//

#ifndef SPACECRAFT_STATE_H
#define SPACECRAFT_STATE_H
#include <utility>

#include "attitude/angular_coordinates.h"
#include "pv_coordinates.h"
#include "simulation/simulation.h"
#include "state_vector.h"

namespace naomi
{

class state_provider
{

public:
  virtual ~state_provider() = default;
  virtual state_type get_state() = 0;
  virtual void update_state(const state_type& state) = 0;
  virtual std::size_t get_size() = 0;
};


class mass_provider
{

};

class spacecraft_state
{
  pv_coordinates _pv_coordinates;
  angular_coordinates _angular_coordinates = angular_coordinates::identity();
  double _mass;

public:
  spacecraft_state(const state_type& state,
                   angular_coordinates angular_coordinates,
                   const double mass)
      : _pv_coordinates(state)
      , _angular_coordinates(std::move(angular_coordinates))
      , _mass(mass)
  {
  }

  spacecraft_state(const state_type& state,
                 const double mass)
    : _pv_coordinates(state)
    , _mass(mass)
  {
  }

  void update(const state_type& state)
  {
    const pv_coordinates new_pv(state);
    _pv_coordinates = new_pv;
  }

  [[nodiscard]] auto get_mass() const -> double
  {
    return _mass;
  }

  [[nodiscard]] auto get_pv_coordinates() -> pv_coordinates
  {
    return _pv_coordinates;
  }

  void set_pv_coordinates(const pv_coordinates& pv)
  {
    _pv_coordinates = pv;
  }
};

class spacecraft_component: simulation_component<spacecraft_state>{};

}
#endif //SPACECRAFT_STATE_H
