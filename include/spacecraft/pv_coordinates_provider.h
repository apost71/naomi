//
// Created by alex_ on 7/31/2024.
//

#ifndef PV_COORDINATES_PROVIDER_H
#define PV_COORDINATES_PROVIDER_H
#include "forces/force_model.h"
#include "frames/transforms.h"
#include "state_provider.h"

using namespace naomi;

class pv_coordinates_provider final
    :
  public state_provider,
  public integrated_provider
{
  pv_coordinates _state;
  std::shared_ptr<naomi::forces::equations_of_motion> _eoms = nullptr;
public:
  explicit pv_coordinates_provider(pv_coordinates initial_state):
    _state(std::move(initial_state)){}
  ~pv_coordinates_provider() override = default;
  pv_coordinates get_pv_coordinates() override { return _state; }

  state_type get_integrated_state() override
  {
    return _state.to_vec();
  }

  std::size_t get_size() override
  {
    return 9;
  }

  void set_integrated_state(const state_type& state) override
  {
    _state = pv_coordinates(state);
  }

  std::shared_ptr<forces::equations_of_motion> get_eoms() override
  {
    return _eoms;
  }

  void apply_control(const state_type& control) override
  {
    const auto curr_state = _state.to_vec();
    const state_type updated_state = curr_state + control;
    _state = pv_coordinates(updated_state);
  }
};

#endif //PV_COORDINATES_PROVIDER_H
