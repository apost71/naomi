//
// Created by alex_ on 7/31/2024.
//

#ifndef STATE_PROVIDER_H
#define STATE_PROVIDER_H

#include "forces/force_model.h"
#include "pv_coordinates.h"

class state_provider
{
public:
  virtual ~state_provider() = default;
  [[nodiscard]] virtual pv_coordinates get_pv_coordinates() = 0;
  virtual void apply_control(const naomi::vector_type& control) = 0;
};

class integrated_provider
{
public:
  virtual ~integrated_provider() = default;
  virtual std::shared_ptr<naomi::forces::equations_of_motion> get_eoms() = 0;
  [[nodiscard]] virtual std::size_t get_size() = 0;
  [[nodiscard]] virtual naomi::vector_type get_integrated_state() = 0;
  virtual void set_integrated_state(const naomi::vector_type& state) = 0;
};



#endif //STATE_PROVIDER_H
