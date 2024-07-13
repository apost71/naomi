//
// Created by alex on 7/10/2024.
//

#ifndef ATTITUDE_LAW_H
#define ATTITUDE_LAW_H
#include "spacecraft/state_vector.h"

namespace naomi::attitude
{

class additional_state_provider
{
// protected:
//   state_type m_state;
public:
  // explicit additional_state_provider(const state_type& state): m_state(state){}
  virtual ~additional_state_provider() = default;
  [[nodiscard]] virtual state_type get_derivative(const state_type& state) const = 0;
  [[nodiscard]] virtual std::size_t get_size() const = 0;

};

class attitude_law: public additional_state_provider
{
public:
  explicit attitude_law(): additional_state_provider(){}
  ~attitude_law() override = default;
};
}
#endif //ATTITUDE_LAW_H
