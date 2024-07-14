//
// Created by alex on 7/10/2024.
//

#ifndef ATTITUDE_LAW_H
#define ATTITUDE_LAW_H
#include <utility>

#include "spacecraft/state_vector.h"

namespace naomi::attitude
{

class additional_state_provider
{
  state_type m_state;
public:
  explicit additional_state_provider(state_type state): m_state(std::move(state)){}
  virtual ~additional_state_provider() = default;
  [[nodiscard]] virtual state_type get_state() const
  {
    return m_state;
  };
  [[nodiscard]] virtual state_type get_derivative(const state_type& state) const = 0;
  [[nodiscard]] virtual std::size_t get_size()
  {
      return m_state.size();
  };

};

class attitude_law: public additional_state_provider
{
public:
  explicit attitude_law(const state_type& attitude): additional_state_provider(attitude){}
  ~attitude_law() override = default;
};
}
#endif //ATTITUDE_LAW_H
