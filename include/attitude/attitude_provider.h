//
// Created by alex on 7/10/2024.
//

#ifndef ATTITUDE_LAW_H
#define ATTITUDE_LAW_H


#include "naomi.h"


namespace naomi::attitude
{

class additional_state_provider
{

protected:
  std::string _name;
  state_type m_state;
public:
  explicit additional_state_provider(std::string name, state_type state): _name(std::move(name)), m_state(std::move(state)){}
  virtual ~additional_state_provider() = default;
  [[nodiscard]] virtual std::string get_name() const
  {
    return _name;
  }
  [[nodiscard]] virtual state_type get_state() const
  {
    return m_state;
  };
  [[nodiscard]] virtual state_type get_derivative(const state_type& state) const = 0;
  [[nodiscard]] virtual std::size_t get_size()
  {
      return m_state.size();
  };
  // virtual void initialize(const spacecraft_state& state) {}

};


class attitude_provider
{
public:
  virtual ~attitude_provider() = default;
  virtual quaternion_type get_rotation() = 0;
  virtual state_type get_angular_momentum() = 0;
  virtual state_type get_angular_velocity() = 0;
  // virtual std::shared_ptr<additional_state_provider> get_additional_state_provider()
  // {
  //   return nullptr;
  // }
  virtual void apply_force(const arma::vec3& forces) = 0;
};
}
#endif //ATTITUDE_LAW_H
