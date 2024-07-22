//
// Created by alex on 7/10/2024.
//

#ifndef ATTITUDE_LAW_H
#define ATTITUDE_LAW_H
#include <utility>

#include "simulation/simulation.h"
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


class attitude_provider: public additional_state_provider
{
protected:
  arma::mat33 m_inertia_matrix;

public:
  attitude_provider(const state_type& attitude, const arma::mat33& inertia_matrix): additional_state_provider(attitude), m_inertia_matrix(inertia_matrix) {}
  ~attitude_provider() override = default;
  virtual quaternion_type get_rotation() = 0;
  virtual state_type get_angular_momentum() = 0;
  virtual state_type get_angular_velocity() = 0;
};
}
#endif //ATTITUDE_LAW_H
