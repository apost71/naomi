//
// Created by alex_ on 8/1/2024.
//

#ifndef CONSTANT_ATTITUDE_PROVIDER_H
#define CONSTANT_ATTITUDE_PROVIDER_H

#include "attitude_provider.h"

namespace naomi::attitude
{
class constant_attitude_provider final : public attitude_provider
{
  quaternion_type _attitude;
public:
  constant_attitude_provider() = default;
  explicit constant_attitude_provider(const quaternion_type& q): _attitude(q){}
  ~constant_attitude_provider() override = default;
  quaternion_type get_rotation() override { return _attitude; }
  state_type get_angular_momentum() override
  {
    return {0, 0, 0};
  };
  state_type get_angular_velocity() override
  {
    return {0, 0, 0};
  };
  void apply_force(const arma::vec3& forces) override {}
};
}
#endif //CONSTANT_ATTITUDE_PROVIDER_H
