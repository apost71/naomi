//
// Created by alex on 7/4/2024.
//

#ifndef TWO_BODY_FORCE_MODEL_H
#define TWO_BODY_FORCE_MODEL_H
#include "attitude/attitude_provider.h"
#include "bodies/celestial_body.h"
#include "force_model.h"

namespace naomi::forces
{
using namespace bodies;
using namespace attitude;
class two_body_force_model : public force_model
{
  std::shared_ptr<celestial_body> m_central_body;
  std::vector<std::shared_ptr<additional_state_provider>> m_additional_state_providers;


public:
  explicit two_body_force_model(
      const std::shared_ptr<celestial_body>& central_body)
      : m_central_body(central_body)
  {
  }

  two_body_force_model(
    const std::shared_ptr<celestial_body>& central_body, const std::vector<std::shared_ptr<additional_state_provider>>& additional_state_providers)
    : m_central_body(central_body), m_additional_state_providers(additional_state_providers)
  {
  }

  ~two_body_force_model() override = default;

  void operator()(const state_type& x,
                  state_type& dxdt,
                  double t) const override
  {
    arma::vec3 pos = x.subvec(0, 2);
    arma::vec3 vel = x.subvec(3, 5);
    arma::vec3 rddot = -1.0 * m_central_body->get_potential_partial(pos);
    dxdt(arma::span(0, 2)) = vel;
    dxdt(arma::span(3, 5)) = rddot;
  }
};
}

#endif //TWO_BODY_FORCE_MODEL_H