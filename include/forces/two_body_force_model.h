//
// Created by alex on 7/4/2024.
//

#ifndef TWO_BODY_FORCE_MODEL_H
#define TWO_BODY_FORCE_MODEL_H
#include "bodies/celestial_body.h"
#include "force_model.h"

namespace naomi::forces
{
using namespace bodies;
class two_body_force_model : public force_model
{
  std::shared_ptr<celestial_body> m_central_body;


public:
  explicit two_body_force_model(
      const std::shared_ptr<celestial_body>& central_body)
      : m_central_body(central_body)
  {
  }

  ~two_body_force_model() override = default;

  void operator()(const vector_type& x,
                  vector_type& dxdt,
                  double t) const override
  {
    arma::vec3 pos = x.subvec(0, 2);
    const arma::vec3 vel = x.subvec(3, 5);
    const arma::vec3 rddot = -1.0 * m_central_body->get_potential_partial(pos);
    dxdt(arma::span(0, 2)) = vel;
    dxdt(arma::span(3, 5)) = rddot;
    dxdt(arma::span(6, 8)) = {0, 0, 0};
  }
};

class two_body_force_model_eoms : public equations_of_motion
{
  std::shared_ptr<celestial_body> m_central_body;


public:
  explicit two_body_force_model_eoms(
      const std::shared_ptr<celestial_body>& central_body)
      : m_central_body(central_body)
  {
  }

  ~two_body_force_model_eoms() override = default;

  [[nodiscard]] vector_type get_derivative(const vector_type& state, double t) const override
  {
    arma::vec3 pos = state.subvec(0, 2);
    const arma::vec3 vel = state.subvec(3, 5);
    const arma::vec3 rddot = -1.0 * m_central_body->get_potential_partial(pos);
    auto dxdt = arma::vec(9);
    dxdt(arma::span(0, 2)) = vel;
    dxdt(arma::span(3, 5)) = rddot;
    dxdt(arma::span(6, 8)) = {0, 0, 0};
    return dxdt;
  }
};
}

#endif //TWO_BODY_FORCE_MODEL_H
