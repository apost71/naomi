//
// Created by alex_ on 6/2/2024.
//

#ifndef CELESTIAL_BODY_H
#define CELESTIAL_BODY_H

#include <Eigen/Dense>
#include <symengine/expression.h>
#include <armadillo>

#include <symengine/lambda_double.h>
#include <symengine/llvm_double.h>

namespace naomi::bodies
{
class celestial_body
{
protected:
  double m_mu;
  double m_soi;
  double m_eq_radius;
  std::vector<double> m_higher_order_terms;
  SymEngine::Expression m_potential_exp;
  SymEngine::RCP<const SymEngine::Basic> m_potential_partial_x;
  SymEngine::RCP<const SymEngine::Basic> m_potential_partial_y;
  SymEngine::RCP<const SymEngine::Basic> m_potential_partial_z;
  SymEngine::LLVMDoubleVisitor m_potential_partial_x_visitor;
  SymEngine::LLVMDoubleVisitor m_potential_partial_y_visitor;
  SymEngine::LLVMDoubleVisitor m_potential_partial_z_visitor;

public:
  virtual ~celestial_body() = default;
  explicit celestial_body(const double mu, const double soi, const double eq_radius, const std::initializer_list<double> higher_order_terms = {}): m_mu(mu), m_soi(soi), m_eq_radius(eq_radius), m_higher_order_terms(higher_order_terms)
  {
    SymEngine::RCP<const SymEngine::Basic> x, y, z, c;
    x = SymEngine::symbol("x");
    y = SymEngine::symbol("y");
    z = SymEngine::symbol("z");
    c = SymEngine::symbol("c");

    auto r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
    auto sin_ph = SymEngine::mul(z, pow(r, -1));
    m_potential_exp = -m_mu * pow(r, -1) - c * pow(r, -3) * (1 - 3 * pow(sin_ph, 2));
    m_potential_partial_x = m_potential_exp.diff(x);
    m_potential_partial_y = m_potential_exp.diff(y);
    m_potential_partial_z = m_potential_exp.diff(z);
    m_potential_partial_x_visitor.init({x, y, z, c}, *m_potential_partial_x);
    m_potential_partial_y_visitor.init({x, y, z, c}, *m_potential_partial_y);
    m_potential_partial_z_visitor.init({x, y, z, c}, *m_potential_partial_z);
  }
  virtual SymEngine::Expression get_potential()
  {
    return m_potential_exp;
  }

  virtual double get_potential(arma::vec& pos)
  {
    SymEngine::Expression x("x");
    SymEngine::Expression y("y");
    SymEngine::Expression z("z");
    SymEngine::Expression c("c");

    auto v = m_potential_exp.subs({
        {x, SymEngine::real_double(pos[0])},
        {y, SymEngine::real_double(pos[1])},
        {z, SymEngine::real_double(pos[2])},
        {c, SymEngine::real_double(m_mu * m_higher_order_terms[0] * pow(m_eq_radius, 2) / 2)}
    });
    return static_cast<double>(v);
  }

  virtual arma::vec get_potential_partial(arma::vec& pos)
  {
    double c = m_mu * m_higher_order_terms[0] * pow(m_eq_radius, 2) / 2;
    const double du_x =
        m_potential_partial_x_visitor.call({pos[0], pos[1], pos[2], c});
    const double du_y =
        m_potential_partial_y_visitor.call({pos[0], pos[1], pos[2], c});
    const double du_z = m_potential_partial_z_visitor.call({pos[0], pos[1], pos[2], c});

    arma::vec result = { du_x, du_y, du_z };
    return result;
  };
  virtual Eigen::Vector3d get_potential_partial_derivative(Eigen::Vector3d position) = 0;
  virtual arma::vec get_potential_partial_derivative(arma::vec position) = 0;

  [[nodiscard]] auto get_mu() const -> double
  {
    return m_mu;
  }

  [[nodiscard]] auto get_sphere_of_influence() const -> double
  {
    return m_soi;
  }
};
}
#endif //CELESTIAL_BODY_H
