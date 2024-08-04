//
// Created by alex on 7/20/2024.
//

#ifndef LQR_CONTROLLER_H
#define LQR_CONTROLLER_H
#include <utility>

#include "controller.h"
#include "frames/transforms.h"

namespace naomi::control
{
class pid_controller: public controller
{
  arma::mat m_kp;
  arma::mat m_ki;
  arma::mat m_kd;
  vector_type m_e;
  double m_t = 0;
  vector_type m_integral;

public:
  pid_controller(arma::mat Kp, arma::mat Ki, arma::mat Kd):
    m_kp(std::move(Kp)), m_ki(std::move(Ki)), m_kd(std::move(Kd)){}

  ~pid_controller() override = default;

  void initialize(const vector_type& state, const vector_type& attitude, const double t) override
  {
    const auto desired = get_desired_state(state, attitude, t);
    m_t = t;
    m_e = desired - state;
  }

  control_input get_control_input(const vector_type& state, const vector_type& attitude, const double t) override
  {
    const auto desired = get_desired_state(state, attitude, t);
    const auto e = desired - state;
    const auto P = m_kp * e;
    m_integral = m_integral * m_ki * e * (t - m_t);
    const auto D = m_kd*(e - m_e)/(t - m_t);
    vector_type control_inp = P + m_integral + D;
    m_e = e;
    m_t = t;
    return {control_inp, attitude};
  }
};

class nadir_pointing_pid_controller final : public pid_controller
{
public:
  nadir_pointing_pid_controller(const arma::mat& Kp,
                                const arma::mat& Ki,
                                const arma::mat& Kd)
      : pid_controller(Kp, Ki, Kd)
  {
  }

  vector_type get_desired_state(const vector_type& state, const vector_type& attitude, double t) override
  {
    const auto trans = eci2ric(state);
    const arma::vec3 desired_vec_ric = {0, 0, 1};
    const arma::vec3 desired_vec_eci = trans.t() * desired_vec_ric;
    return {boost::math::double_constants::pi, desired_vec_eci[0], desired_vec_eci[1], desired_vec_eci[2]};
  }
};
}

#endif //LQR_CONTROLLER_H
