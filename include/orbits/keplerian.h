//
// Created by alex_ on 6/14/2024.
//

#ifndef KEPLERIAN_H
#define KEPLERIAN_H
#include "cartesian.h"
#include "constants.h"
#include <fmt/core.h>
#include <cmath>

#include <boost/math/constants/constants.hpp>
#include <boost/math/tools/tuple.hpp>
#include <boost/math/tools/roots.hpp>
#include "integrators/integrator.h"

namespace naomi::orbits
{
enum class AnomalyType{ TRUE, MEAN, ECCENTRIC };

template <class T>
struct eccentric_anomaly_functor
{
  eccentric_anomaly_functor(T const& ecc, T const& ma) : m_ecc(ecc), m_ma(ma)
  { // Constructor stores value to be 'cube-rooted'.
  }
  boost::math::tuple<T, T> operator()(T const& psi)
  { // z is estimate so far.
    return boost::math::make_tuple(
      psi - m_ecc*sin(psi) - m_ma,
      1 - m_ecc * cos(psi)
    );
  }
private:
  T m_ecc;
  T m_ma;
};


class keplerian_orbit
{
  double m_sma;
  double m_ecc;
  double m_inc;
  double m_raan;
  double m_aop;
  double m_ma;
  double m_ea;
  double m_anomaly;
  AnomalyType m_anomaly_type;

public:
  keplerian_orbit(
    const double sma,
    const double ecc = 0,
    const double inc = 0,
    const double raan = 0,
    const double aop = 0,
    const double anomaly = 0,
    const AnomalyType anomaly_type = AnomalyType::MEAN
    ): m_sma(sma), m_ecc(ecc), m_inc(inc), m_raan(raan), m_aop(aop),
        m_anomaly(anomaly), m_anomaly_type(anomaly_type) {}

  static keplerian_orbit from_cartesian(cartesian_orbit& cart)
  {
    arma::vec3 r = cart.get_position();
    arma::vec3 v = cart.get_velocity();
    arma::vec3 h_vec = cross(r, v);
    arma::vec3 ih = normalise(h_vec);
    arma::vec3 rn = normalise(r);
    arma::vec3 e_vec = cross(v, h_vec) / constants::EARTH_MU - rn;
    arma::vec3 ie = normalise(e_vec);
    arma::vec3 ip = cross(ih, ie);

    // RAAN
    arma::vec3 ih_cross_k = cross(constants::PLUS_K, ih);
    arma::vec3 no_hat = normalise(ih_cross_k);
    auto cos_big_om = dot(no_hat, constants::PLUS_I);
    auto sin_big_om = dot(no_hat, constants::PLUS_J);
    auto raan = atan2(sin_big_om, cos_big_om);

    // Inclination
    auto cos_i = dot(constants::PLUS_K, ih);
    auto inc = acos(cos_i);


    // Argument of Perigee
    auto cos_om = dot(no_hat, ie);
    arma::vec3 no_cross_ie = cross(no_hat, ie);
    auto sin_om = dot(ih, no_cross_ie);
    auto aop = atan2(sin_om, cos_om);
    if (aop < 0) {
      aop = 2 * boost::math::double_constants::pi + aop;
    }

    // Eccentricity
    auto e = norm(e_vec);

    // Semi-Major Axis
    auto p = dot(h_vec, h_vec) / constants::EARTH_MU;
    auto sma = p / (1-pow(e, 2));

    // Mean Anomaly
    auto th = acos(1/e*(p/norm(cart.get_position()) -1));
    auto tan_psi = sqrt((1 - e ) / (1 + e))  * tan(th/2);
    auto psi = atan(tan_psi) * 2;
    auto ma = psi - e * sin(psi);
    return {sma, e, inc, raan, aop, ma, AnomalyType::MEAN};
  }

  static double v_norm(arma::vec3 v)
  {
    return sqrt(pow(v[0], 2) + pow(v[1], 2) + pow(v[2], 2));
  }

  static double compute_eccentricity(arma::vec3& r, arma::vec3& v, double mu)
  {
    typedef arma::vec3 t;

    auto rn = norm(r, 2);
    auto vn =  norm(v, 2);
    auto v_r = dot(r/rn, v);
    fmt::print("rn = {}\nvn = {}\nv_r = {}\n\n", rn, vn, v_r);

    t h = cross(r, v);
    h.print("h = ");
    auto hn = norm(h);
    fmt::print("hn = {}\n\n", hn);
    t vxh = cross(v, h);
    t vxh_mu = vxh / mu;
    t norm_r = r / rn;
    vxh.print("vxh = ");
    vxh_mu.print("vxh_mu = ");
    norm_r.print("norm_r = ");
    t e = vxh_mu - norm_r;
    e.print("e = ");
    double en = norm(e, 2);
    return en;
  }

  auto to_cartesian() const -> arma::vec6
  {
    auto psi = get_eccentric_anomaly();
    auto th = 2 * atan(sqrt((1 + m_ecc)/(1 - m_ecc)) * tan(psi/2));
    auto p = m_sma * (1 - pow(m_ecc, 2));
    auto r = p / (1 + m_ecc*cos(th));
    auto h = sqrt(p * constants::EARTH_MU);
    arma::vec3 r_peri = {r * cos(th), r*sin(th), 0};
    arma::vec3 v_peri = {-(constants::EARTH_MU/h) * sin(th), (constants::EARTH_MU/h)*(m_ecc + cos(th)), 0};
    arma::mat a = {
      {cos(m_aop), sin(m_aop), 0},
      {-sin(m_aop), cos(m_aop), 0},
      {0, 0, 1}
    };
    arma::mat b = {
      { 1, 0, 0},
        {0, cos(m_inc), sin(m_inc)},
        {0, -sin(m_inc), cos(m_inc)},
    };
    arma::mat c = {
      {cos(m_raan), sin(m_raan), 0},
      {-sin(m_raan), cos(m_raan), 0},
      {0, 0, 1}
    };
    auto o_pg = (a * b * c).eval();
    auto o_gp = o_pg.t().eval();
    arma::vec3 r_eci = o_gp * r_peri;
    arma::vec3 v_eci = o_gp * v_peri;
    return arma::join_cols(r_eci, v_eci);
  }

  [[nodiscard]] auto get_orbital_period() const -> double
  {
    return get_orbital_period(m_sma);
  }

  auto static get_orbital_period(const double sma) -> double
  {
    return 2.0*boost::math::double_constants::pi / sqrt(constants::EARTH_MU / pow(sma, 3));
  }

  [[nodiscard]] auto get_a() const -> double
  {
    return m_sma;
  }

  [[nodiscard]] auto get_i(bool degrees = true) const -> double
  {
    if (degrees) {
      return rad2deg(m_inc);
    }
    return m_inc;
  }

  [[nodiscard]] auto get_e() const -> double
  {
    return m_ecc;
  }

  [[nodiscard]] auto get_raan(bool degrees = true) const -> double
  {
    if (degrees) {
      return rad2deg(m_raan);
    }
    return m_raan;
  }

  [[nodiscard]] auto get_aop(bool degrees = true) const -> double
  {
    if (degrees) {
      return rad2deg(m_aop);
    }
    return m_aop;
  }

  [[nodiscard]] auto get_anomaly(bool degrees = true) const -> double
  {
    if (degrees) {
      return rad2deg(m_anomaly);
    }
    return m_anomaly;
  }

  [[nodiscard]] auto get_eccentric_anomaly() const -> double
  {
    if (m_anomaly_type == AnomalyType::MEAN) {
      int digits = std::numeric_limits<double>::digits; // Maximum possible binary digits accuracy for type T.
      return boost::math::tools::newton_raphson_iterate(eccentric_anomaly_functor<double>(m_ecc, m_anomaly), 0.1, -2*boost::math::double_constants::pi, 2*boost::math::double_constants::pi, digits);
    }
    if (m_anomaly_type == AnomalyType::TRUE) {
      auto tan_psi = sqrt((1 - m_ecc ) / (1 + m_ecc))  * tan(m_anomaly/2);
      auto psi = atan(tan_psi) * 2;
      return psi;
    }
    return 0.0;
  }

  double rad2deg(double rad) const
  {
    auto result = rad * (180.0/boost::math::double_constants::pi);
    return result;
  }

  double deg2rad(double deg) const
  {
    return deg * (boost::math::double_constants::pi/180.0);
  }

  [[nodiscard]] auto to_vec(bool degrees = false) const -> arma::vec6 {
    if (degrees) {
      return {m_sma, m_ecc, rad2deg(m_inc), rad2deg(m_raan), rad2deg(m_aop), rad2deg(m_anomaly)};
    }
    return {m_sma, m_ecc, m_inc, m_raan, m_aop, m_anomaly};
  }

  [[nodiscard]] auto fn(const double& psi) const -> boost::math::tuple<double, double> {
    return boost::math::make_tuple(
      psi - m_ecc*sin(psi) - m_anomaly,
      1 - m_ecc * cos(psi)
    );
  }

  bool operator==(const keplerian_orbit& other) const
  {
    return
      m_sma == other.m_sma
      && m_ecc == other.m_ecc
      && m_inc == other.m_inc
      && m_raan == other.m_raan
      && m_aop == other.m_aop
      && m_anomaly == other.m_anomaly;
  }
};
}

#endif //KEPLERIAN_H
