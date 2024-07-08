//
// Created by alex_ on 5/28/2024.
//

#ifndef TWO_BODY_H
#define TWO_BODY_H

// #include "spacecraft/state_vector.h"
// #include "bodies/celestial_body.h"
// #include "spacecraft/spacecraft.h"
//
// class two_body_system
// {
//   std::vector<spacecraft> m_spacecraft;
//   // combined_state_vector m_states;
//   double m_mu;
//   std::shared_ptr<celestial_body> m_central_body;
//
// public:
//
//   two_body_system(const spacecraft& spacecraft, const double& mu, std::shared_ptr<celestial_body>& central_body):
//     m_spacecraft({spacecraft}), m_mu(mu), m_central_body(central_body)
//   {
//     // for (auto & i : m_spacecraft) {
//     //   m_states.add_state_vector( i.get_state());
//     // }
//   }
//
//   two_body_system(const std::vector<spacecraft>& spacecraft, const double& mu, std::shared_ptr<celestial_body>& central_body):
//       m_spacecraft(spacecraft), m_mu(mu), m_central_body(central_body)
//   {
//     // for (int i = 0; i < m_spacecraft.size(); ++i) {
//     //   m_states.add_state_vector( m_spacecraft[i].get_state());
//     // }
//   }
//
//   auto get_state() -> state_type&
//   {
//     return m_spacecraft[0].get_state();
//   }
//
//   auto set_state(state_type& state)
//   {
//     m_spacecraft[0].set_state(state);
//   }
//
//
//   auto get_spacecrafts() -> std::vector<spacecraft>
//   {
//     return m_spacecraft;
//   }
//
//
//   void operator()( const state_type& x , state_type& dxdt, const double t) const
//   {
//     arma::vec3 pos = x.subvec(0, 2);
//     arma::vec3 vel = x.subvec(3, 5);
//     arma::vec3 rddot = m_central_body->get_potential_partial(pos);
//     for (std::size_t i = 0; i < pos.size(); ++i) {
//       rddot[i] *= -1.0;
//     }
//     dxdt(arma::span(0, 2)) = vel;
//     dxdt(arma::span(3, 5)) = rddot;
//   }
// };
#endif //TWO_BODY_H
