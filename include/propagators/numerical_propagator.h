//
// Created by alex_ on 5/24/2024.
//

#ifndef NUMERICAL_PROPAGATOR_H
#define NUMERICAL_PROPAGATOR_H

#include "integrators/integrator.h"
#include "spacecraft/spacecraft.h"
#include "forces/force_model.h"

namespace naomi::numeric
{
using namespace events;
using namespace maneuvers;
template <typename Stepper>
class numerical_propagator
{
public:
  numerical_propagator(const numerical_propagator& other)
      : m_integrator(other.m_integrator)
      , m_system(other.m_system)
      , m_spacecrafts(other.m_spacecrafts)
      , m_event_detectors(other.m_event_detectors)
      , m_t(other.m_t)
  {
  }
  numerical_propagator(numerical_propagator&& other) noexcept
      : m_integrator(std::move(other.m_integrator))
      , m_system(std::move(other.m_system))
      , m_spacecrafts(std::move(other.m_spacecrafts))
      , m_event_detectors(std::move(other.m_event_detectors))
      , m_t(other.m_t)
  {
  }
  numerical_propagator& operator=(const numerical_propagator& other)
  {
    if (this == &other)
      return *this;
    m_integrator = other.m_integrator;
    m_system = other.m_system;
    m_spacecrafts = other.m_spacecrafts;
    m_event_detectors = other.m_event_detectors;
    m_t = other.m_t;
    return *this;
  }
  numerical_propagator& operator=(numerical_propagator&& other) noexcept
  {
    if (this == &other)
      return *this;
    m_integrator = std::move(other.m_integrator);
    m_system = std::move(other.m_system);
    m_spacecrafts = std::move(other.m_spacecrafts);
    m_event_detectors = std::move(other.m_event_detectors);
    m_t = other.m_t;
    return *this;
  }

private:
  integrator<Stepper> m_integrator;
  std::shared_ptr<force_model> m_system;
  std::map<std::string, std::shared_ptr<spacecraft>> m_spacecrafts;
  std::vector<std::shared_ptr<event_detector>> m_event_detectors = {};
  double m_t = 0.0;

public:
  ~numerical_propagator() = default;
  numerical_propagator() = default;



  void initialize(const std::shared_ptr<force_model>& system, const std::map<std::string, std::shared_ptr<spacecraft>>& spacecrafts)
  {
    m_system = system;
    m_spacecrafts = spacecrafts;
    for (const auto & [fst, sc] : m_spacecrafts) {
      if (sc->get_maneuver_plan() != nullptr) m_event_detectors.emplace_back(sc->get_maneuver_plan());
    }
  }

  std::vector<std::shared_ptr<event_detector>> check_events(const state_and_time_type& prev, const state_and_time_type& curr, double t)
  {
    std::vector<std::shared_ptr<event_detector>> active_events;
    for(const std::shared_ptr<event_detector>& e: m_event_detectors) {
      if ((*e)(prev, curr)) {
        active_events.push_back(e);
      }
    }
    return active_events;
  }

  auto make_system(const std::shared_ptr<force_model>& force_model,
                   const std::shared_ptr<spacecraft>& spacecraft)
  {
    state_type initial_state = get_initial_state(spacecraft);
    return std::make_pair(
        initial_state,
        [force_model, spacecraft](const auto& x, auto& dxdt, double t)
        {
          (*force_model)(x, dxdt, t);
          std::size_t start_idx = 8;
          auto addl_state_providers =
              spacecraft->get_additional_state_providers();
          for (const auto& addl_states : addl_state_providers) {
            const auto size = addl_states->get_size();
            const auto end_idx = start_idx + size;
            auto state = x(arma::span(start_idx + 1, end_idx));
            const auto d_state = addl_states->get_derivative(x);
            dxdt(arma::span(start_idx + 1, end_idx)) = d_state;
            start_idx = end_idx;
          }
        });
  }

  static std::vector<double> get_integration_times(const double t_start,
                                                   const double t_end)
  {
    std::vector<double> times;
    double t_curr = t_start;
    double t_step = 2;
    while (t_curr + t_step < t_end) {
      times.push_back(t_curr);
      t_curr += t_step;
    }
    times.push_back(t_end);
    return times;
  }

  static auto get_initial_state(const std::shared_ptr<spacecraft>& spacecraft)
  {
    std::size_t size = 6;
    arma::vec states = {spacecraft->get_state()};
    for (const auto& p: spacecraft->get_additional_state_providers()) {
      size += p->get_size();
      states = join_cols(states, p->get_state());
    }

    return states;
  }

  void propagate_to(const std::shared_ptr<spacecraft>& spacecraft, double dt)
  {
    auto system = make_system(m_system, spacecraft);
    double end = dt;
    auto times = get_integration_times(m_t, end);
    for (std::size_t i = 0; i < times.size() - 1; i++) {
      double start_t = times[i];
      double end_t = times[i + 1];
      state_type state = join_cols(spacecraft->get_state(), spacecraft->get_attitude());
      state_and_time_type prev_state = {state, start_t};
      start_t = m_integrator.integrate( system.second, state , start_t , end_t , 0.1 );
      state_and_time_type new_state = {state, start_t};
      auto active_events = check_events(prev_state, new_state, start_t);
      for (std::shared_ptr<event_detector> e: active_events) {
        // TODO: This wont work with multiple events
        auto event = m_integrator.find_event_time(system.second, times[i], times[i+1], e, prev_state, 0.1);
        start_t = event.first;
        spacecraft->set_state(event.second(arma::span(0, 5)));
        spacecraft->set_attitude(event.second(arma::span(6, 9)));
        e->handle_event(spacecraft, start_t);
        state = join_cols(spacecraft->get_state(), spacecraft->get_attitude());
        std::cout << "event occurred at: " << start_t << "state: " << state << "\n";
      }
      m_integrator.integrate( system.second, state , start_t , end_t , 0.1 );
      // spacecraft->update(state);
      spacecraft->set_state(state(arma::span(0, 5)));
      spacecraft->set_attitude(state(arma::span(6, 9)));
    }
  }

  double propagate_to(const double dt)
  {
    const double end = dt;
    for (const auto & [scid, sc]: m_spacecrafts) {
      propagate_to(sc, dt);
    }
    m_t = end;
    return m_t;
  }

  void propagate_by(const std::shared_ptr<spacecraft>& spacecraft, double dt)
  {
    propagate_to(spacecraft, m_t + dt);
  }

  double propagate_by(const double dt)
  {
    return propagate_to(m_t + dt);
  }
};



typedef
  boost::numeric::odeint::runge_kutta_dopri5<
    state_type,
    double,
    state_type,
    double,
    boost::numeric::odeint::vector_space_algebra> rk_dopri5_stepper;
}


#endif //NUMERICAL_PROPAGATOR_H
