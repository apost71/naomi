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

class propagation_state
{
  state_type _state;
  std::map<
    std::string,
    std::pair<arma::span, std::shared_ptr<additional_state_provider>>
  > _additional_state_providers;

  std::map<std::string, std::pair<arma::span, std::shared_ptr<additional_state_provider>>> map_providers(
    const std::vector<std::shared_ptr<additional_state_provider>>& additional_providers
  )
  {
    const std::size_t start_idx = 8;
    std::map<
      std::string,
      std::pair<arma::span, std::shared_ptr<additional_state_provider>>
    > provider_map;
    for (const auto& provider : additional_providers) {
      const auto size = provider->get_size();
      const auto end_idx = start_idx + size;
      const auto prov_span = arma::span(start_idx + 1, end_idx);
      provider_map[provider->get_name()] = {prov_span, provider};
    }
    return provider_map;
  }

public:
  propagation_state(const state_type& state,
      const std::vector<std::shared_ptr<additional_state_provider>>&
          additional_state_providers):
    _state(state), _additional_state_providers(map_providers(additional_state_providers)) {}
};
template <typename Stepper>
class numerical_propagator
{
public:
  numerical_propagator(const numerical_propagator& other)
      : m_integrator(other.m_integrator)
      , m_system(other.m_system)
      , m_spacecrafts(other.m_spacecrafts)
      , _system_eoms(other._system_eoms)
      , m_event_detectors(other.m_event_detectors)
      , m_t(other.m_t)
  {
  }
  numerical_propagator(numerical_propagator&& other) noexcept
      : m_integrator(std::move(other.m_integrator))
      , m_system(std::move(other.m_system))
      , m_spacecrafts(std::move(other.m_spacecrafts))
      , _system_eoms(std::move(other._system_eoms))
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
    _system_eoms = other._system_eoms;
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
    _system_eoms = std::move(other._system_eoms);
    m_event_detectors = std::move(other.m_event_detectors);
    m_t = other.m_t;
    return *this;
  }

private:
  integrator<Stepper> m_integrator;
  std::shared_ptr<force_model> m_system;
  std::shared_ptr<equations_of_motion> _system_eoms;
  std::map<std::string, std::shared_ptr<spacecraft>> m_spacecrafts;
  std::vector<std::shared_ptr<event_detector>> m_event_detectors = {};
  double m_t = 0.0;

public:
  ~numerical_propagator() = default;
  numerical_propagator() = default;

  void initialize(const std::shared_ptr<equations_of_motion>& system_eoms, const std::map<std::string, std::shared_ptr<spacecraft>>& spacecrafts)
  {
    _system_eoms = system_eoms;
    m_spacecrafts = spacecrafts;
    for (const auto & [fst, sc] : m_spacecrafts) {
      if (sc->get_maneuver_plan() != nullptr) m_event_detectors.emplace_back(sc->get_maneuver_plan());
    }
  }

  std::vector<std::pair<arma::span, std::shared_ptr<additional_state_provider>>> map_providers(
    const std::vector<std::shared_ptr<additional_state_provider>>& additional_providers
  )
    {
      const std::size_t start_idx = 8;
      std::vector<
        std::pair<arma::span, std::shared_ptr<additional_state_provider>>
      > providers;
      for (const auto& provider : additional_providers) {
        const auto size = provider->get_size();
        const auto end_idx = start_idx + size;
        const auto prov_span = arma::span(start_idx + 1, end_idx);
        providers.emplace_back(prov_span, provider);
      }
      return providers;
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
    auto system_eoms = _system_eoms;
    const auto provider_map = spacecraft->get_state().get_provider_mapping();
    auto initial_state = spacecraft->get_state().get_integrated_state();
    return [system_eoms, provider_map](const auto& x, auto& dxdt, double t)
        {
          for (const auto& [fst, snd] : provider_map) {
            if (auto eoms = snd->get_eoms(); eoms == nullptr) {
              dxdt(fst) = system_eoms->get_derivative(x(fst), t);
            } else {
              dxdt(fst) = eoms->get_derivative(x(fst), t);
            }
          }
        };
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

  void propagate_to(const std::shared_ptr<spacecraft>& spacecraft, double dt)
  {
    auto system = make_system(m_system, spacecraft);
    const double end = dt;
    auto times = get_integration_times(m_t, end);
    for (std::size_t i = 0; i < times.size() - 1; i++) {
      double start_t = times[i];
      double end_t = times[i + 1];
      state_type state = spacecraft->get_state().get_integrated_state();
      state_and_time_type prev_state = {state, start_t};
      start_t = m_integrator.integrate( system, state , start_t , end_t , 0.1 );
      state_and_time_type new_state = {state, start_t};
      auto active_events = check_events(prev_state, new_state, start_t);
      for (std::shared_ptr<event_detector> e: active_events) {
        // TODO: This wont work with multiple events
        auto event = m_integrator.find_event_time(system, times[i], times[i+1], e, prev_state, 0.1);
        start_t = event.first;
        spacecraft->get_state().set_integrated_state(event.second);
        e->handle_event(spacecraft, start_t);
        spacecraft->update(event.first);
        state = spacecraft->get_state().get_integrated_state();
        std::cout << "event occurred at: " << start_t << "state: " << state << "\n";
      }
      m_integrator.integrate( system, state , start_t , end_t , 0.1 );
      spacecraft->get_state().set_integrated_state(state);
      spacecraft->update(end_t);
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
