//
// Created by alex_ on 6/24/2024.
//

#ifndef EVENT_DETECTOR_H
#define EVENT_DETECTOR_H
#include <armadillo>

#include "event_handler.h"
#include "naomi.h"


namespace naomi::events
{
enum EventDetectorTrigger {INCREASING, DECREASING, ALL};

typedef std::pair<state_type, double> state_and_time_type;

class event_detector
{
  EventDetectorTrigger m_trigger = ALL;
  double m_max_check_interval = 10;
  double m_abs_tol = 1e-6;
  double m_rel_tol = 1e-6;
  std::vector<std::shared_ptr<event_handler>> m_handlers;

protected:
  bool m_is_active = true;

public:
  virtual ~event_detector() = default;
  // event_detector() = default;
  explicit event_detector(const EventDetectorTrigger trigger): m_trigger(trigger){}
  event_detector(
    const EventDetectorTrigger trigger,
    const double max_check_interval
    ): m_trigger(trigger), m_max_check_interval(max_check_interval){}
  event_detector(
    const EventDetectorTrigger trigger,
    const double max_check_interval,
    const double abs_tol
    ): m_trigger(trigger), m_max_check_interval(max_check_interval), m_abs_tol(abs_tol){}
  event_detector(
    const EventDetectorTrigger trigger,
    const double max_check_interval,
    const double abs_tol,
    const double rel_tol
    ): m_trigger(trigger), m_max_check_interval(max_check_interval), m_abs_tol(abs_tol), m_rel_tol(rel_tol){}


  [[nodiscard]] virtual double g(const state_and_time_type& sv) const = 0;

  [[nodiscard]] bool operator()(const state_and_time_type& initial, const state_and_time_type& final) const
  {
    if (! m_is_active) return false;

    const double init_val = g(initial);
    const double final_val = g(final);
    const bool is_root = init_val * final_val <= 0;

    if (m_trigger == INCREASING) {
      return is_root && init_val <= 0;
    }
    if (m_trigger == DECREASING) {
      return is_root && final_val <= 0;
    }
    return init_val * final_val <= 0;
  }

  virtual void handle_event(const std::shared_ptr<spacecraft>& sc, double t)
  {
    for (const auto & handler: m_handlers) {
      handler->handle_event(sc, t);
    }
  }
};

struct event
{
  double time_occurred;
  std::shared_ptr<event_detector> detector;

};

class apside_detector final : public event_detector
{
public:
  explicit apside_detector(const EventDetectorTrigger trigger)
      : event_detector(trigger)
  {
  }

  [[nodiscard]] double g(const state_and_time_type& sv) const override
  {
    const arma::vec3 pos = sv.first(arma::span(0, 2));
    const arma::vec3 vel = sv.first(arma::span(3, 5));
    return dot(pos, vel);
  }
};

class time_detector final : public event_detector
{
  double m_time;
public:
  explicit time_detector(double time)
      : event_detector(DECREASING), m_time(time)
  {
  }

  [[nodiscard]] double g(const state_and_time_type& sv) const override
  {
    return m_time - sv.second;
  }

};

struct event_detector_condition
{
  std::vector<apside_detector> m_detectors;
  std::vector<double> m_detector_g;

  explicit event_detector_condition(apside_detector& ad)
  {
    m_detectors.push_back(ad);
  }

};
}

#endif //EVENT_DETECTOR_H
