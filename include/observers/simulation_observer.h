//
// Created by alex on 7/6/2024.
//

#ifndef SIMULATION_OBSERVER_H
#define SIMULATION_OBSERVER_H
#include <memory>
#include <utility>

#include "propagators/event_detector.h"
#include "spacecraft/pv_coordinates.h"

namespace naomi::observers
{
using namespace events;
template<typename system_t>
class simulation_observer
{
  double m_obs_interval;
  double m_prev_obs_time = 0;
  double m_next_obs_time = 0;

public:
  virtual ~simulation_observer() = default;
  explicit simulation_observer(const double obs_interval): m_obs_interval(obs_interval){}
  virtual void initialize(const std::shared_ptr<system_t>& system){}
  void observe_state(const std::shared_ptr<system_t>& system)
  {
    handle_observe_state(system);
    m_prev_obs_time = system->get_current_time();
    m_next_obs_time = m_prev_obs_time + m_obs_interval;
  }
  virtual void handle_observe_state(const std::shared_ptr<system_t>& system) = 0;
  virtual void observe_event(const std::shared_ptr<event_detector>& event){}
  virtual void terminate(const std::shared_ptr<system_t>& system){}
  [[nodiscard]] double get_next_update() const
  {
    return m_next_obs_time;
  }
};


template<class system_t>
class results_csv_writer_observer: public simulation_observer<system_t>
{
  std::string m_filepath;
  std::fstream m_fout;
  std::fstream::openmode m_openmode = std::ios::out;

public:
  results_csv_writer_observer(const double obs_interval, std::string filepath):
    simulation_observer<system_t>(obs_interval), m_filepath(std::move(filepath)){}

  results_csv_writer_observer(const double obs_interval, std::string filepath, std::fstream::openmode openmode):
    simulation_observer<system_t>(obs_interval), m_filepath(std::move(filepath)), m_openmode(openmode){}

  void initialize(const std::shared_ptr<system_t>& system) override
  {
    m_fout.open(m_filepath, m_openmode);
    m_fout << "scid,x,y,z,vx,vy,vz,q0,q1,q2,q3\n";
    this->observe_state(system);
  }

  void handle_observe_state(const std::shared_ptr<system_t>& system) override
  {
    auto spacecrafts = system->get_spacecrafts();
    for (const auto& [id, sc]: spacecrafts) {
      pv_coordinates pv = sc->get_pv_coordinates();
      auto pos = pv.get_position();
      auto vel = pv.get_velocity();
      auto attitude = sc->get_attitude();
      m_fout << id << ",";
      m_fout << pos[0]  << "," << pos[1] << "," << pos[2] << ",";
      m_fout << pos[3]  << "," << pos[4] << "," << pos[5] << ",";
      m_fout << attitude[0] << "," << attitude[1] << "," << attitude[2] << "," << attitude[3] << "\n";
    }
  }

  void terminate(const std::shared_ptr<system_t>& system) override
  {
    m_fout.close();
  }
};
}

#endif //SIMULATION_OBSERVER_H
