//
// Created by alex_ on 6/16/2024.
//

#ifndef SYSTEM_H
#define SYSTEM_H

#include <utility>
#include <spacecraft/spacecraft.h>
#include "forces/force_model.h"
#include <fmt/core.h>

namespace naomi
{
using namespace forces;

/**  A `physical_system` contains the force model definitions of the system to
 * be simulated, the spacecrafts that exist within the system, and the
 * propagator that will be used to integrate the dynamics.
 *
 * @brief A definition of the physical system to be simulated.
 * @tparam Propagator The propagator to use for simulation, right now the only
 * type allowed is `NumericalPropagator<T>` where T defines the stepper (see
 * boost docs).  Future work hopefully enables more customizability.
 */
template <typename Propagator>
class physical_system
{
  std::map<std::string, std::shared_ptr<spacecraft>> m_spacecrafts;
  std::shared_ptr<equations_of_motion> _system_eoms;
  Propagator m_propagator;

  double m_t = 0;
  // A system is made up of a central body and some number of additional perturbing bodies
  // It can have one or many spacecraft
  // fidelity defined at the system level
public:
  /**
   * @brief
   * @param spacecraft
   * @param system_eoms
   */
  physical_system(const spacecraft& spacecraft, const std::shared_ptr<equations_of_motion>& system_eoms): physical_system({spacecraft}, system_eoms){
  }

  /**
   * @brief
   * @param spacecraft
   * @param system_eoms
   */
  physical_system(const std::shared_ptr<spacecraft>& spacecraft, const std::shared_ptr<equations_of_motion>& system_eoms): physical_system({spacecraft}, system_eoms){
  }

  /**
   * @brief
   * @param spacecrafts
   * @param system_eoms
   */
  physical_system(const std::initializer_list<spacecraft>& spacecrafts, const std::shared_ptr<equations_of_motion>& system_eoms):
    _system_eoms(system_eoms)
  {
    for (const spacecraft& s: spacecrafts) {
      m_spacecrafts[s.get_identifier()] = std::make_shared<spacecraft>(s);
    }
    m_propagator.initialize(_system_eoms, m_spacecrafts);
  }

  /**
   * @brief
   * @param spacecrafts
   * @param system_eoms
   */
  physical_system(const std::initializer_list<std::shared_ptr<spacecraft>>& spacecrafts, const std::shared_ptr<equations_of_motion>& system_eoms):
  _system_eoms(system_eoms)
  {
    for (const std::shared_ptr<spacecraft>& s: spacecrafts) {
      m_spacecrafts[s->get_identifier()] = s;
    }
    m_propagator.initialize(_system_eoms, m_spacecrafts);
  }

  /**
   * @brief 
   * @param scid
   * @return 
   */
  auto get_spacecraft(const std::string& scid) -> std::shared_ptr<spacecraft>
  {
    return m_spacecrafts[scid];
  }

  /**
   * @brief
   * @return
   */
  auto get_spacecrafts() -> std::map<std::string, std::shared_ptr<spacecraft>>
  {
    return m_spacecrafts;
  }

  /**
   * @brief
   * @return
   */
  [[nodiscard]] auto get_current_time() const -> double
  {
    return m_t;
  }

  /**
   * @brief
   * @param dt
   */
  void simulate_by(double dt)
  {
    m_t = m_propagator.propagate_by(dt);
  }

  /**
   * @brief
   * @param dt
   * @return
   */
  double simulate_to(double dt)
  {
    // fmt::print("Simulating to: {}\n", dt);
    m_t = m_propagator.propagate_to(dt);
    return m_t;
  }
};
}

#endif //SYSTEM_H
