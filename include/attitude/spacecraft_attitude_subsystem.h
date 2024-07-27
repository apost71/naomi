//
// Created by alex_ on 7/24/2024.
//

#ifndef ATTITUDE_SUBSYSTEM_H
#define ATTITUDE_SUBSYSTEM_H
#include "spacecraft/spacecraft_subsystem.h"

class spacecraft_attitude_subsystem: public spacecraft_subsystem
{
  std::shared_ptr<attitude::attitude_provider> _provider;
public:
  spacecraft_attitude_subsystem(
    const std::shared_ptr<attitude::attitude_provider>& provider)
      : spacecraft_subsystem("AttitudeSubsystem")
        , _provider(provider)
  {
  }
  std::vector<std::shared_ptr<attitude::additional_state_provider>>
  get_additional_state_providers() override;
  void initialize(const spacecraft_state& state) override;
  void update(const spacecraft_state& state) override;
  void terminate(const spacecraft_state& state) override;
  ~spacecraft_attitude_subsystem() override;
};
#endif //ATTITUDE_SUBSYSTEM_H
