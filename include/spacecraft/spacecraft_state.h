//
// Created by alex on 7/21/2024.
//

#ifndef SPACECRAFT_STATE_H
#define SPACECRAFT_STATE_H
#include <utility>

#include "attitude/attitude_provider.h"
#include "state_provider.h"

namespace naomi
{
using namespace attitude;

class mass_provider
{

};

class spacecraft_state
{
  std::shared_ptr<state_provider> _state_provider;
  std::shared_ptr<attitude_provider> _attitude_provider;
  double _mass;
  std::vector<
      std::pair<
        arma::span,
        std::shared_ptr<integrated_provider>
      >
    > _integrated_state_idxs = {};

  std::vector<std::shared_ptr<integrated_provider>> get_integrated_providers()
      const
  {
    std::vector<std::shared_ptr<integrated_provider>> integrated_providers;
    const std::shared_ptr<integrated_provider> int_state_provider =
        std::dynamic_pointer_cast<integrated_provider>(_state_provider);
    const std::shared_ptr<integrated_provider> int_attitude_provider =
        std::dynamic_pointer_cast<integrated_provider>(_attitude_provider);

    if (int_state_provider != nullptr) {
      integrated_providers.push_back(int_state_provider);
    }
    if (int_attitude_provider != nullptr) {
      integrated_providers.push_back(int_attitude_provider);
    }
    return integrated_providers;
  }

public:
  spacecraft_state(
      const std::shared_ptr<state_provider>& state_provider,
      const std::shared_ptr<attitude_provider>& attitude_provider,
      const double mass)
      : _state_provider(state_provider)
      , _attitude_provider(attitude_provider)
      , _mass(mass)
  {
    _integrated_state_idxs = get_provider_mapping();
    _integrated_state_idxs = get_provider_mapping();
  }

  std::shared_ptr<state_provider> get_state_provider()
  {
    return _state_provider;
  }

  std::shared_ptr<attitude_provider> get_attitude_provider()
  {
    return _attitude_provider;
  }

  std::vector<std::pair<arma::span, std::shared_ptr<integrated_provider>>> get_provider_mapping()
  {
    auto integrated_providers = get_integrated_providers();
    std::size_t start_idx = 0;
    state_type state;
    std::vector<
      std::pair<
        arma::span,
        std::shared_ptr<integrated_provider>
      >
    > new_idxs;
    for (const auto& provider : integrated_providers) {
      const auto size = provider->get_size();
      const auto end_idx = start_idx + size;
      const auto prov_span = arma::span(start_idx, end_idx - 1);
      new_idxs.emplace_back(prov_span, provider);
      start_idx = end_idx;
    }
    _integrated_state_idxs = new_idxs;
    return _integrated_state_idxs;
  }

  state_type get_integrated_state()
  {
    state_type state;
    for (const auto& [spn, prv] : _integrated_state_idxs) {
      auto int_state = prv->get_integrated_state();
      state.insert_rows(spn.a, int_state);
    }
    return state;
  }

  void set_integrated_state(const state_type& state)
  {
    for (const auto& [spn, prv] : _integrated_state_idxs) {
      prv->set_integrated_state(state(spn));
    }
  }
};

}
#endif //SPACECRAFT_STATE_H
