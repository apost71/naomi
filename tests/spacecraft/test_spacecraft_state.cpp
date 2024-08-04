//
// Created by alex_ on 7/28/2024.
//

#include <armadillo>

#include <gtest/gtest.h>

#include "attitude/constant_attitude_provider.h"
#include "attitude/torque_free_provider.h"
#include "spacecraft/pv_coordinates_provider.h"
#include "spacecraft/spacecraft_state.h"

typedef std::shared_ptr<attitude_provider> att_provider_ptr;
typedef std::shared_ptr<state_provider> state_provider_ptr;

using namespace std;

TEST(TestSpacecraftState, TestIntegratedStates)
{
  att_provider_ptr attitude_prov = make_shared<constant_attitude_provider>();
  state_provider_ptr state_prov = make_shared<pv_coordinates_provider>(pv_coordinates(vector_type({1, 2, 3, 4, 5, 6})));
  auto state = spacecraft_state(state_prov, attitude_prov, 100.0);
  const auto vec = state.get_integrated_state();
  EXPECT_EQ(vec.size(), 9);
}

TEST(TestSpacecraftState, TestIntegratedStatesAttitude)
{
  const vector_type state_vec = {1, 2, 3, 4, 5, 6};
  const att_provider_ptr attitude_prov =
      make_shared<torque_free_attitude_provider>(torque_free_attitude_provider(
          arma::mat33(), {1, 2, 3, 4}, pv_coordinates(state_vec)));
  const state_provider_ptr state_prov = make_shared<pv_coordinates_provider>(pv_coordinates(vector_type({1, 2, 3, 4, 5, 6})));
  auto state = spacecraft_state(state_prov, attitude_prov, 100.0);
  auto vec = state.get_integrated_state();
  EXPECT_EQ(vec.size(), 19);
}

TEST(TestSpacecraftState, TestSetIntegratedState)
{
  const vector_type state_vec = {1, 2, 3, 4, 5, 6};
  const att_provider_ptr attitude_prov =
      make_shared<torque_free_attitude_provider>(torque_free_attitude_provider(
          arma::mat33(), {1, 2, 3, 4}, pv_coordinates(state_vec)));
  const state_provider_ptr state_prov = make_shared<pv_coordinates_provider>(pv_coordinates(vector_type({1, 2, 3, 4, 5, 6})));
  auto state = spacecraft_state(state_prov, attitude_prov, 100.0);
  const vector_type expected_vec = {2, 3, 4, 5, 6, 7, 1, 1, 1, 2, 3, 4, 5, 2, 3, 4, 1, 1, 1};
  state.set_integrated_state(expected_vec);
  const auto updated_vec = state.get_integrated_state();
  EXPECT_TRUE(arma::approx_equal(updated_vec, expected_vec, "absdiff", 1e-6));
}
