#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>
#include "boost/numeric/odeint/stepper/runge_kutta4.hpp"
#include <vector>
#include <string>



#ifdef _WIN32
  #define NAOMI_EXPORT __declspec(dllexport)
#else
  #define NAOMI_EXPORT
#endif

class state_vector;
NAOMI_EXPORT void naomi();
NAOMI_EXPORT void naomi_print_vector(const std::vector<std::string> &strings);

typedef state_vector state_type;