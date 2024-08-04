//
// Created by alex on 7/4/2024.
//

#ifndef FORCE_MODEL_H
#define FORCE_MODEL_H

#include "naomi.h"

namespace naomi::forces
{
class force_model
{

public:
  virtual void operator()( const vector_type& x , vector_type& dxdt, double t) const = 0;
  virtual ~force_model() = default;

};

class equations_of_motion
{

public:
  [[nodiscard]] virtual vector_type get_derivative( const vector_type& state, double t) const = 0;
  virtual ~equations_of_motion() = default;

};

}

#endif //FORCE_MODEL_H
