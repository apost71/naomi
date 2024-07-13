//
// Created by alex on 7/12/2024.
//

#ifndef ROTATION_H
#define ROTATION_H
#include <memory>
#include <armadillo>

namespace naomi::attitude
{
class rotation
{
public:
  virtual ~rotation() = default;
  virtual std::shared_ptr<rotation> apply_to(const std::shared_ptr<rotation>& r) = 0;
  virtual arma::vec3 apply_to(const arma::vec3& r) = 0;
  virtual arma::mat33 get_dcm() = 0;

};
}
#endif //ROTATION_H
