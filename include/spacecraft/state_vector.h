//
// Created by alex_ on 5/30/2024.
//

#ifndef VECTOR_H
#define VECTOR_H
#include <memory>
#include <fmt/ostream.h>
#include <boost/numeric/odeint.hpp>

class state_vector
{
  arma::vec3 m_position;
  arma::vec3 m_velocity;
  arma::vec3 m_acceleration;


public:
  state_vector() = default;
  explicit state_vector(const arma::vec3& position): m_position(position){}
  state_vector(const arma::vec3& position, const arma::vec3& velocity): m_position(position), m_velocity(velocity){}
  state_vector(const arma::vec3&position, const arma::vec3& velocity, const arma::vec3& acceleration): m_position(position), m_velocity(velocity), m_acceleration(acceleration){}
  explicit state_vector(std::vector<double>& position)
  {
    if (position.size() != 3) {
      throw std::runtime_error("Position vector must be size 3");
    }
    m_position = {position[0], position[1], position[2]};
  }
  state_vector(const state_vector&) = default;
  state_vector& operator=(const state_vector&) = default;
  ~state_vector() = default;
  state_vector(state_vector&& other) noexcept // move constructor
    : m_position(std::move(other.m_position)),
      m_velocity(std::move(other.m_velocity)),
      m_acceleration(std::move(other.m_acceleration)) {}
  state_vector& operator=(state_vector&& other) noexcept // move assignment
{
    std::swap(m_position, other.m_position);
    std::swap(m_velocity, other.m_velocity);
    std::swap(m_acceleration, other.m_acceleration);
    return *this;
}

  state_vector(std::vector<double>& position, std::vector<double>& velocity)
  {
    if (position.size() != 3) {
      throw std::runtime_error("Position vector must be size 3");
    }
    if (velocity.size() != 3) {
      throw std::runtime_error("Velocity vector must be size 3");
    }
    m_position = {position[0], position[1], position[2]};
    m_velocity = {velocity[0], velocity[1], velocity[2]};
  }

  state_vector(std::vector<double>& position, std::vector<double>& velocity, std::vector<double>& acceleration)
  {
    if (position.size() != 3) {
      throw std::runtime_error("Position vector must be size 3");
    }
    if (velocity.size() != 3) {
      throw std::runtime_error("Velocity vector must be size 3");
    }
    if (acceleration.size() != 3) {
      throw std::runtime_error("Acceleration vector must be size 3");
    }
    m_position = {position[0], position[1], position[2]};
    m_velocity = {velocity[0], velocity[1], velocity[2]};
    m_acceleration = {acceleration[0], acceleration[0], acceleration[0]};
  }

  auto get_position() -> arma::vec3&
  {
    return m_position;
  }

  auto get_position() const  -> arma::vec3
  {
    return m_position;
  }

  auto get_velocity() -> arma::vec3&
  {
    return m_velocity;
  }

  auto get_velocity() const  -> arma::vec3
  {
    return m_velocity;
  }

  auto get_acceleration() -> arma::vec3
  {
    return m_acceleration;
  }

  auto set_position(arma::vec3& position)
  {
    m_position = position;
  }

  auto set_velocity(arma::vec3& velocity)
  {
    m_velocity = velocity;
  }

  auto set_acceleration(arma::vec3& acceleration)
  {
    m_acceleration = acceleration;
  }

  auto operator*(const double d) const
  {
    state_vector new_sv(*this);
    new_sv.m_position *= d;
    new_sv.m_velocity *= d;
    new_sv.m_acceleration *= d;

    return new_sv;
  }

  auto operator*=(const double d)
  {
    m_position *= d;
    m_velocity *= d;
    m_acceleration *= d;

    return *this;
  }


  auto operator+(const state_vector& other) const
  {
    state_vector new_sv(*this);
    new_sv.m_position += other.m_position;
    new_sv.m_velocity += other.m_velocity;
    new_sv.m_acceleration += other.m_acceleration;

    return new_sv;
  }

  auto operator+=(const state_vector& other)
  {
    m_position += other.m_position;
    m_velocity += other.m_velocity;
    m_acceleration += other.m_acceleration;

    return *this;
  }

  friend auto operator<<(std::ostream &output, const state_vector &state ) -> std::ostream & {
    output << "Position: " << state.m_position << "\n";
    output << "Velocity: " << state.m_velocity << "\n";
    output << "Acceleration: " << state.m_acceleration;
    return output;
  }

  friend auto operator*(const double d, const state_vector &state ) -> state_vector {
    auto result = state * d;
    return result;
  }

};


class combined_state_vector
{
  std::vector<std::shared_ptr<double>> m_positions;
  std::vector<std::shared_ptr<double>> m_velocities;
  std::vector<std::shared_ptr<double>> m_accelerations;

public:
  combined_state_vector()
      : m_positions({})
      , m_velocities({})
      , m_accelerations({})
  {}

  combined_state_vector(const combined_state_vector& other)
    : m_positions(other.m_positions)
    , m_velocities(other.m_velocities)
    , m_accelerations(other.m_accelerations)
  {}

  auto size() const -> size_t
  {
    return m_positions.size();
  }

  auto resize(size_t n)
  {
    m_positions.resize(n, std::make_shared<double>(0.0));
    m_velocities.resize(n, std::make_shared<double>(0.0));
    m_accelerations.resize(n, std::make_shared<double>(0.0));

  }
  void add_state_vector(state_vector& other)
  {
    // for (int i = 0; i < 3; ++i) {
    //   m_positions.push_back(other.get_position()[i]);
    //   m_velocities.push_back(other.get_velocity()[i]);
    //   m_accelerations.push_back(other.get_acceleration()[i]);
    // }
  }

  void set_positions(const std::vector<double>& positions)
  {
    std::vector<std::shared_ptr<double>> pos(positions.size());
    for (int i = 0; i < positions.size(); ++i) {
      pos[i] = std::make_shared<double>(positions[i]);
    }
    m_positions = pos;
  }

  void set_positions(const std::vector<std::shared_ptr<double>>& positions)
  {
    m_positions = positions;
  }

  void set_velocities(const std::vector<double>& velocities)
  {
    std::vector<std::shared_ptr<double>> vel(velocities.size());
    for (int i = 0; i < velocities.size(); ++i) {
      vel[i] = std::make_shared<double>(velocities[i]);
    }
    m_velocities = vel;
  }

  void set_velocities(const std::vector<std::shared_ptr<double>>& velocities)
  {
    m_positions = velocities;
  }

  void set_accelerations(const std::vector<double>& accelerations)
  {
    std::vector<std::shared_ptr<double>> acc(accelerations.size());
    for (int i = 0; i < accelerations.size(); ++i) {
      acc[i] = std::make_shared<double>(accelerations[i]);
    }
    m_accelerations = acc;
  }

  void set_accelerations(const std::vector<std::shared_ptr<double>>& accelerations)
  {
    m_positions = accelerations;
  }

  auto get_positions() const -> std::vector<std::shared_ptr<double>>
  {
    return m_positions;
  }

  auto get_velocities() const -> std::vector<std::shared_ptr<double>>
  {
    return m_velocities;
  }

  auto get_accelerations() const -> std::vector<std::shared_ptr<double>>
  {
    return m_accelerations;
  }

  auto operator*(const double d) const
  {
    combined_state_vector new_csv(*this);
    for (std::shared_ptr<double>& val : new_csv.get_positions()) {
      *val *= d;
    }
    for (std::shared_ptr<double>& val: new_csv.get_velocities()) {
      *val *= d;
    }
    for (std::shared_ptr<double>& val : new_csv.get_accelerations()) {
      *val *= d;
    }

    return new_csv;
  }

  auto operator+(const double d) const
  {
    combined_state_vector new_csv(*this);
    for (std::shared_ptr<double>& val : new_csv.get_positions()) {
      *val += d;
    }
    for (std::shared_ptr<double>& val : new_csv.get_velocities()) {
      *val += d;
    }
    for (std::shared_ptr<double>& val : new_csv.get_accelerations()) {
      *val += d;
    }

    return new_csv;
  }

  auto operator+=(const double d) -> combined_state_vector&
  {
    for (int i = 0; i < m_positions.size(); ++i) {
      *m_positions[i] += d;
      *m_velocities[i] += d;
      *m_accelerations[i] += d;
    }

    return *this;
  }

  auto operator+(const combined_state_vector& other) const
  {
    combined_state_vector new_csv(*this);
    for (int i = 0; i < m_positions.size(); ++i) {
      if (other.m_positions[i] == nullptr) {
        *new_csv.m_positions[i] += 0.0;
      } else {
        *new_csv.m_positions[i] += *other.m_positions[i];
      }

      if (other.m_velocities[i] == nullptr) {
        *new_csv.m_velocities[i] += 0.0;
      } else {
        *new_csv.m_velocities[i] += *other.m_velocities[i];

      }

      if (other.m_accelerations[i] == nullptr) {
        *new_csv.m_accelerations[i] += 0.0;
      } else {
        *new_csv.m_accelerations[i] += *other.m_accelerations[i];
      }
    }
    return new_csv;
  }

  auto operator+=(const combined_state_vector& other) -> combined_state_vector&
  {

    for (int i = 0; i < m_positions.size(); ++i) {
      *m_positions[i] += *other.m_positions[i];
      *m_velocities[i] += *other.m_velocities[i];
      *m_accelerations[i] += *other.m_accelerations[i];
    }

    return *this;
  }

  friend auto operator*(const double d, const combined_state_vector& csv) -> combined_state_vector
  {
    combined_state_vector new_csv(csv);
    for (int i = 0; i < new_csv.m_positions.size(); ++i) {
      *new_csv.m_positions[i] *= d;
      *new_csv.m_velocities[i] *= d;
      *new_csv.m_accelerations[i] *= d;
    }

    return new_csv;
  }

  friend auto operator*=(const double d, combined_state_vector& csv) -> combined_state_vector
  {
    for (int i = 0; i < csv.m_positions.size(); ++i) {
      *csv.m_positions[i] *= d;
      *csv.m_velocities[i] *= d;
      *csv.m_accelerations[i] *= d;
    }

    return csv;
  }
};

namespace boost {
namespace numeric {
namespace odeint {
template<>
struct vector_space_norm_inf<state_vector> {
  typedef double result_type;

  result_type operator()(state_vector sv) const {
    auto r = sv.get_position();
    auto v = sv.get_velocity();
    auto s = arma::join_cols(r, v);
    return arma::norm(s, "inf");
  }
};
}}}

// typedef combined_state_vector state_type;
//
// namespace boost { namespace numeric { namespace odeint {
//
// template< >
// struct is_resizeable< state_type >
// { // declare resizeability
//   typedef boost::true_type type;
//   const static bool value = type::value;
// };
//
// template< >
// struct same_size_impl< state_type , state_type >
// { // define how to check size
//   static bool same_size( const state_type &v1 ,
//                          const state_type &v2 )
//   {
//     return v1.size() == v2.size();
//   }
// };
//
// template< >
// struct resize_impl< state_type , state_type >
// { // define how to resize
//   static void resize( state_type &v1 ,
//                       const state_type &v2 )
//   {
//     v1.resize( v2.size() );
//   }
// };
//
// } } }

#endif //VECTOR_H
