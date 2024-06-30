//
// Created by alex_ on 6/10/2024.
//

#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <naomi.h>

#include "boost/numeric/odeint.hpp"
#include "propagators/event_detector.h"
#include "spacecraft/state_vector.h"

template< class TStepper>
class integrator
{
    TStepper m_stepper;

public:
    ~ integrator() = default;
    integrator() = default;
    integrator(const integrator&) = delete;
    integrator(integrator&&) = delete;
    auto operator=(const integrator&) -> integrator& = delete;
    auto operator=(const integrator&&) -> integrator& = delete;

    void integrate(std::function<void(state_vector&, state_vector&, double)> system, state_vector& state, double start_time, double end_time, double step_size)
    {
        auto stepper = make_dense_output(1.0e-6, 1.0e-6, m_stepper);
        auto ode_range = boost::numeric::odeint::make_adaptive_range(stepper, system, state,
                                                                 start_time, end_time, step_size);
        auto cond = apside_detector();
        auto found_iter = std::find_if(ode_range.first, ode_range.second, cond);
        if(found_iter == ode_range.second)
        {
          // no threshold crossing -> return time after t_end and ic
          std::cout << "No apside crossing\n";
        }
        double t0 = stepper.previous_time();
        double t1 = stepper.current_time();
        double t_m;
        state_type x_m;
        // use odeint's resizing functionality to allocate memory for x_m
        boost::numeric::odeint::adjust_size_by_resizeability(x_m, state,
                                                             typename boost::numeric::odeint::is_resizeable<state_type>::type());
        while(std::abs(t1 - t0) > 1e-6) {
          t_m = 0.5 * (t0 + t1);  // get the mid point time
          stepper.calc_state(t_m, x_m); // obtain the corresponding state
          if (cond(x_m))
            t1 = t_m;  // condition changer lies before midpoint
          else
            t0 = t_m;  // condition changer lies after midpoint
        }
        // we found the interval of size eps, take it's midpoint as final guess
        t_m = 0.5 * (t0 + t1);
        stepper.calc_state(t_m, x_m);
        // integrate_const(m_stepper, system, state, start_time, end_time, step_size);
    }
};

// template<class System, class Condition>
// std::pair<double, state_type>
// find_condition(state_type &x0, System sys, Condition cond,
//                const double t_start, const double t_end, const double dt,
//                const double precision = 1E-6) {
//
//     // integrates an ODE until some threshold is crossed
//     // returns time and state at the point of the threshold crossing
//     // if no threshold crossing is found, some time > t_end is returned
//
//     auto stepper = make_dense_output(1.0e-6, 1.0e-6,
//                                                              boost::numeric::odeint::runge_kutta_dopri5<state_type>());
//
//     auto ode_range = boost::numeric::odeint::make_adaptive_range(std::ref(stepper), sys, x0,
//                                                                  t_start, t_end, dt);
//
//     // find the step where the condition changes
//     auto found_iter = std::find_if(ode_range.first, ode_range.second, cond);
//
//     if(found_iter == ode_range.second)
//     {
//         // no threshold crossing -> return time after t_end and ic
//         return std::make_pair(t_end + dt, x0);
//     }
//
//     // the dense out stepper now covers the interval where the condition changes
//     // improve the solution by bisection
//     double t0 = stepper.previous_time();
//     double t1 = stepper.current_time();
//     double t_m;
//     state_type x_m;
//     // use odeint's resizing functionality to allocate memory for x_m
//     boost::numeric::odeint::adjust_size_by_resizeability(x_m, x0,
//                                                          typename boost::numeric::odeint::is_resizeable<state_type>::type());
//     while(std::abs(t1 - t0) > precision) {
//         t_m = 0.5 * (t0 + t1);  // get the mid point time
//         stepper.calc_state(t_m, x_m); // obtain the corresponding state
//         if (cond(x_m))
//             t1 = t_m;  // condition changer lies before midpoint
//         else
//             t0 = t_m;  // condition changer lies after midpoint
//     }
//     // we found the interval of size eps, take it's midpoint as final guess
//     t_m = 0.5 * (t0 + t1);
//     stepper.calc_state(t_m, x_m);
//     return std::make_pair(t_m, x_m);
// }




#endif //INTEGRATOR_H
