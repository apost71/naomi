//
// Created by alex on 6/30/2024.
//
/*
 * adaptive_iterator.cpp
 *
 * Copyright 2012-2013 Karsten Ahnert
 * Copyright 2012 Mario Mulansky
 *
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or
 * copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#include <algorithm>
#include <array>
#include <cassert>
#include <iostream>
#include <iterator>
#include <utility>

#include <boost/numeric/odeint/iterator/adaptive_iterator.hpp>
#include <boost/numeric/odeint/iterator/adaptive_time_iterator.hpp>
#include <boost/numeric/odeint/stepper/generation.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta_cash_karp54.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta_dopri5.hpp>
#include <boost/range/adaptor/filtered.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range/numeric.hpp>
#include <gtest/gtest.h>

#define tab "\t"

using namespace std;
using namespace boost::numeric::odeint;

const double sigma = 10.0;
const double R = 28.0;
const double b = 8.0 / 3.0;

struct lorenz
{
    template< class State , class Deriv >
    void operator()( const State &x , Deriv &dxdt , double t ) const
    {
        dxdt[0] = sigma * ( x[1] - x[0] );
        dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
        dxdt[2] = -b * x[2] + x[0] * x[1];
    }
};

#include <typeinfo>

TEST(BoostAdaptiveRange, TestAdaptiveRange)
{
    typedef std::array< double , 3 > state_type;

    /*
     * Controlled steppers with time iterator
     */

    // // std::for_each
    // {
    //     auto stepper = make_controlled( 1.0e-6 , 1.0e-6 , runge_kutta_cash_karp54< state_type >() );
    //     state_type x = {{ 10.0 , 10.0 , 10.0 }};
    //     std::for_each( make_adaptive_time_iterator_begin( stepper , lorenz() , x , 0.0 , 1.0 , 0.01 ) ,
    //                    make_adaptive_time_iterator_end( stepper , lorenz() , x ) ,
    //                    []( const std::pair< const state_type&, double > &x ) {
    //                        std::cout << x.second << tab << x.first[0] << tab << x.first[1] << tab << x.first[2] << "\n"; } );
    // }

  {
      auto stepper = make_controlled( 1.0e-6 , 1.0e-6 , runge_kutta_cash_karp54< state_type >() );
      state_type x = {{ 10.0 , 10.0 , 10.0 }};
      auto iter = boost::find_if( make_adaptive_time_range( stepper , lorenz() , x , 0.0 , 1.0 , 0.01 ) ,
                                  []( const std::pair< const state_type & , double > &x ) {
                                      return ( x.first[0] < 0.0 ); } );
      cout << iter->second << "\t" << iter->first[0] << "\t" << iter->first[1] << "\t" << iter->first[2] << "\n";
  }


}