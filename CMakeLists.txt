cmake_minimum_required(VERSION 3.28)
project(naomi CXX)


find_package(fmt)
find_package(Boost)
find_package(Eigen3)
find_package(symengine)
find_package(Armadillo)

add_library(naomi src/naomi.cpp
        include/bodies/celestial_body.h
        include/bodies/earth.h
        include/constants.h
        include/integrators/integrator.h
        include/spacecraft/spacecraft.h
        include/spacecraft/state_vector.h
        include/propagators/numerical_propagator.h
        include/orbits/orbits.h
        include/orbits/keplerian.h
        include/orbits/cartesian.h
        include/systems/system.h
        include/maneuvers/hohmann_transfer.h
        include/maneuvers/maneuver.h
        include/frames/transforms.h
        include/propagators/event_detector.h)
target_include_directories(naomi PUBLIC include)
target_link_libraries(naomi
        fmt::fmt
        Eigen3::Eigen
        symengine::symengine
        Armadillo::Armadillo
        Boost::headers
        Boost::atomic
        Boost::container
        Boost::context
        Boost::coroutine
        Boost::date_time
        Boost::exception
        Boost::filesystem
        Boost::graph
        Boost::iostreams
        Boost::locale
        Boost::log
        Boost::log_setup
        Boost::math_c99
        Boost::math_c99f
        Boost::math_c99l
        Boost::math_tr1
        Boost::math_tr1f
        Boost::math_tr1l
        Boost::nowide
        Boost::program_options
        Boost::random
        Boost::regex
        Boost::serialization
        Boost::system
        Boost::thread
        Boost::timer
        Boost::type_erasure
        Boost::wave
        Boost::wserialization
)

set_target_properties(naomi PROPERTIES PUBLIC_HEADER "include/naomi.h")

if (NOT BUILD_TESTING STREQUAL OFF)
    add_subdirectory(tests)
endif()

install(TARGETS naomi)
