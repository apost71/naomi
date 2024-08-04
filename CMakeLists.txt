cmake_minimum_required(VERSION 3.22)
project(naomi C CXX)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")


include(GNUInstallDirs)


find_package(fmt)
find_package(Boost)
find_package(Eigen3)
find_package(SymEngine 0.1.0 REQUIRED CONFIG
        PATH_SUFFIXES lib/cmake/symengine CMake/)
find_package(Armadillo)
#set(CMAKE_BUILD_TYPE ${SYMENGINE_BUILD_TYPE})
set(CMAKE_CXX_FLAGS_RELEASE ${SYMENGINE_CXX_FLAGS_RELEASE})
set(CMAKE_CXX_FLAGS_DEBUG ${SYMENGINE_CXX_FLAGS_DEBUG})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SYMENGINE_CXX_FLAGS} -std=c++17")

add_library(naomi src/naomi.cpp
        include/bodies/celestial_body.h
        include/bodies/earth.h
        include/constants.h
        include/integrators/integrator.h
        include/spacecraft/spacecraft.h
        include/propagators/numerical_propagator.h
        include/orbits/orbits.h
        include/orbits/keplerian.h
        include/orbits/cartesian.h
        include/systems/system.h
        include/maneuvers/hohmann_transfer.h
        include/maneuvers/maneuver.h
        include/frames/transforms.h
        include/propagators/event_detector.h
        include/forces/force_model.h
        include/forces/two_body_force_model.h
        include/propagators/event_handler.h
        include/maneuvers/maneuver_plan.h
        include/control/controller.h
        include/observers/simulation_observer.h
        include/simulation/simulation.h
        src/observers/results_csv_writer_observer.cpp
        include/forces/two_body_rot_force_model.h
        include/attitude/attitude_provider.h
        include/attitude/torque_free.h
        include/math/vector_utils.h
        include/math/volume_integrals.h
        include/attitude/rotation.h
        include/attitude/euler_angles.h
        include/spacecraft/body_shape.h
        include/control/pid_controller.h
        include/math/quaternion.h
        include/spacecraft/spacecraft_state.h
        include/math/transport_theorem.h
        include/spacecraft/pv_coordinates.h
        include/frames/frame.h
        include/spacecraft/spacecraft_subsystem.h
        include/attitude/spacecraft_attitude_subsystem.h
        include/attitude/angular_coordinates.h
        include/spacecraft/state_provider.h
        include/spacecraft/pv_coordinates_provider.h
        include/attitude/torque_free_provider.h
        include/attitude/constant_attitude_provider.h
        include/attitude/constant_attitude_provider.h)
target_compile_features(naomi PUBLIC cxx_std_17)
include_directories(${SYMENGINE_INCLUDE_DIRS})
target_include_directories(naomi PUBLIC include )
target_link_libraries(naomi
        fmt::fmt
        Eigen3::Eigen
        ${SYMENGINE_LIBRARIES}
        ${ARMADILLO_LIBRARIES}
        ${BOOST_LIBRARIES}
)

set_target_properties(naomi PROPERTIES PUBLIC_HEADER "include/naomi.h")

if(BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

if (NOT BUILD_TESTING STREQUAL OFF)
    add_subdirectory(tests)
endif()

install(TARGETS naomi)

