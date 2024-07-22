cmake_minimum_required(VERSION 3.22)
project(NaomiTest CXX)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

enable_testing()

find_package(GTest REQUIRED CONFIG)


add_executable(test_naomi test_naomi.cpp
        orbits/test_orbit_utils.cpp
        maneuvers/test_hohmann_transfer.cpp
        frames/test_frame_transforms.cpp
        boost_examples.cpp
        simulation/test_simulation.cpp
        observers/test_observer.cpp
        attitude/test_torque_free.cpp
        attitude/test_euler_angles.cpp
        spacecraft/test_body_shape.cpp)
target_link_libraries(test_naomi naomi GTest::gtest GTest::gtest_main)

include(GoogleTest)
gtest_discover_tests(test_naomi)