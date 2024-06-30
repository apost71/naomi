cmake_minimum_required(VERSION 3.28)
project(NaomiTest CXX)

find_package(GTest REQUIRED CONFIG)


add_executable(test_naomi test_naomi.cpp
        orbits/test_orbit_utils.cpp
        maneuvers/test_hohmann_transfer.cpp
        frames/test_frame_transforms.cpp)
target_link_libraries(test_naomi naomi GTest::gtest GTest::gtest_main)

include(GoogleTest)
gtest_discover_tests(test_naomi)