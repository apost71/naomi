cmake_minimum_required(VERSION 3.22)
project(naomi-simulation LANGUAGES C CXX)

string(COMPARE EQUAL "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}" is_top_level)

if(is_top_level)
    find_package(naomi REQUIRED)
endif()

set(sources simulation_main.cpp
        two_body/simple_two_body_propagation.h
        maneuvers/hohmann_transfer_example.h
)
#source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" FILES ${sources})

add_executable(naomi-simulation ${sources})
#target_sources(naomi-simulation PRIVATE ${sources})
target_link_libraries(naomi-simulation naomi)
