cmake_minimum_required(VERSION 3.25)
project(naomi-simulation LANGUAGES CXX)

string(COMPARE EQUAL "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}" is_top_level)

if(is_top_level)
    find_package(naomi REQUIRED)
endif()

set(sources main.cpp
        two_body/simple_two_body_propagation.h)
source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" FILES ${sources})

add_executable(naomi-simulation ${sources})
#target_sources(naomi-simulation PRIVATE ${sources})
target_link_libraries(naomi-simulation PRIVATE naomi)
