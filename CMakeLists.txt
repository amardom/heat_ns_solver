cmake_minimum_required(VERSION 3.20)
project(heat_ns_solver LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(solver src/solver.cpp)
target_include_directories(solver PUBLIC src)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(solver PRIVATE -Wall -Wextra -Wpedantic)
endif()
add_executable(simulate src/main.cpp)
target_link_libraries(simulate PRIVATE solver)

include(CTest)
if(BUILD_TESTING)
    find_package(GTest QUIET)
    if(NOT GTest_FOUND)
        include(FetchContent)
        FetchContent_Declare(googletest
            URL https://github.com/google/googletest/archive/refs/tags/v1.15.2.tar.gz)
        set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
        FetchContent_MakeAvailable(googletest)
    endif()
    add_executable(solver_tests
        tests/heat_equation_test.cpp
        tests/fluid_test.cpp
        tests/poiseuille_test.cpp
        tests/couette_test.cpp
        tests/heated_wall_test.cpp
        tests/stokes_test.cpp
        tests/stability_test.cpp)
    target_link_libraries(solver_tests PRIVATE solver GTest::gtest_main)
    include(GoogleTest)
    gtest_discover_tests(solver_tests)
endif()
