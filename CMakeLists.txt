cmake_minimum_required(VERSION 3.16)
project(cpp_pricer CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(pricer INTERFACE)
target_include_directories(pricer INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)

enable_testing()
add_executable(test_pricer tests/test_pricer.cpp)
target_link_libraries(test_pricer PRIVATE pricer)
add_test(NAME test_pricer COMMAND test_pricer)

add_executable(pricer_run src/main.cpp)
target_link_libraries(pricer_run PRIVATE pricer)
