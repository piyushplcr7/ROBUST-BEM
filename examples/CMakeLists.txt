cmake_minimum_required(VERSION 2.8)

set(CMAKE_CXX_FLAGS  "-std=c++11 -DNDEBUG")

include_directories(/home/piyush/ETH/work/Code/BEM/2DParametricBEM/include)
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(/home/piyush/ETH/work/Code/BEM/CppHilbert/Library/source)

add_subdirectory(/home/piyush/ETH/work/Code/BEM/2DParametricBEM ${CMAKE_BINARY_DIR}/external)
link_directories(${CMAKE_BINARY_DIR}/external)

add_subdirectory(/home/piyush/ETH/work/Code/BEM/CppHilbert/Library ${CMAKE_BINARY_DIR}/CppHilbert)
link_directories(${CMAKE_BINARY_DIR}/CppHilbert)

add_executable(main2 model_problem_2.cpp)
add_executable(main3 model_problem_3.cpp)
add_executable(dist distance_plotting.cpp)
add_executable(main_sol model_problem_soln.cpp)
add_executable(main_sol2 model_problem_2_soln.cpp)
add_executable(main_sol22 model_problem_2_soln_2.cpp)
add_executable(main_sol3 model_problem_3_soln.cpp)
add_executable(gentest gen_weight_test.cpp)

target_link_libraries(main2 parametrizations single_layer quadrature)
target_link_libraries(main3 robustbem parametrizations single_layer quadrature CppHilbert)
target_link_libraries(main_sol2 parametrizations single_layer quadrature robustbem)
target_link_libraries(main_sol22 parametrizations single_layer quadrature robustbem)
target_link_libraries(main_sol3 parametrizations single_layer quadrature CppHilbert robustbem)
target_link_libraries(gentest quadrature robustbem)
target_link_libraries(dist quadrature robustbem)
