cmake_minimum_required(VERSION 2.8)

set(CMAKE_CXX_FLAGS  "-std=c++11 -DNDEBUG")

include_directories(/home/piyush/ETH/work/Code/BEM/2DParametricBEM/include)
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(/home/piyush/ETH/work/Code/BEM/CppHilbert/Library/source)

#add_subdirectory(/home/piyush/ETH/work/Code/BEM/2DParametricBEM ${CMAKE_BINARY_DIR}/external)
#link_directories(${CMAKE_BINARY_DIR}/external)

#add_subdirectory(/home/piyush/ETH/work/Code/BEM/CppHilbert/Library ${CMAKE_BINARY_DIR}/CppHilbert)
#link_directories(${CMAKE_BINARY_DIR}/CppHilbert)

add_library(robustbem shifted_logweight_quadrature.cpp parabolic_parametrization.cpp line_near_singularity.cpp point_near_singularity.cpp)
