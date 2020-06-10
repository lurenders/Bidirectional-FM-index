cmake_minimum_required(VERSION 2.6)
project(BWT)

add_executable(BWT main.cpp BWT.cpp BiBWT.cpp customtypedefs.cpp)
add_executable(BWT_build buildIndex.cpp customtypedefs.cpp)
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall -pedantic -mpopcnt -O3")

install(TARGETS BWT RUNTIME DESTINATION bin)
install(TARGETS BWT_build RUNTIME DESTINATION bin)


