cmake_minimum_required(VERSION 3.6)

project(ccr)

SET(CMAKE_CXX_COMPILER "g++")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -march=native -mavx -O3")
SET(EXECUTABLE_OUTPUT_PATH "bin/")

set(exName10 "kcc") #kCListClique
set(exName20 "mce") #maximal clique

add_subdirectory(src/graph)
add_subdirectory(src/clique)
add_subdirectory(src/tools)
add_subdirectory(src/MClique)

add_executable(${exName10} src/run/kclistCliqueRunner.cpp)
target_link_libraries(${exName10} getArgs)
target_link_libraries(${exName10} Graph)
target_link_libraries(${exName10} kCListClique)

add_executable(${exName20} src/run/maximalCliqueLargeClique.cpp)
target_link_libraries(${exName20} getArgs)
target_link_libraries(${exName20} Graph)
target_link_libraries(${exName20} MCE)