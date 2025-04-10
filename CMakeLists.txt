cmake_minimum_required(VERSION 3.6)

project(ccr)

SET(CMAKE_CXX_COMPILER "g++")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -std=c++14 -fopenmp -march=native -mavx -O3")
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


set(exName21 "KCCparallel") #parallel ccrkce clique
add_executable(${exName21} src/run/parallelKCC.cpp)
target_link_libraries(${exName21} getArgs)
target_link_libraries(${exName21} Graph)
target_link_libraries(${exName21} parallelKcc)


set(exName22 "MCEparallel") #parallel ccrmce clique
add_executable(${exName22} src/run/parallelMCE.cpp)
target_link_libraries(${exName22} getArgs)
target_link_libraries(${exName22} Graph)
target_link_libraries(${exName22} parallelMCE)