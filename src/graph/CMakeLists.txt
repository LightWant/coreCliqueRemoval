include_directories(../tools)

add_library(Graph graph.h graph.cpp)
target_link_libraries(Graph fastIO types)
