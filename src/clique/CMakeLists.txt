include_directories(../graph ../tools)


add_library(kCListClique kCListClique_v6.cpp)
target_link_libraries(kCListClique Graph linearSet types)
