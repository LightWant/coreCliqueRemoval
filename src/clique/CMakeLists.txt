include_directories(../graph ../tools)


add_library(kCListClique kCListClique_v6.cpp)
target_link_libraries(kCListClique Graph linearSet types)


add_library(parallelKcc parallelkCListClique.h parallelkCListClique.cpp)
target_link_libraries(parallelKcc Graph linearSet types)
