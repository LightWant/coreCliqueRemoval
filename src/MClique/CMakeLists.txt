include_directories(../graph ../tools)

add_library(MCE mce.h mcecrrv3.cpp)
target_link_libraries(MCE Graph linearSet types)

add_library(parallelMCE parallelCCRMCE.h parallelCCRMCE.cpp)
target_link_libraries(parallelMCE Graph linearSet types)
