#include "../graph/graph.h"
#include "../tools/getArgs.hpp"
#include "../MClique/parallelCCRMCE.h"
#include <iostream>
#include <chrono>
#include <iomanip>
#include <omp.h>

void printUsage() {
    std::cout << "bin/run " << std::endl;
    std::cout << "-f graph file directory(edge.bin & idx.bin)" << std::endl;
    std::cout << "-f_txt graph file text file, each edge exists one time" << std::endl;
    std::cout << "-f_txtD graph file text file, each edge exists two times" << std::endl;
    std::cout << "-k" << std::endl;
}

int main(int argc, char * argv[])
{
    argsController ac(argc, argv);

    if(!ac.exist("-f_txt") && !ac.exist("-f") && !ac.exist("-f_txtD")) {
        printUsage();
        return 0;
    }
    
    bool noVUM = false;
    if(ac.exist("noUVM")) noVUM = true;

    Graph g;
    if(ac.exist("-f_txt")) g.readFromText(ac["-f_txt"], noVUM);
    else if(ac.exist("-f_txtD")) g.readFromTextDoubleEdges(ac["-f_txtD"]);
    else g.readFromBin(ac["-f"]);


    std::cout << "load graph: n " << g.n << " m " << g.m << " maxD " << g.maxD << std::endl;

    auto s1 = std::chrono::steady_clock::now();

    g.changeToCoreOrder();

    auto s2 = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(s2 - s1);
    std::cout << "changeToCoreOrder:" << duration.count() << "ms" << std::endl;
    std::cout << "coreNumber:" << g.coreNumber << std::endl;


    
    // g.initHash();

    

    int numT = 5;
    if(ac.exist("-t")) numT =  std::stoi(ac["-t"]);
    std::cout << "numThreads:" << numT << std::endl;

    omp_set_num_threads(numT);

    double tS = omp_get_wtime();

    parallelMCE pM(std::move(g), numT);
    pM.useHash();

    double tEInit = omp_get_wtime();
    std::cout << "initTime:" << tEInit - tS << "s" << std::endl;
    
    ui cnt = 0;
    
    if(ac.exist("-e")) {
        cnt = pM.runCoreCliqueRemovalV3EdgeParallel();
    }
    else if(ac.exist("-steal")) {
        cnt = pM.runCoreCliqueRemovalV3WorkStealing();
    }
    else {
        cnt = pM.runCoreCliqueRemovalV3();
    }

    double tE = omp_get_wtime();
    std::cout << "time:" << tE - tEInit << "s" << std::endl;

    std::cout << "Mclique:" << cnt << std::endl;
    
    return 0;
}