#include "../graph/graph.h"
#include "../tools/getArgs.hpp"
#include "../MClique/mce.h"
#include <iostream>
#include <chrono>
#include <iomanip>

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

    MCE pM(std::move(g));

    auto t1 = std::chrono::steady_clock::now();
    
    ui cnt = 0;

    // if(ac.exist("-ccr")) {
    //     cnt = pM.runCoreCliqueRemoval();
    // }
    // else if(ac.exist("-v2")) {//delete
    //     cnt = pM.runCoreCliqueRemovalV2();
    // }
    // else if(ac.exist("-v3")) {//CoreCliqueRemoval
    //     cnt = pM.runCoreCliqueRemovalV3();
    // }
    // else if(ac.exist("-v4")) {//CoreCliqueRemoval+graph deduction
    //     cnt = pM.runCoreCliqueRemovalV4();
    // }
    // else if(ac.exist("-v5")) {//based on v3, the pivot to make T1T2 be as small as possible
    //     cnt = pM.runCoreCliqueRemovalV5();
    // }
    // else cnt = pM.run();//BK_degeneracy

    cnt = pM.runCoreCliqueRemovalV3();

    auto t2 = std::chrono::steady_clock::now();
    auto durationt = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cout << "time:" << durationt.count() << "ms" << std::endl;

    std::cout << "Mclique:" << cnt << std::endl;
    
    return 0;
}