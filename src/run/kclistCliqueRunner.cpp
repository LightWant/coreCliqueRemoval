#include "../graph/graph.h"
#include "../tools/getArgs.hpp"
#include "../clique/kCListClique.h"
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

    int k = 7;
    if(ac.exist("-k")) k =  std::stoi(ac["-k"]);
    std::cout << "k:" << k << std::endl;

    kclistClique pP(std::move(g), k);

    auto t1 = std::chrono::steady_clock::now();
    
    double cnt = 0;
    // if(ac.exist("-peel")) cnt = pP.runPeel();
    // else if(ac.exist("-peel_core")) cnt = pP.runPeelCore();
    // else if(ac.exist("-v2")) cnt = pP.runV2();
    // else if(ac.exist("-v3")) cnt = pP.runV3();
    // else if(ac.exist("-v4")) cnt = pP.run_v4();
    // else if(ac.exist("-v5")) cnt = pP.run_v5();//the weight of a node is the count of neighbors in the remaining node set except the found clique
    // else if(ac.exist("-v6")) cnt = pP.run_v6();//core-peeling based large clique
    // else if(ac.exist("-v7")) cnt = pP.run_v7();
    // else if(ac.exist("-s")) {
    //     ull sps = std::stoll(ac["-s"]);
    //     printf("SampleSize:%llu\n", sps);
    //     if(ac.exist("-sv2")) cnt = pP.runCCPath_V2SavingLC(sps);
    //     else cnt = pP.runCCPath(sps);
    // }
    // else cnt = pP.run();

    cnt = pP.run_v6();

    auto t2 = std::chrono::steady_clock::now();
    auto durationt = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cout << "time:" << durationt.count() << "ms" << std::endl;

    std::cout << std::fixed << std::setprecision(0) <<
             k << "-clique:" << cnt << std::endl;
    
    return 0;
}