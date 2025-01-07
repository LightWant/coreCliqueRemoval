#ifndef KCLISTCLIQUE_H
#define KCLISTCLIQUE_H

#include "../graph/graph.h"
#include "../tools/types.hpp"

class kclistClique {
private:
    Graph g;
    ui k = 0;
    double answer = 0;
    std::vector<std::vector<ui>> nodes;
    std::vector<std::vector<ui>> adj;
    std::vector<std::vector<ui>> edAdj;
    std::vector<ui> level;
    std::vector<ui> clique;
    std::vector<ui> cliqueNei;
    std::vector<ui> deg;

    void buildSg(ui u);
    void listing(ui deep, ui edClique);

    std::vector<ui> vis;
    void listingV2(ui deep, ui edClique);
    void listingV3(ui deep, ui edClique);

private://for ccpath and V3
    void allocaMemoryForSubgraph();
    void findALargeClique(ui u); 
    void reSet_findALargeClique(ui);

private://for v3
    void findALargeClique_greedyByDegree(ui u); 

private://for ccpath
    std::vector<ui> color, colorsOfNei;
    ui maxC;
    std::vector<ui> pIdx, eIdx;
    void coloring();
    void reOrderingByColor();

    std::vector<std::vector<ull>> dp;
    std::vector<ull> sWeight;
    double computeDp();

    std::vector<ui> timesOfU;
    ull sample(ui u, ull s, double cntPath, ull & realSz);

    std::vector<ui> LCs, pLCs;//O(m)

private://for v5, extracts multiple cliques
    ui cntOfLargeCliques;
    std::vector<std::vector<ui>> largeCliques;
    std::vector<std::vector<ui>> edLargeCliques;
    void listingV5(ui deep);

private:
    ui maxK = 31;
    ui maxSize = 3000;
    double ** CN = nullptr, *bf3 = nullptr;
    void computeC() {
        ui maxD = maxSize, maxD2 = maxK;
        CN = new double*[maxD];
        bf3 = new double[maxD2 * maxD];
        for(int i = 0; i < maxD; i++) {
            CN[i] = bf3 + i * maxD2;
        }
        CN[0][0] = 1;
        CN[1][0] = 1;
        CN[1][1] = 1;
        for(int i = 2; i < maxD; i++) {
            CN[i][0] = 1;
            if(i < maxD2) CN[i][i] = 1;
            for(int j = 1; j < i && j < maxD2; j++) {
                CN[i][j] = CN[i - 1][j - 1] + CN[i - 1][j];
            }
        }
    }

public:
    kclistClique(Graph && g, ui k) :g(g), k(k) {
        
        maxK = k + 1;
        answer = 0;
        

        maxSize = g.coreNumber + 1;

        computeC();


        nodes.resize(maxK);
        for(ui i = 0; i < maxK; i++) {
            nodes[i].resize(maxSize);
        }

        

        clique.resize(g.coreNumber);
        clique.clear();
        
        printf("kCListClique.h\n");
    }
    ~kclistClique() { 
        if(bf3 != nullptr) {
            delete [] bf3; delete [] CN; 
        }
    }

    double run();
    double runPeel();//kClistClique_Peel.cpp
    double runPeelCore();
    double runV2();//build a subgraph for each ego network
    double runV3();// impletented in V2.cpp
    double run_v4();//greedy large degree node
    double run_v5();//remove multiple cliques, gready large degree
    double run_v6();//core-peeling based
    double run_v7();//remove multiple cliques, core-peeling based

    double runCCPath(ull sps, bool useListing=false);
    double runCCPath_V2SavingLC(ull sps, bool useListing=false);
};

#endif