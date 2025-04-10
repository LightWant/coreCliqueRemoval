#pragma once

#include "../graph/graph.h"
#include "../tools/types.hpp"
#include "../tools/linearSet.hpp"
#include <omp.h>

class parallelMCE {
private:
    Graph g;
    ui cnt = 0;
    ui threshold = 0;
    int nThreads = 1;

    std::vector<LinearSet> Vs;
    std::vector<std::vector<std::vector<ui>>> adjs;
    std::vector<std::vector<std::vector<ui>>> edAdjs, stAdjs;
    std::vector<std::vector<ui>> reIds;
    std::vector<std::vector<bool>> inCs;

    // std::vector<std::vector<bool>> viss;

    std::vector<LinearSet> Cv3s, Rv3s, Pv3s, Xv3s;
    std::vector<ui> numOfNonCNodess;
    void BKTomitaRecCCRV3(ui deep, ui ER, ui EC, ui EP, 
        ui SX, ui EX, ui * ans, int threadId);


public:
    parallelMCE(Graph && g, int nThreads) : g(g), nThreads(nThreads) {
        Vs.resize(nThreads);
        for(int i = 0; i < nThreads; i++) {
            Vs[i].resize(g.n);
        }

        Cv3s.resize(nThreads);
        Pv3s.resize(nThreads);
        Rv3s.resize(nThreads);
        Xv3s.resize(nThreads);
        for(int i = 0; i < nThreads; i++) {
            Cv3s[i].resize(g.maxD);
            Pv3s[i].resize(g.maxD);
            Rv3s[i].resize(g.maxD);
            Xv3s[i].resize(g.maxD);
        }

        edAdjs.resize(nThreads);
        stAdjs.resize(nThreads);
        adjs.resize(nThreads);
        reIds.resize(nThreads);
        inCs.resize(nThreads);
        // viss.resize(nThreads);
        for(int ii = 0; ii < nThreads; ii++) {
            edAdjs[ii].resize(g.coreNumber+5);//search depth
            for(ui i = 0; i < g.coreNumber+5; i++) edAdjs[ii][i].resize(g.maxD);
            stAdjs[ii].resize(g.coreNumber+1);
            for(ui i = 0; i < g.coreNumber; i++) stAdjs[ii][i].resize(g.maxD);
            adjs[ii].resize(g.maxD);
            for(ui i = 0; i < g.maxD; i++) {
                adjs[ii][i].resize(g.coreNumber);
            }
            reIds[ii].resize(g.n, g.n);
            // viss[ii].resize(g.n);

            inCs[ii].resize(g.maxD, false);

            
        }

        numOfNonCNodess.resize(nThreads);
        
    }

    ui runCoreCliqueRemovalV3();

    ui runCoreCliqueRemovalV3EdgeParallel();

    ui runCoreCliqueRemovalV3WorkStealing();

    void useHash() {
        g.initHash();
    }
};