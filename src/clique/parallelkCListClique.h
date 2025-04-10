#ifndef PARALLELKCLISTCLIQUE_H
#define PARALLELKCLISTCLIQUE_H

#include "../graph/graph.h"
#include "../tools/types.hpp"
#include <omp.h>

class parallelKclistClique {
private:
    Graph g;
    ui k = 0;
    int nThreads = 1;
    double answer = 0;
    std::vector<std::vector<std::vector<ui>>> nodess;
    std::vector<std::vector<std::vector<ui>>> adjs;
    std::vector<std::vector<std::vector<ui>>> edAdjs;
    std::vector<std::vector<ui>> levels;
    std::vector<std::vector<ui>> cliques;
    std::vector<ui> cliqueNei;
    std::vector<ui> deg;

    void buildSg(ui u);
    void listingParallel(ui deep, ui edClique, double * ans, std::vector<std::vector<ui>> & nodes,
        std::vector<std::vector<ui>> & edAdj,
        std::vector<ui> & level,
        std::vector<std::vector<ui>> & adj,
        std::vector<ui> & clique);

    std::vector<std::vector<ui>> viss;
    

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
    parallelKclistClique(Graph && g, ui k, int nThreads) :g(g), k(k), nThreads(nThreads) {
        maxK = k + 1;
        answer = 0;
        

        maxSize = g.coreNumber + 1;

        computeC();


        nodess.resize(nThreads);
        #pragma omp parallel for
        for(ui ii = 0; ii < nThreads; ii++) {
            nodess[ii].resize(maxK);
            for(ui i = 0; i < maxK; i++) {
                nodess[ii][i].resize(maxSize);
            }
        }
        
        levels.resize(nThreads);
        #pragma omp parallel for
        for(ui i = 0; i < nThreads; i++) levels[i].resize(g.n);
        

        cliques.resize(nThreads);
        #pragma omp parallel for
        for(ui i = 0; i < nThreads; i++) {
            cliques[i].resize(g.coreNumber);
            cliques[i].clear();
        }

        viss.resize(nThreads);
        #pragma omp parallel for
        for(ui i = 0; i < nThreads; i++) {
            viss[i].resize(g.n);
            viss[i].clear();
        }

        edAdjs.resize(nThreads);
        adjs.resize(nThreads);
        #pragma omp parallel for
        for(ui ii = 0; ii < nThreads; ii++) {
            edAdjs[ii].resize(k);
            for(ui i = 0; i < k; i++) edAdjs[ii][i].resize(g.n);
            adjs[ii].resize(g.n);
            for(ui i = 0; i < g.n; i++) {
                adjs[ii][i].resize(g.pIdx[i+1] - g.pIdx2[i]);
                for(ui j = 0; j < g.pIdx[i+1] - g.pIdx2[i]; j++) {
                    adjs[ii][i][j] = g.pEdge[g.pIdx2[i] + j];
                }
            }
        }

        
        
        
        printf("parallelKCListClique.h\n");
    }
    ~parallelKclistClique() { 
        if(bf3 != nullptr) {
            delete [] bf3; delete [] CN; 
        }
    }

    void usingHash() {
        g.initHash(); printf("initHash\n"); 
    }

    double run_v1();//core-peeling based

    double run_v3();

    double run_v2(); //edge parallel core-peeling based
};

#endif