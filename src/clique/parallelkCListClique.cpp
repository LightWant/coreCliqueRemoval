#include "parallelkCListClique.h"
#include "../tools/listLinearHeap.hpp"
#include <algorithm>
#include <cassert>
#include <random>
// #define DEBUG
// #define DDEBUG

// #define PRINT_SG

#ifdef PRINT_SG
#define PRINT_AVG_SG
double sRaw = 0.0, sCCr = 0.0, cntNodes__ = 0;
#endif


double parallelKclistClique::run_v3() {
    printf("parallelKclistClique.cpp::run_v3 work stealing\n");

    std::vector<std::vector<ui>> id_ss(nThreads);
    for(ui i = 0; i < nThreads; i++) id_ss[i].resize(g.coreNumber+1);
	std::vector<std::vector<ui>> degrees(nThreads);
    for(ui i = 0; i < nThreads; i++) degrees[i].resize(g.coreNumber);
	std::vector<std::vector<bool>> removeds(nThreads);
    for(ui i = 0; i < nThreads; i++) removeds[i].resize(g.coreNumber);

    std::vector<std::vector<ui>> pIdxs(nThreads);
    for(ui i = 0; i < nThreads; i++) pIdxs[i].resize(g.coreNumber+1);
    std::vector<std::vector<ui>> pEdges(nThreads);
    for(ui i = 0; i < nThreads; i++) pEdges[i].resize(g.coreNumber * g.coreNumber);



    std::vector<ui> pTs(nThreads+1), pTe(nThreads);
    pTs[0] = 0;

    std::vector<ui> weights(g.n);
    #pragma omp for schedule(dynamic, 1) nowait
    for(ui u = 0; u < g.n; u++) {
        if(g.pIdx[u+1] - g.pIdx2[u] < k-1) continue;
        // for(ui j = g.pIdx2[u]; j < g.pIdx[u+1]; j++) {
        //     for(ui l = j + 1; l < g.pIdx[u+1]; l++) {
        //         if(g.connectHash(g.pEdge[j], g.pEdge[l])) {
        //             weights[u]++;
        //         }
        //     }
        // }
        weights[u] = g.pIdx[u+1] - g.pIdx2[u];
    }
    ui chunk = 0;
    for(ui u = 0; u < g.n; u++) chunk += weights[u];
    chunk /=  nThreads;
    ui tmp = 0, p = 0;
    pTs[p] = 0;
    for(ui u = 0; u < g.n; u++) {
        tmp += weights[u];
        if(tmp >= chunk) {
            pTe[p] = u;
            pTs[p+1] = u;
            p++;
            tmp = 0;
        }
    }
    pTe[nThreads-1] = g.n;

    // for(int i = 1; i <= nThreads; i++) {
    //     pTs[i] = pTs[i-1] + g.n / nThreads;
    // }
    // for(int i = 0; i < nThreads; i++) {
    //     pTe[i] = pTs[i+1];
    // }

    // std::vector<ui> coreCliques(g.m);
    // std::vector<ui> pccIdx(g.n+1);
// double timeOfThreads[200];

    
    #pragma omp parallel reduction(+:answer)
    {
        int threadId = omp_get_thread_num();
        std::vector<ui> & level = levels[threadId];
        std::vector<ui> & id_s = id_ss[threadId];
        std::vector<ui> & degree = degrees[threadId];
        std::vector<bool> & removed = removeds[threadId];
        std::vector<ui> & pIdx = pIdxs[threadId];
        std::vector<ui> & pEdge = pEdges[threadId];
        std::vector<ui> & clique = cliques[threadId];
        std::vector<std::vector<ui>> & nodes = nodess[threadId];
        std::vector<ui> & vis = viss[threadId];
        std::vector<std::vector<ui>> & adj = adjs[threadId];
        std::vector<std::vector<ui>> & edAdj = edAdjs[threadId];

        ListLinearHeap heap(g.coreNumber, g.coreNumber);
        heap.reserve(g.coreNumber, g.coreNumber);
        std::vector<ui> newCliqueNei;
        double ansLocal = 0.0;


        // double tS = omp_get_wtime();


        // #pragma omp for schedule(static, 32)
        auto computeU = [&](ui u) {
            #ifdef DEBUG
printf("    st %u\n", u);
for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
    printf("%u ", g.pEdge[i]);
}
printf("\n");
#endif
// printf("    st %u, t %d\n", u, threadId);fflush(stdout);
            for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
                level[g.pEdge[i]] = i-g.pIdx2[u]+1;
    #ifdef DDEBUG
    if(u == 4035354) {
    printf("level %u %u, %u\n", g.pEdge[i], i-g.pIdx2[u]+1, level[g.pEdge[i]]);
    }
    #endif
            }
            ui n = g.pIdx[u+1] - g.pIdx2[u];
            for(ui i = 0; i <= n; i++) pIdx[i] = 0;
            
            for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
                ui v = g.pEdge[i];
    #ifdef DDEBUG
    if(u == 4035354 ) {
    printf("st %u\n", i-g.pIdx2[u]);
    }
    #endif
                for(ui j = g.pIdx2[v]; j < g.pIdx[v+1]; j++) {
                    ui w = g.pEdge[j];
                    if(level[w] > 0) {
    #ifdef DDEBUG
    
    if(u == 4035354 ) {
    
    if(i-g.pIdx2[u]+1 == 1 || level[w]==1) 
    printf("%u(%u)-%u, %u, %u %u %u\n", 
        i-g.pIdx2[u]+1, level[v], level[w], i-g.pIdx2[u], u,v,w);
    }
    #endif
    
                        pIdx[i-g.pIdx2[u]+1]++;
                        pIdx[level[w]]++;
                    }
                }
            }
                    
#ifdef DDEBUG

if(u == 4035354 ) {
printf("deg:");
for(ui i = 0; i < n; i++) {
    printf("%u-%u: ", i, pIdx[i+1]);
}
printf("\n");
}
#endif
            // #endif
            for(ui i = 1; i <= n; i++) pIdx[i] += pIdx[i-1];
            for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
                ui v = g.pEdge[i];
                for(ui j = g.pIdx2[v]; j  < g.pIdx[v+1]; j++) {
                    ui w = g.pEdge[j];
                    if(level[w]) {
                        pEdge[pIdx[i-g.pIdx2[u]]++] = level[w]-1;
                        pEdge[pIdx[level[w]-1]++] = i-g.pIdx2[u];
                    }
                }
            }
            for(ui i = n; i > 0; i--) pIdx[i] = pIdx[i-1];
            pIdx[0] = 0;
#ifdef DDEBUG
if(u == 4035354 ) {
printf("sg:\n", u);
for(ui i = 0; i < n; i++) {
    printf("%u-%u/%u-%u: ", i, pIdx[i+1]-pIdx[i], pIdx[i+1], pIdx[i]);
    for(ui j = pIdx[i]; j < pIdx[i+1]; j++) {
        printf("%u ", pEdge[j]);
    }
    printf("\n");
}
printf("\n");
}
#endif
#ifdef DEBUG
printf("sg:\n", u);
for(ui i = 0; i < n; i++) {
    printf("%u-%u: ", g.pEdge[g.pIdx2[u] + i], pIdx[i+1]-pIdx[i]);
    for(ui j = pIdx[i]; j < pIdx[i+1]; j++) {
        printf("%u ", g.pEdge[g.pIdx2[u] + pEdge[j]]);
    }
    printf("\n");
}
printf("\n");
#endif
            clique.clear();
            std::vector<ui> & nxtC = nodes[1];
            nxtC.clear();
            
            ui queue_n = 0, new_size = 0, threshold = k-2;
            for(ui i = 0;i < n;i ++) removed[i] = false;
            for(ui i = 0;i < n;i ++) {
                if(pIdx[i+1]-pIdx[i] < threshold) id_s[queue_n ++] = i;
                degree[i] = pIdx[i+1]-pIdx[i];
            }
            for(ui i = 0;i < queue_n;i ++) {
                ui uu = id_s[i]; degree[uu] = 0;
                removed[uu] = true;
                for(ui j = pIdx[uu];j < pIdx[uu+1];j ++) if(degree[pEdge[j]] > 0) {
                    if((degree[pEdge[j]] --) == k-2) id_s[queue_n ++] = pEdge[j];
                }
            }
                    
            for(ui i = 0;i < n;i ++) {
                if(degree[i] >= threshold) id_s[queue_n + (new_size ++)] = i;
                else {
                    removed[i] = true;
                }
            }
            assert(queue_n + new_size == n);
            // printf("a\n");fflush(stdout);
            if(new_size == 0) {
                for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
                    level[g.pEdge[i]] = 0;
                }
                // continue;
                return;
            }
            else {
                heap.initRaw(new_size, new_size-1, id_s.data()+queue_n, degree.data());
    
    // printf("b");fflush(stdout);
                for(ui i = 0;i < new_size;i ++) {
                    ui v, key;
    
                    heap.pop_min(v, key);
    
                    id_s[queue_n + i] = v;
                    if(key + i + 1 == new_size) {
    
                        ui x_size = i+1;
                        heap.get_ids(id_s.data()+queue_n, x_size);
    
                        assert(x_size == new_size);
        // printf("	last key %u, x_size %u\n", key, x_size);
                        for(ui j = i;j < new_size;j ++) {
                            clique.pb(g.pEdge[g.pIdx2[u]+id_s[queue_n+j]]);
                        }
                        for(ui j = 0;j < i;j ++) {
                            nxtC.pb(g.pEdge[g.pIdx2[u]+id_s[queue_n+j]]);
                        }
                        break;
                    }
                    removed[v] = 1;
        // printf("%u %u\n", u, key);
        // fflush(stdout);
    
                    for(ui j = pIdx[v];j < pIdx[v+1];j ++) if(removed[pEdge[j]] == false) {
                        heap.decrement(pEdge[j], 1);
                    }
                }
    // printf("c\n");fflush(stdout);
            }
    
            for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
                level[g.pEdge[i]] = 0;
            }
            
                    
            for(auto u : clique) vis[u] = 1;
            for(auto v: nxtC) level[v] = 1;
                    
            #ifdef DEBUG
            printf("nxtC: ");for(auto v:nxtC)printf("%u ", v); printf("\n");
            printf("clique: ");for(auto v:clique)printf("%u ", v); printf("\n");
            
            #endif
            
            for(auto v:nxtC) {
                ui ed = adj[v].size();
                for(ui i = 0; i < ed; ) {
                    if(level[adj[v][i]] != 1) std::swap(adj[v][i], adj[v][--ed]);
                    else i++;
                }
                edAdj[1][v] = ed;
            }
            #ifdef DEBUG
            printf("subgraph %u\n", u);
            for(auto v:nxtC) {
                printf("%u:", v);
                for(ui i = 0; i < edAdj[1][v]; i++) printf("%u ", adj[v][i]);
                printf("\n");
            }
            #endif
       
            listingParallel(1, clique.size(), &ansLocal, nodes, edAdj, level, adj, clique);
            
            #ifdef DEBUG
            printf("ans %.0f\n", ansLocal);
            #endif
            for(auto u : clique) vis[u] = 0;
            for(ui v : nxtC) level[v] = 0;
        };

        ui u;
        while((u = __sync_fetch_and_add(&pTs[threadId], 1)) < pTe[threadId]) {
            if(g.pIdx[u+1] - g.pIdx2[u] < k-1) continue;
            computeU(u);
        }

        for(ui jjj = 0; jjj < nThreads; jjj++) {
            if(jjj == threadId) continue;
            if(pTs[jjj] >= pTe[jjj]) continue;
            while((u = __sync_fetch_and_add(&pTs[jjj], 1)) < pTe[jjj]) {
                if(g.pIdx[u+1] - g.pIdx2[u] < k-1) continue;
                computeU(u);
            }
        }

        // #pragma omp for schedule(dynamic, 1) nowait
        // for(ui u = 0; u < g.n; u++) {
        //     if(g.pIdx[u+1] - g.pIdx2[u] < k-1) continue;
        //     computeU(u);
        // }
        // double tE = omp_get_wtime();
        // timeOfThreads[threadId] = tE-tS;
        
        answer += ansLocal;
        
        // std::cout << threadId << ':' << tE-tS << std::endl;
    }

    // for(int i = 0; i < nThreads; i++) {
    //     std::cout << i << ':' << timeOfThreads[i] << std::endl;
    // }

    return answer;
}

double parallelKclistClique::run_v2() {
    printf("parallelKclistClique.cpp::run_v2\n");
#ifdef DEBUG
g.print();
#endif

    // usingHash();
 
    std::vector<std::vector<ui>> id_ss(nThreads);
    for(ui i = 0; i < nThreads; i++) id_ss[i].resize(g.coreNumber+1);
	std::vector<std::vector<ui>> degrees(nThreads);
    for(ui i = 0; i < nThreads; i++) degrees[i].resize(g.coreNumber);
	std::vector<std::vector<bool>> removeds(nThreads);
    for(ui i = 0; i < nThreads; i++) removeds[i].resize(g.coreNumber);

    std::vector<std::vector<ui>> pIdxs(nThreads);
    for(ui i = 0; i < nThreads; i++) pIdxs[i].resize(g.coreNumber+1);
    std::vector<std::vector<ui>> pEdges(nThreads);
    for(ui i = 0; i < nThreads; i++) pEdges[i].resize(g.coreNumber * g.coreNumber);

    // std::vector<ui> coreCliques(g.m);
    // std::vector<ui> pccIdx(g.n+1);
// double timeOfThreads[200];
    std::vector<ui> coreCliquesBuffer(g.m/2);
    std::vector<std::pair<ui, ui>> nxtCandidates(g.m/2);

    std::vector<ui> pSt(g.n), pClique(g.n+1), pCandidate(g.n+1);
    for(ui u = 0; u < g.n; u++) {
        pClique[u+1] = g.pIdx[u+1]-g.pIdx2[u];
    }
    for(ui u = 0; u < g.n; u++) {
        pClique[u+1] += pClique[u];
    }
    for(ui u = 0; u < g.n; u++) {
        pCandidate[u] = pSt[u] = pClique[u];
    }
    ui sumOfEdgeCandidates = 0;
    
    
    #pragma omp parallel reduction(+:answer) shared(sumOfEdgeCandidates)
    {
        int threadId = omp_get_thread_num();
        std::vector<ui> & level = levels[threadId];
        std::vector<ui> & id_s = id_ss[threadId];
        std::vector<ui> & degree = degrees[threadId];
        std::vector<bool> & removed = removeds[threadId];
        std::vector<ui> & pIdx = pIdxs[threadId];
        std::vector<ui> & pEdge = pEdges[threadId];
        std::vector<ui> & clique = cliques[threadId];
        std::vector<std::vector<ui>> & nodes = nodess[threadId];
        std::vector<ui> & vis = viss[threadId];
        std::vector<std::vector<ui>> & adj = adjs[threadId];
        std::vector<std::vector<ui>> & edAdj = edAdjs[threadId];

        ListLinearHeap heap(g.coreNumber, g.coreNumber);
        heap.reserve(g.coreNumber, g.coreNumber);
        std::vector<ui> newCliqueNei;
        double ansLocal = 0.0;


        // double tS = omp_get_wtime();


        // #pragma omp for schedule(static, 32)
        #pragma omp for 
        for(ui u = 0; u < g.n; u++) if(g.pIdx[u+1] - g.pIdx2[u] >= k-1) {
#ifdef DEBUG
printf("    st %u\n", u);
for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
    printf("%u ", g.pEdge[i]);
}
printf("\n");
#endif
// printf("    st %u, t %d\n", u, threadId);fflush(stdout);
            for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
                level[g.pEdge[i]] = i-g.pIdx2[u]+1;
    #ifdef DDEBUG
    if(u == 4035354) {
    printf("level %u %u, %u\n", g.pEdge[i], i-g.pIdx2[u]+1, level[g.pEdge[i]]);
    }
    #endif
            }
            ui n = g.pIdx[u+1] - g.pIdx2[u];
            for(ui i = 0; i <= n; i++) pIdx[i] = 0;
            
            for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
                ui v = g.pEdge[i];
    #ifdef DDEBUG
    if(u == 4035354 ) {
    printf("st %u\n", i-g.pIdx2[u]);
    }
    #endif
                for(ui j = g.pIdx2[v]; j < g.pIdx[v+1]; j++) {
                    ui w = g.pEdge[j];
                    if(level[w] > 0) {
    #ifdef DDEBUG
    
    if(u == 4035354 ) {
    
    if(i-g.pIdx2[u]+1 == 1 || level[w]==1) 
    printf("%u(%u)-%u, %u, %u %u %u\n", 
        i-g.pIdx2[u]+1, level[v], level[w], i-g.pIdx2[u], u,v,w);
    }
    #endif
    
                        pIdx[i-g.pIdx2[u]+1]++;
                        pIdx[level[w]]++;
                    }
                }
            }
                    
#ifdef DDEBUG

if(u == 4035354 ) {
printf("deg:");
for(ui i = 0; i < n; i++) {
    printf("%u-%u: ", i, pIdx[i+1]);
}
printf("\n");
}
#endif
            // #endif
            for(ui i = 1; i <= n; i++) pIdx[i] += pIdx[i-1];
            for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
                ui v = g.pEdge[i];
                for(ui j = g.pIdx2[v]; j  < g.pIdx[v+1]; j++) {
                    ui w = g.pEdge[j];
                    if(level[w]) {
                        pEdge[pIdx[i-g.pIdx2[u]]++] = level[w]-1;
                        pEdge[pIdx[level[w]-1]++] = i-g.pIdx2[u];
                    }
                }
            }
            for(ui i = n; i > 0; i--) pIdx[i] = pIdx[i-1];
            pIdx[0] = 0;
#ifdef DDEBUG
if(u == 4035354 ) {
printf("sg:\n", u);
for(ui i = 0; i < n; i++) {
    printf("%u-%u/%u-%u: ", i, pIdx[i+1]-pIdx[i], pIdx[i+1], pIdx[i]);
    for(ui j = pIdx[i]; j < pIdx[i+1]; j++) {
        printf("%u ", pEdge[j]);
    }
    printf("\n");
}
printf("\n");
}
#endif
#ifdef DEBUG
printf("sg:\n", u);
for(ui i = 0; i < n; i++) {
    printf("%u-%u: ", g.pEdge[g.pIdx2[u] + i], pIdx[i+1]-pIdx[i]);
    for(ui j = pIdx[i]; j < pIdx[i+1]; j++) {
        printf("%u ", g.pEdge[g.pIdx2[u] + pEdge[j]]);
    }
    printf("\n");
}
printf("\n");
#endif
            clique.clear();
            std::vector<ui> & nxtC = nodes[1];
            nxtC.clear();
            
            ui queue_n = 0, new_size = 0, threshold = k-2;
            for(ui i = 0;i < n;i ++) removed[i] = false;
            for(ui i = 0;i < n;i ++) {
                if(pIdx[i+1]-pIdx[i] < threshold) id_s[queue_n ++] = i;
                degree[i] = pIdx[i+1]-pIdx[i];
            }
            for(ui i = 0;i < queue_n;i ++) {
                ui uu = id_s[i]; degree[uu] = 0;
                removed[uu] = true;
                for(ui j = pIdx[uu];j < pIdx[uu+1];j ++) if(degree[pEdge[j]] > 0) {
                    if((degree[pEdge[j]] --) == k-2) id_s[queue_n ++] = pEdge[j];
                }
            }
                    
            for(ui i = 0;i < n;i ++) {
                if(degree[i] >= threshold) id_s[queue_n + (new_size ++)] = i;
                else {
                    removed[i] = true;
                }
            }
            assert(queue_n + new_size == n);
            // printf("a\n");fflush(stdout);
            if(new_size == 0) {
                for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
                    level[g.pEdge[i]] = 0;
                }
                continue;
            }
            else {
                heap.initRaw(new_size, new_size-1, id_s.data()+queue_n, degree.data());
    
    // printf("b");fflush(stdout);
                for(ui i = 0;i < new_size;i ++) {
                    ui v, key;
    
                    heap.pop_min(v, key);
    
                    id_s[queue_n + i] = v;
                    if(key + i + 1 == new_size) {
    
                        ui x_size = i+1;
                        heap.get_ids(id_s.data()+queue_n, x_size);
    
                        // assert(x_size == new_size);
        // printf("	last key %u, x_size %u\n", key, x_size);
                        for(ui j = i;j < new_size;j ++) {
                            // clique.pb(g.pEdge[g.pIdx2[u]+id_s[queue_n+j]]);
#ifdef DEBUG
printf("findClique, u %u, v %u, pClique %u\n", 
        u, g.pEdge[g.pIdx2[u]+id_s[queue_n+j]], pClique[u]);fflush(stdout);
#endif
                            coreCliquesBuffer[pClique[u]++] = g.pEdge[g.pIdx2[u]+id_s[queue_n+j]];
                        }
                        // pCandidate[u+1] = i;
                        for(ui j = 0; j < i; j++) {
                            nxtCandidates[pCandidate[u]].first = u;
                            nxtCandidates[pCandidate[u]].second = g.pEdge[g.pIdx2[u]+id_s[queue_n+j]];
                            pCandidate[u]++;
                        }

if(new_size-i >= k - 1) {
#ifdef DEBUG
printf("add (%u,%u) %.0f \n", new_size-i, k-1,  CN[new_size-i][k-1]);
#endif

                        ansLocal += CN[new_size-i][k-1];
}
                        // for(ui j = 0;j < i;j ++) {
                            // nxtC.pb(g.pEdge[g.pIdx2[u]+id_s[queue_n+j]]);
                        // }
                        break;
                    }
                    removed[v] = 1;
        // printf("%u %u\n", u, key);
        // fflush(stdout);
    
                    for(ui j = pIdx[v];j < pIdx[v+1];j ++) if(removed[pEdge[j]] == false) {
                        heap.decrement(pEdge[j], 1);
                    }
                }
    // printf("c\n");fflush(stdout);
            }
    
            for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
                level[g.pEdge[i]] = 0;
            }


        }
        
        // #pragma omp single
        // {
        //     for(ui u = 1; u <= g.n; u++) {
        //         pCandidate[u] += pCandidate[u-1];
        //     }
        //     sumOfEdgeCandidates = pCandidate[g.n];
        // }

        #pragma omp single
        {   
            ui e = 0;
            for(ui u = 0; u < g.n; u++) {
                for(ui i = pSt[u]; i < pCandidate[u]; i++) {
                    nxtCandidates[e++] = nxtCandidates[i]; 
                }
            }
            sumOfEdgeCandidates = e;
        }
        

        // #pragma omp barrier
        // printf("edgeCanNum:%u \n", sumOfEdgeCandidates);
#ifdef DEBUG
printf("coreClique:\n");
for(ui u = 0; u < g.n; u++) {
    printf("%u:", u);
    for(ui i = pSt[u]; i < pClique[u]; i++) {
        printf("%u ", coreCliquesBuffer[i]);
    }
    printf("\n");
}
printf("pCandidte:\n");
for(ui u = 0; u < g.n; u++) {
    printf("%u-%u ", u, pCandidate[u]);
}
printf("\n");
#endif

        // #pragma omp for schedule(dynamic, 1)
        // for(ui u = 0; u < g.n; u++) if(g.pIdx[u+1] - g.pIdx2[u] >= k-1) {

        //     for(ui i = pSt[u]; i < pClique[u]; i++) {
        //         level[coreCliquesBuffer[i]] = 1;
        //     }

        //     for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
        //         if(level[g.pEdge[i]] == 1) continue;
        //         nxtCandidates[pCandidate[u]].first = u;
        //         nxtCandidates[pCandidate[u]].second = g.pEdge[i];
        //         pCandidate[u]++;
        //     }

        //     for(ui i = pSt[u]; i < pClique[u]; i++) {
        //         level[coreCliquesBuffer[i]] = 0;
        //     }
        // }
        
#ifdef DEBUG
printf("canEdegs:\n");
for(ui e = 0; e < sumOfEdgeCandidates; e++) {
    ui u = nxtCandidates[e].first;
    ui v = nxtCandidates[e].second;
    printf("%u-%u ", u, v);
}
printf("\n");
#endif

        #pragma omp for schedule(dynamic, 1) nowait
        for(ui e = 0; e < sumOfEdgeCandidates; e++) {
            ui u = nxtCandidates[e].first;
            ui v = nxtCandidates[e].second;

            clique.clear();
            for(ui i = pSt[u]; i < pClique[u]; i++) {
                ui c = coreCliquesBuffer[i];
                if(g.connectHash(v, c)) {
                    clique.push_back(c);
                }
            }
            for(auto u : clique) vis[u] = 1;

            ui deep = 2;
            std::vector<ui> & nxtC = nodes[deep];
            nxtC.clear();
            for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
                ui w = g.pEdge[i];
                if(w <= v || vis[w]) continue;
                if(g.connectHash(v, w)) {
                    nxtC.push_back(w);
                }
            }
            for(auto v: nxtC) level[v] = 1;

            for(auto v:nxtC) {
                ui ed = adj[v].size();
                for(ui i = 0; i < ed; ) {
                    if(level[adj[v][i]] != 1) std::swap(adj[v][i], adj[v][--ed]);
                    else i++;
                }
                edAdj[deep][v] = ed;
            }

#ifdef DEBUG
printf("run %u-%u:, ansLOCAl %.0f\n", u, v, ansLocal);
for(auto v: clique) printf("%u ", v); printf("\n");
for(auto v: nxtC) printf("%u ", v); printf("\n");
#endif
            if(nxtC.size() == 0) {
                if(clique.size() >= k-deep) ansLocal += CN[clique.size()][k-2];
            }
            else{
                listingParallel(deep, clique.size(), &ansLocal, nodes, edAdj, level, adj, clique);
            }            
#ifdef DEBUG
printf("afterrun:, ansLOCAl %.0f, t %u, %u-%u\n", 
    ansLocal, threadId, u, v);
#endif
            for(auto u : clique) vis[u] = 0;
            for(ui v : nxtC) level[v] = 0;
        }

        

        // double tE = omp_get_wtime();
        // timeOfThreads[threadId] = tE-tS;
        
        answer += ansLocal;
        
        // std::cout << threadId << ':' << tE-tS << std::endl;
    }

    // for(int i = 0; i < nThreads; i++) {
    //     std::cout << i << ':' << timeOfThreads[i] << std::endl;
    // }

    return answer;
}



double parallelKclistClique::run_v1() {
    printf("parallelKclistClique.cpp::run_v1\n");
#ifdef DEBUG
g.print();
#endif

    // usingHash();
 
    std::vector<std::vector<ui>> id_ss(nThreads);
    for(ui i = 0; i < nThreads; i++) id_ss[i].resize(g.coreNumber+1);
	std::vector<std::vector<ui>> degrees(nThreads);
    for(ui i = 0; i < nThreads; i++) degrees[i].resize(g.coreNumber);
	std::vector<std::vector<bool>> removeds(nThreads);
    for(ui i = 0; i < nThreads; i++) removeds[i].resize(g.coreNumber);

    std::vector<std::vector<ui>> pIdxs(nThreads);
    for(ui i = 0; i < nThreads; i++) pIdxs[i].resize(g.coreNumber+1);
    std::vector<std::vector<ui>> pEdges(nThreads);
    for(ui i = 0; i < nThreads; i++) pEdges[i].resize(g.coreNumber * g.coreNumber);

    // std::vector<ui> coreCliques(g.m);
    // std::vector<ui> pccIdx(g.n+1);
// double timeOfThreads[200];

    
    #pragma omp parallel reduction(+:answer)
    {
        int threadId = omp_get_thread_num();
        std::vector<ui> & level = levels[threadId];
        std::vector<ui> & id_s = id_ss[threadId];
        std::vector<ui> & degree = degrees[threadId];
        std::vector<bool> & removed = removeds[threadId];
        std::vector<ui> & pIdx = pIdxs[threadId];
        std::vector<ui> & pEdge = pEdges[threadId];
        std::vector<ui> & clique = cliques[threadId];
        std::vector<std::vector<ui>> & nodes = nodess[threadId];
        std::vector<ui> & vis = viss[threadId];
        std::vector<std::vector<ui>> & adj = adjs[threadId];
        std::vector<std::vector<ui>> & edAdj = edAdjs[threadId];

        ListLinearHeap heap(g.coreNumber, g.coreNumber);
        heap.reserve(g.coreNumber, g.coreNumber);
        std::vector<ui> newCliqueNei;
        double ansLocal = 0.0;


        // double tS = omp_get_wtime();


        // #pragma omp for schedule(static, 32)
        #pragma omp for schedule(dynamic, 1) nowait
        for(ui u = 0; u < g.n; u++) if(g.pIdx[u+1] - g.pIdx2[u] >= k-1) {
#ifdef DEBUG
printf("    st %u\n", u);
for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
    printf("%u ", g.pEdge[i]);
}
printf("\n");
#endif
// printf("    st %u, t %d\n", u, threadId);fflush(stdout);
            for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
                level[g.pEdge[i]] = i-g.pIdx2[u]+1;
    #ifdef DDEBUG
    if(u == 4035354) {
    printf("level %u %u, %u\n", g.pEdge[i], i-g.pIdx2[u]+1, level[g.pEdge[i]]);
    }
    #endif
            }
            ui n = g.pIdx[u+1] - g.pIdx2[u];
            for(ui i = 0; i <= n; i++) pIdx[i] = 0;
            
            for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
                ui v = g.pEdge[i];
    #ifdef DDEBUG
    if(u == 4035354 ) {
    printf("st %u\n", i-g.pIdx2[u]);
    }
    #endif
                for(ui j = g.pIdx2[v]; j < g.pIdx[v+1]; j++) {
                    ui w = g.pEdge[j];
                    if(level[w] > 0) {
    #ifdef DDEBUG
    
    if(u == 4035354 ) {
    
    if(i-g.pIdx2[u]+1 == 1 || level[w]==1) 
    printf("%u(%u)-%u, %u, %u %u %u\n", 
        i-g.pIdx2[u]+1, level[v], level[w], i-g.pIdx2[u], u,v,w);
    }
    #endif
    
                        pIdx[i-g.pIdx2[u]+1]++;
                        pIdx[level[w]]++;
                    }
                }
            }
                    
#ifdef DDEBUG

if(u == 4035354 ) {
printf("deg:");
for(ui i = 0; i < n; i++) {
    printf("%u-%u: ", i, pIdx[i+1]);
}
printf("\n");
}
#endif
            // #endif
            for(ui i = 1; i <= n; i++) pIdx[i] += pIdx[i-1];
            for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
                ui v = g.pEdge[i];
                for(ui j = g.pIdx2[v]; j  < g.pIdx[v+1]; j++) {
                    ui w = g.pEdge[j];
                    if(level[w]) {
                        pEdge[pIdx[i-g.pIdx2[u]]++] = level[w]-1;
                        pEdge[pIdx[level[w]-1]++] = i-g.pIdx2[u];
                    }
                }
            }
            for(ui i = n; i > 0; i--) pIdx[i] = pIdx[i-1];
            pIdx[0] = 0;
#ifdef DDEBUG
if(u == 4035354 ) {
printf("sg:\n", u);
for(ui i = 0; i < n; i++) {
    printf("%u-%u/%u-%u: ", i, pIdx[i+1]-pIdx[i], pIdx[i+1], pIdx[i]);
    for(ui j = pIdx[i]; j < pIdx[i+1]; j++) {
        printf("%u ", pEdge[j]);
    }
    printf("\n");
}
printf("\n");
}
#endif
#ifdef DEBUG
printf("sg:\n", u);
for(ui i = 0; i < n; i++) {
    printf("%u-%u: ", g.pEdge[g.pIdx2[u] + i], pIdx[i+1]-pIdx[i]);
    for(ui j = pIdx[i]; j < pIdx[i+1]; j++) {
        printf("%u ", g.pEdge[g.pIdx2[u] + pEdge[j]]);
    }
    printf("\n");
}
printf("\n");
#endif
            clique.clear();
            std::vector<ui> & nxtC = nodes[1];
            nxtC.clear();
            
            ui queue_n = 0, new_size = 0, threshold = k-2;
            for(ui i = 0;i < n;i ++) removed[i] = false;
            for(ui i = 0;i < n;i ++) {
                if(pIdx[i+1]-pIdx[i] < threshold) id_s[queue_n ++] = i;
                degree[i] = pIdx[i+1]-pIdx[i];
            }
            for(ui i = 0;i < queue_n;i ++) {
                ui uu = id_s[i]; degree[uu] = 0;
                removed[uu] = true;
                for(ui j = pIdx[uu];j < pIdx[uu+1];j ++) if(degree[pEdge[j]] > 0) {
                    if((degree[pEdge[j]] --) == k-2) id_s[queue_n ++] = pEdge[j];
                }
            }
                    
            for(ui i = 0;i < n;i ++) {
                if(degree[i] >= threshold) id_s[queue_n + (new_size ++)] = i;
                else {
                    removed[i] = true;
                }
            }
            assert(queue_n + new_size == n);
            // printf("a\n");fflush(stdout);
            if(new_size == 0) {
                for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
                    level[g.pEdge[i]] = 0;
                }
                continue;
            }
            else {
                heap.initRaw(new_size, new_size-1, id_s.data()+queue_n, degree.data());
    
    // printf("b");fflush(stdout);
                for(ui i = 0;i < new_size;i ++) {
                    ui v, key;
    
                    heap.pop_min(v, key);
    
                    id_s[queue_n + i] = v;
                    if(key + i + 1 == new_size) {
    
                        ui x_size = i+1;
                        heap.get_ids(id_s.data()+queue_n, x_size);
    
                        assert(x_size == new_size);
        // printf("	last key %u, x_size %u\n", key, x_size);
                        for(ui j = i;j < new_size;j ++) {
                            clique.pb(g.pEdge[g.pIdx2[u]+id_s[queue_n+j]]);
                        }
                        for(ui j = 0;j < i;j ++) {
                            nxtC.pb(g.pEdge[g.pIdx2[u]+id_s[queue_n+j]]);
                        }
                        break;
                    }
                    removed[v] = 1;
        // printf("%u %u\n", u, key);
        // fflush(stdout);
    
                    for(ui j = pIdx[v];j < pIdx[v+1];j ++) if(removed[pEdge[j]] == false) {
                        heap.decrement(pEdge[j], 1);
                    }
                }
    // printf("c\n");fflush(stdout);
            }
    
            for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
                level[g.pEdge[i]] = 0;
            }
            
                    
            for(auto u : clique) vis[u] = 1;
            for(auto v: nxtC) level[v] = 1;
                    
            #ifdef DEBUG
            printf("nxtC: ");for(auto v:nxtC)printf("%u ", v); printf("\n");
            printf("clique: ");for(auto v:clique)printf("%u ", v); printf("\n");
            
            #endif
            
            for(auto v:nxtC) {
                ui ed = adj[v].size();
                for(ui i = 0; i < ed; ) {
                    if(level[adj[v][i]] != 1) std::swap(adj[v][i], adj[v][--ed]);
                    else i++;
                }
                edAdj[1][v] = ed;
            }
            #ifdef DEBUG
            printf("subgraph %u\n", u);
            for(auto v:nxtC) {
                printf("%u:", v);
                for(ui i = 0; i < edAdj[1][v]; i++) printf("%u ", adj[v][i]);
                printf("\n");
            }
            #endif
       
            listingParallel(1, clique.size(), &ansLocal, nodes, edAdj, level, adj, clique);
            
            #ifdef DEBUG
            printf("ans %.0f\n", ansLocal);
            #endif
            for(auto u : clique) vis[u] = 0;
            for(ui v : nxtC) level[v] = 0;
        }
        // double tE = omp_get_wtime();
        // timeOfThreads[threadId] = tE-tS;
        
        answer += ansLocal;
        
        // std::cout << threadId << ':' << tE-tS << std::endl;
    }

    // for(int i = 0; i < nThreads; i++) {
    //     std::cout << i << ':' << timeOfThreads[i] << std::endl;
    // }

    return answer;
}


void parallelKclistClique::listingParallel(
    ui deep, ui edClique, double * ans,
    std::vector<std::vector<ui>> & nodes,
    std::vector<std::vector<ui>> & edAdj,
    std::vector<ui> & level,
    std::vector<std::vector<ui>> & adj,
    std::vector<ui> & clique

) {
// printf("\nclique:");
// for(ui i = 0; i < edClique; i++) 
//     printf("%u ", clique[i]); 
// printf("\n");
// printf("nodes[%u]:", deep);
// for(ui i = 0; i < nodes[deep].size(); i++) {
//     ui u = nodes[deep][i];
//     printf("%u: ", u); 
//     for(ui j = 0; j < edAdj[deep][u]; j++) {
//         ui v = adj[u][j]; printf("%u ", v);
//     }
//     printf("\n");
// }
// printf("\n");

    std::vector<ui> & C = nodes[deep];
    if(deep == k-1) {
        (*ans) += C.size();
        (*ans) += edClique;
        return;
    }

    if(edClique >= k - deep) {
        // for(ui i = 0; i < CN[edClique][k - deep]; i++) {
        //     answer++;
        // }
        (*ans) += CN[edClique][k - deep];
    }
        

    std::vector<ui> & nxtC = nodes[deep + 1];
    
    for(ui i = 0; i < C.size(); i++) {
        ui u = C[i];
        if(edAdj[deep][u]+deep+edClique+1 < k) continue;

        nxtC.clear();
        for(ui j = 0; j < edAdj[deep][u]; j++) {
            ui v = adj[u][j];
            nxtC.push_back(v);
            level[v] = deep + 1;
        }

        for(ui j = 0; j < nxtC.size(); j++) {
            ui v = nxtC[j];
            ui & ed = edAdj[deep + 1][v];
            ed = edAdj[deep][v];

            for(ui l = 0; l < ed; ) {
                ui w = adj[v][l];
                if(level[w] == deep + 1) l++;
                else std::swap(adj[v][l], adj[v][--ed]);
            }
        }

        ui newEdClique = edClique;
        for(ui j = 0; j < newEdClique; ) {
            if(!g.connectHash(u, clique[j])) {
                std::swap(clique[j], clique[--newEdClique]);
            }
            else j++;
        }

        listingParallel(deep + 1, newEdClique, ans,  nodes, edAdj, level, adj, clique);

        for(auto v : nxtC) level[v] = deep;
    }
}