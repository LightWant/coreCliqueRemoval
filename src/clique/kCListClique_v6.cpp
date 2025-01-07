#include "kCListClique.h"
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

// #define PRINT_ARBORICITY

#ifdef PRINT_ARBORICITY
std::vector<std::pair<ui, ui>> edgesRandom;
std::vector<ui> r;
std::vector<ui> Pava, reordered;

#define PRING_AVG_ARBORICITY
double sRaw = 0.0, sCCr = 0.0, cntNodes__ = 0;
#endif

// #define PRINT_RATIO_OF_KCLIQUES
#ifdef PRINT_RATIO_OF_KCLIQUES
double numberOfKcliquesInCoreClique = 0;
#endif

double kclistClique::run_v6() {
    printf("kClistClique_v6.cpp\n");
#ifdef DEBUG
g.print();
#endif

#ifdef PRINT_ARBORICITY
edgesRandom.resize(g.coreNumber*g.coreNumber);
r.resize(g.coreNumber);
Pava.resize(g.coreNumber);
reordered.resize(g.coreNumber);
#endif

    g.initHash();
    printf("initHash\n");

    edAdj.resize(k);
    for(ui i = 0; i < k; i++) edAdj[i].resize(g.n);
    adj.resize(g.n);
    for(ui i = 0; i < g.n; i++) {
        adj[i].resize(g.pIdx[i+1] - g.pIdx2[i]);
        for(ui j = 0; j < g.pIdx[i+1] - g.pIdx2[i]; j++) {
            adj[i][j] = g.pEdge[g.pIdx2[i] + j];
        }
    }

    level.resize(g.n);
    vis.resize(g.n);    

    std::vector<ui> newCliqueNei;

    ListLinearHeap heap(g.coreNumber, g.coreNumber);
    heap.reserve(g.coreNumber, g.coreNumber);
    // std::vector<ui> keys(g.coreNumber);
    // std::vector<ui> heads(g.coreNumber);
    // std::vector<ui> nxts(g.coreNumber);
 
    std::vector<ui> id_s(g.coreNumber+1);
	std::vector<ui> degree(g.coreNumber);
	std::vector<bool> removed(g.coreNumber);

    std::vector<ui> pIdx(g.coreNumber+1);
    std::vector<ui> pEdge(g.coreNumber * g.coreNumber);

#ifdef PRINT_SG
#ifndef PRINT_AVG_SG
printf("PRINT_SG\n");
#endif
#endif
    for(ui u = 0; u < g.n; u++) if(g.pIdx[u+1] - g.pIdx2[u] >= k-1) {
        // nodes[0].clear();
        // for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
        //     nodes[0].push_back(g.pEdge[i]);
        //     // level[g.pEdge[i]] = 1;
        // }
#ifdef DEBUG
printf("    st %u\n", u);
for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
    printf("%u ", g.pEdge[i]);
}
printf("\n");
#endif
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
        
        // for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
        //     for(ui j = i+1; j < g.pIdx[u + 1]; j++) {
        //         if(g.connectHash(g.pEdge[i], g.pEdge[j])) {
        //             pIdx[i-g.pIdx2[u]+1]++;
        //             pIdx[j-g.pIdx2[u]+1]++;
        //         }
        //     }
        // }
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
// printf("a%u ", u);fflush(stdout);
            //ListLinearHeap *heap = new ListLinearHeap(n, new_size-1);
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
        // for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
        //     if(vis[g.pEdge[i]] == 0) {
        //         nxtC.push_back(g.pEdge[i]);
        //         level[g.pEdge[i]] = 1;
        //     }
        // }
        // for(auto u : clique) level[u] = 1;
        for(auto v: nxtC) level[v] = 1;
        // for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
        //     level[g.pEdge[i]] = 1;
        // }
        
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

#ifdef PRINT_ARBORICITY

double a0 = 0, a1 = 0;
for(ui xx = 0; xx < 2; xx++) {


ui t = 10;
ui i = 0;


edgesRandom.clear();

if(xx == 1) {
    for(auto v: nxtC) level[v] = i++;
    for(auto v:nxtC) {
        for(ui i = 0; i < edAdj[1][v]; i++) {
            edgesRandom.push_back({level[v], level[adj[v][i]]});
        }
    }
    n = nxtC.size();
}
else {
    for(ui i = 0; i < n; i++) {
        for(ui j = pIdx[i]; j < pIdx[i+1]; j++) {
            if(pEdge[j] > i) {
                edgesRandom.push_back({i, pEdge[j]});
            }
        }
    }
}

std::random_device rd;
// std::mt19937 gg(rd());
std::mt19937 gg(31);
typedef std::uniform_int_distribution<ui> distr_t;
typedef typename distr_t::param_type param_t;
distr_t D;
if(edgesRandom.size()>1)
for(ui i = edgesRandom.size()-1; i > 0; i--) {
    ui j = D(gg, param_t(0, i));
    std::swap(edgesRandom[i], edgesRandom[j]);
}

for(ui i = 0; i < n; i++) r[i] = 0;
ui ansE = 0, ansN = 2;

while(true) {
    for(ui i = 0; i < t; i++) {
        for(ui j = 0; j < edgesRandom.size(); j++) {
            ui a = edgesRandom[j].first;
            ui b = edgesRandom[j].second;
            if(r[a] < r[b]) r[a]++;
            else r[b]++;
        }
    }

    for(ui i = 0; i < n; i++) reordered[i] = i;

    std::sort(reordered.begin(), reordered.begin()+n, 
        [&](ui a, ui b) { return r[a] > r[b]; }
    );
    for(ui i = 0; i < n; i++) level[reordered[i]] = i;
    for(ui i = 0; i < n; i++) Pava[i] = 0;
    for(ui j = 0; j < edgesRandom.size(); j++) {
        ui a = edgesRandom[j].first;
        ui b = edgesRandom[j].second;
        if(level[a] > level[b]) Pava[level[a]]++;
        else Pava[level[b]]++;
    }
    ui sum = 0;
    ui sizeDensestSubgraph=0;
    ui eDensestSubgraph=0;
    double rho = 0.0;
    for (ui i = 0; i < n; ++i) {
        sum += Pava[i];
        if ((double)sum / (i + 1) + 1e-5 >= rho) {
            rho = (double)sum / (i + 1);
            eDensestSubgraph = sum;
            sizeDensestSubgraph = i + 1;
        }
    }
// if(n==42) {
// printf("level:"); for(ui i = 0; i < n; i++) printf("%u ", level[i]);printf("\n");
// printf("pava:"); for(ui i = 0; i < n; i++) printf("%u ", Pava[i]);printf("\n");
// }
    if(eDensestSubgraph==ansE && sizeDensestSubgraph == ansN) {
        break;
    }
    else if(rho > 1.0*ansE/ansN + 1e-5) {
        ansE = eDensestSubgraph;
        ansN = sizeDensestSubgraph;
    }
    else if(eDensestSubgraph == 0 && t > 10) {
        break;
    }

    t *= 2;
    if(t > 1000) {
        // printf("%f %f\n", rho, 1.0*ansE/ansN);fflush(stdout);
        break;
    }
}
// if(n==42) {
// printf("%u %u %u %u\n", n, ui(std::ceil(1.0*ansE/(ansN-1))), ansE, ansN);
// printf("sg:\n", u);
// for(ui i = 0; i < n; i++) {
//     printf("%u-%u: ", i, pIdx[i+1]-pIdx[i]);
//     for(ui j = pIdx[i]; j < pIdx[i+1]; j++) {
//         printf("%u ", pEdge[j]);
//     }
//     printf("\n");
// }
// printf("\n");
// }
if(xx == 0) a0 = ui(std::ceil(1.0*ansE/(ansN-1)));
else a1 = ui(std::ceil(1.0*ansE/(ansN-1)));
}
#ifndef PRING_AVG_ARBORICITY
printf("%u %u\n", n, ui(std::ceil(1.0*ansE/(ansN-1))));
#endif

if(a0 > 0) {
sRaw += a1/a0;
cntNodes__+=1;
}
// fflush(stdout);
for(ui i = 0; i < n; i++) level[i] = 0;
for(auto u : clique) vis[u] = 0;
for(ui v : nxtC) level[v] = 0;
continue;
#endif

#ifdef PRINT_SG
double densityAll = 1.0 * pIdx[n] / n;
double densityNxt = 0;
for(auto v:nxtC) densityNxt += edAdj[1][v];
if(nxtC.size() > 0) densityNxt /= nxtC.size();

#ifndef PRINT_AVG_SG
printf("%u %u %f %u %f\n", 
    u, g.pIdx[u+1]-g.pIdx2[u], densityAll, nxtC.size(), densityNxt);
#endif
if(densityAll > 0)
    sRaw += densityNxt/densityAll, cntNodes__ += 1;
#else

#ifdef PRINT_RATIO_OF_KCLIQUES
numberOfKcliquesInCoreClique += CN[clique.size()][k-1];
#else
        listing(1, clique.size());
        // listingV2(1, clique.size());
        // answer += CN[clique.size()][k - 1];
#endif
#endif

#ifdef DEBUG
printf("ans %.0f\n", answer);
#endif
        for(auto u : clique) vis[u] = 0;
        for(ui v : nxtC) level[v] = 0;
        // for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
        //     level[g.pEdge[i]] = 0;
        // }
    }
#ifdef PRINT_SG
printf("avgDensityRaw:%f\n", sRaw / cntNodes__);
// printf("avgDensityCCR:%f\n", sCCr / cntNodes__);
#endif
#ifdef PRINT_ARBORICITY
printf("avgArboricityRaw:%f\n", sRaw / cntNodes__);
#endif

#ifdef PRINT_RATIO_OF_KCLIQUES
printf("numberOfKcliquesInCoreClique:%.0f\n", numberOfKcliquesInCoreClique);
#endif
    return answer;
}


void kclistClique::listing(ui deep, ui edClique) {
#ifdef DDEBUG
printf("    deep %u\n", deep);
printf("C:");
for(auto v: nodes[deep]) printf("%u ", v);printf("\n");
printf("clique:");
for(ui i = 0; i < edClique; i++) printf("%u ", clique[i]);printf("\n");
#endif

    std::vector<ui> & C = nodes[deep];
    if(deep == k-1) {
        // for(auto v:C) answer += edAdj[deep][v];;
        answer += C.size();
        answer += edClique;
// printf("+:%u %u\n", C.size(), edClique);
        return;
    }

    if(edClique >= k - deep) {
        answer += CN[edClique][k - deep];
// printf("+:%.0f\n", CN[edClique][k - deep]);
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

        listing(deep + 1, newEdClique);

        for(auto v : nxtC) level[v] = deep;
    }
}
