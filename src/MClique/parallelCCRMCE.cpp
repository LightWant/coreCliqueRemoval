#include "parallelCCRMCE.h"
#include "../tools/listLinearHeap.hpp"

// #define DEBUG
#ifdef DEBUG
ui xxx=10, xx=10;
#endif

ui parallelMCE::runCoreCliqueRemovalV3WorkStealing() {
    printf("MCE::runCoreCliqueRemovalV3WorkStealing\n");

    std::vector<ui> pTs(nThreads+1), pTe(nThreads);
    pTs[0] = 0;

    std::vector<ui> weights(g.n);
    #pragma omp for schedule(dynamic, 1) nowait
    for(ui u = 0; u < g.n; u++) {
        // if(g.pIdx[u+1] - g.pIdx2[u] < k-1) continue;
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

    #pragma omp parallel reduction(+:cnt)
    {
        double tS = omp_get_wtime();

        int threadId = omp_get_thread_num();
        std::vector<ui> newCliqueNei;

        ListLinearHeap heap(g.coreNumber, g.coreNumber);
        heap.reserve(g.coreNumber, g.coreNumber);

        std::vector<ui> id_s(g.coreNumber+1);
        std::vector<ui> degree(g.coreNumber);
        std::vector<bool> removed(g.coreNumber);

        std::vector<ui> pIdx(g.coreNumber+1);
        std::vector<ui> eIdx(g.coreNumber+1);
        std::vector<ui> pEdge(g.coreNumber * g.coreNumber);

        std::vector<ui> clique(g.coreNumber);
        std::vector<ui> nxtC(g.coreNumber);

        LinearSet & V = Vs[threadId];
        std::vector<std::vector<ui>> & adj = adjs[threadId];
        std::vector<std::vector<ui>> & edAdj = edAdjs[threadId];
        std::vector<std::vector<ui>> & stAdj = stAdjs[threadId];
        std::vector<ui> & reId = reIds[threadId];
        std::vector<bool> & inC = inCs[threadId];

        LinearSet & Cv3 = Cv3s[threadId];
        LinearSet & Rv3 = Rv3s[threadId];
        LinearSet & Pv3 = Pv3s[threadId];
        LinearSet & Xv3 = Xv3s[threadId];

        ui ans = 0;

        // double tS = omp_get_wtime();

        auto computeU = [&](ui u) {
#ifdef DEBUG
xx=u;
if(xx == xxx) 
printf("\n    st %u\n", u);
#endif  
        
            for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
                reId[g.pEdge[i]] = i-g.pIdx2[u]+1;
            }
    
            ui n = g.pIdx[u+1] - g.pIdx2[u];
            for(ui i = 0; i <= n; i++) pIdx[i] = 0;
    
            for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
                ui v = g.pEdge[i];
    
                for(ui j = g.pIdx2[v]; j < g.pIdx[v+1]; j++) {
                    ui w = g.pEdge[j];
                    if(reId[w] < g.n) {
                        pIdx[i-g.pIdx2[u]+1]++;
                        pIdx[reId[w]]++;
                    }
                }
            }
            for(ui i = 1; i <= n; i++) pIdx[i] += pIdx[i-1];
            for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
                ui v = g.pEdge[i];
                for(ui j = g.pIdx2[v]; j  < g.pIdx[v+1]; j++) {
                    ui w = g.pEdge[j];
                    if(reId[w] < g.n) {
                        pEdge[pIdx[i-g.pIdx2[u]]++] = reId[w]-1;
                        pEdge[pIdx[reId[w]-1]++] = i-g.pIdx2[u];
                    }
                }
            }
            for(ui i = n; i > 0; i--) pIdx[i] = pIdx[i-1];
            pIdx[0] = 0;
    
            clique.clear();
            nxtC.clear();
    
            ui queue_n = 0, new_size = 0;
            for(ui i = 0;i < n;i ++) removed[i] = false;
            for(ui i = 0;i < n;i ++) {
                if(pIdx[i+1]-pIdx[i] < threshold) id_s[queue_n ++] = i;
                degree[i] = pIdx[i+1]-pIdx[i];
            }
            for(ui i = 0;i < queue_n;i ++) {
                ui uu = id_s[i]; degree[uu] = 0;
                removed[uu] = true;
                for(ui j = pIdx[uu];j < pIdx[uu+1];j ++) if(degree[pEdge[j]] > 0) {
                    if((degree[pEdge[j]] --) == threshold) id_s[queue_n ++] = pEdge[j];
                }
            }
            
            for(ui i = 0;i < n;i ++) {
                if(degree[i] >= threshold) id_s[queue_n + (new_size ++)] = i;
                else {
                    removed[i] = true;
                }
            }
    
            if(new_size == 0) {
                for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
                    reId[g.pEdge[i]] = g.n;
                }
                // continue;
                return;
            }
            else {
    #ifdef DEBUG
    if(xx == xxx) {
    printf("queue_n:%u\n", queue_n);//should be zero by default
    printf("To peel:");
    for(ui i = 0; i < new_size; i++) {
        printf("%u-%u ", id_s[queue_n + i]+1, 
            degree[id_s[queue_n + i]]);
        
    }
    printf("\n");fflush(stdout);
    }
    #endif
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
                            clique.push_back(g.pEdge[g.pIdx2[u]+id_s[queue_n+j]]);
                        }
                        
                        for(ui j = 0; j < i; j++) {
                            //gready extend
                            ui tu = g.pEdge[g.pIdx2[u]+id_s[queue_n+j]];
                            // bool connectAll = true;
                            // for(ui cu : clique) {
                            //     if(cu < tu) {
                            //         if(!g.connectOut(cu, tu)) {
                            //             connectAll = false; break;
                            //         }
                            //     }
                            //     else if(!g.connectOut(tu, cu)) {
                            //         connectAll = false; break;
                            //     }
                            // }
                            // if(connectAll)  clique.push_back(tu);
                            // else 
                                nxtC.push_back(tu);
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
    
    // printf("findClique\n");fflush(stdout);        
    
            for(ui i = g.pIdx[u]; i < g.pIdx[u + 1]; i++) {
                reId[g.pEdge[i]] = i - g.pIdx[u];
            }
    
    #ifdef DEBUG
    if(xx == xxx) {
    
    printf("nxtC: ");for(auto v:nxtC)printf("%u ", reId[v]); printf("\n");
    printf("clique: ");for(auto v:clique)printf("%u ", reId[v]); printf("\n");
    fflush(stdout);
    }
    #endif
            ui EC = 0;
            ui EX = g.maxD, SX = g.maxD;
            ui ER = 0;
            ui EP = 0;
            for(ui i = g.pIdx[u]; i < g.pIdx2[u]; i++)  {
                Xv3.changeTo(reId[g.pEdge[i]], --SX);
            }
            for(ui v : nxtC)  {
                Pv3.changeTo(reId[v], EP++);
            }
            for(ui v : clique) {
                Cv3.changeTo(reId[v], EC++);
            }
    #ifdef DEBUG
    if(xx == xxx) {
    printf("EC %u, EP %u, EX %u\n", EC, EP, EX);
    fflush(stdout);
    }
    #endif
            //X nodes, the neighbors in P
            for(ui i = g.pIdx[u]; i < g.pIdx2[u]; i++) {
                ui v = g.pEdge[i];
                ui vId = i - g.pIdx[u];
                edAdj[0][vId] = 0;
                for(ui j = g.pIdx2[v]; j < g.pIdx[v+1]; j++) {
                    ui w = g.pEdge[j]; if(reId[w] == g.n) continue;
                    if(Pv3.isIn(reId[w], 0, EP)) {
                        adj[vId][edAdj[0][vId]++] = reId[w];
                    }
                }
            }
            //P/C nodes, the neighbors in P
            for(ui i = g.pIdx2[u]; i < g.pIdx[u+1]; i++) {
                ui vId = i - g.pIdx[u]; edAdj[0][vId] = 0;
            }
            for(ui v : nxtC) {
                ui vId = reId[v];
                for(ui j = g.pIdx2[v]; j < g.pIdx[v+1]; j++) {
                    ui w = g.pEdge[j];
                    ui wId = reId[w];  if(reId[w] == g.n) continue;
                    if(Pv3.isIn(wId, 0, EP)) {
                        adj[vId][edAdj[0][vId]++] = wId;
                        adj[wId][edAdj[0][wId]++] = reId[v];
                    }
                    else if(Cv3.isIn(wId, 0, EC)) {
                        adj[wId][edAdj[0][wId]++] = reId[v];
                    }
                }
            }
            for(ui v : clique) {
                ui vId = reId[v];
                for(ui j = g.pIdx2[v]; j < g.pIdx[v+1]; j++) {
                    ui w = g.pEdge[j];
                    ui wId = reId[w];  if(reId[w] == g.n) continue;
                    if(Pv3.isIn(wId, 0, EP)) {
                        adj[vId][edAdj[0][vId]++] = wId;
                    }
                }
            }
            //X nodes, the neighbors in C
            for(ui i = g.pIdx[u]; i < g.pIdx2[u]; i++) {
                ui v = g.pEdge[i];
                ui vId = i - g.pIdx[u];
                stAdj[0][vId] = g.coreNumber;
    
                for(ui j = g.pIdx2[v]; j < g.pIdx[v+1]; j++) {
                    ui w = g.pEdge[j];  if(reId[w] == g.n) continue;
                    if(Cv3.isIn(reId[w], 0, EC)) {
                        adj[vId][--stAdj[0][vId]] = reId[w];
                    }
                }
            }
            //P/C nodes, the neighbors in C
            for(ui v : nxtC) {
                ui vId = reId[v]; stAdj[0][vId] = g.coreNumber;
            }
            for(ui v : clique) {
                ui vId = reId[v]; stAdj[0][vId] = g.coreNumber;
            }
            for(ui v : nxtC) {
                ui vId = reId[v];
                for(ui j = g.pIdx2[v]; j < g.pIdx[v+1]; j++) {
                    ui w = g.pEdge[j];  if(reId[w] == g.n) continue;
                    if(Cv3.isIn(reId[w], 0, EC)) {
                        adj[vId][--stAdj[0][vId]] = reId[w];
                    }
                }
            }
            for(ui w : clique) {
                ui wId = reId[w];
                inC[wId] = true;
                for(ui j = g.pIdx2[w]; j < g.pIdx[w+1]; j++) {
                    ui v = g.pEdge[j];
                    ui vId = reId[v];
                    if(vId == g.n) continue;
                    if(Pv3.isIn(vId, 0, EP)) {
                        adj[vId][--stAdj[0][vId]] = wId;
                    }
                    // else if(Cv3.isIn(vId, 0, EC)) {
                    //     adj[vId][--stAdj[0][vId]] = wId;
                    //     adj[wId][--stAdj[0][wId]] = vId;
                    // }
                }
            }       
    
    
    #ifdef DEBUG
    ui preA = ans;
    #endif
            numOfNonCNodess[threadId] = 1;
            BKTomitaRecCCRV3(0, ER, EC, EP, SX, EX, &ans, threadId);
    #ifdef DEBUG
    // if(xx == xxx)
    printf("local:%u:%u, ans %u\n", u, ans - preA, ans);fflush(stdout);
    #endif
    
            for(ui w : clique) {
                ui wId = reId[w];
                inC[wId] = false;
            }
    
            for(ui i = g.pIdx[u]; i < g.pIdx[u + 1]; i++) {
                reId[g.pEdge[i]] = g.n;
            }
    
        };

        ui u;
        while((u = __sync_fetch_and_add(&pTs[threadId], 1)) < pTe[threadId]) {
            // if(g.pIdx[u+1] - g.pIdx2[u] < k-1) continue;
            computeU(u);
        }

        for(ui jjj = 0; jjj < nThreads; jjj++) {
            if(jjj == threadId) continue;
            if(pTs[jjj] >= pTe[jjj]) continue;
            while((u = __sync_fetch_and_add(&pTs[jjj], 1)) < pTe[jjj]) {
                // if(g.pIdx[u+1] - g.pIdx2[u] < k-1) continue;
                computeU(u);
            }
        }

        // #pragma omp for schedule(dynamic, 1) nowait
        // for(ui u = 0; u < g.n; u++) {
        //     computeU(u);
        // }
        cnt += ans;
    }

    return cnt;
}

ui parallelMCE::runCoreCliqueRemovalV3EdgeParallel() {
    printf("MCE::runCoreCliqueRemovalV3EdgeParallel\n");
#ifdef DEBUG
g.print();
#endif  

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

    // double threadTT[50];
    #pragma omp parallel reduction(+:cnt)
    {
        double tS = omp_get_wtime();

        int threadId = omp_get_thread_num();
        std::vector<ui> newCliqueNei;

        ListLinearHeap heap(g.coreNumber, g.coreNumber);
        heap.reserve(g.coreNumber, g.coreNumber);

        std::vector<ui> id_s(g.coreNumber+1);
        std::vector<ui> degree(g.coreNumber);
        std::vector<bool> removed(g.coreNumber);

        std::vector<ui> pIdx(g.coreNumber+1);
        std::vector<ui> eIdx(g.coreNumber+1);
        std::vector<ui> pEdge(g.coreNumber * g.coreNumber);

        std::vector<ui> clique(g.coreNumber);
        std::vector<ui> nxtC(g.coreNumber);

        LinearSet & V = Vs[threadId];
        std::vector<std::vector<ui>> & adj = adjs[threadId];
        std::vector<std::vector<ui>> & edAdj = edAdjs[threadId];
        std::vector<std::vector<ui>> & stAdj = stAdjs[threadId];
        std::vector<ui> & reId = reIds[threadId];
        std::vector<bool> & inC = inCs[threadId];

        // std::vector<bool> & vis = viss[threadId];

        LinearSet & Cv3 = Cv3s[threadId];
        LinearSet & Rv3 = Rv3s[threadId];
        LinearSet & Pv3 = Pv3s[threadId];
        LinearSet & Xv3 = Xv3s[threadId];

        ui ans = 0;

        // double tS = omp_get_wtime();
        
        #pragma omp for schedule(dynamic, 1)
        for(ui u = 0; u < g.n; u++) {
        #ifdef DEBUG
        xx=u;
        if(xx == xxx) 
        printf("\n    st %u\n", u);
        #endif  
        
            for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
                reId[g.pEdge[i]] = i-g.pIdx2[u]+1;
            }
    
            ui n = g.pIdx[u+1] - g.pIdx2[u];
            for(ui i = 0; i <= n; i++) pIdx[i] = 0;
    
            for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
                ui v = g.pEdge[i];
    
                for(ui j = g.pIdx2[v]; j < g.pIdx[v+1]; j++) {
                    ui w = g.pEdge[j];
                    if(reId[w] < g.n) {
                        pIdx[i-g.pIdx2[u]+1]++;
                        pIdx[reId[w]]++;
                    }
                }
            }
            for(ui i = 1; i <= n; i++) pIdx[i] += pIdx[i-1];
            for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
                ui v = g.pEdge[i];
                for(ui j = g.pIdx2[v]; j  < g.pIdx[v+1]; j++) {
                    ui w = g.pEdge[j];
                    if(reId[w] < g.n) {
                        pEdge[pIdx[i-g.pIdx2[u]]++] = reId[w]-1;
                        pEdge[pIdx[reId[w]-1]++] = i-g.pIdx2[u];
                    }
                }
            }
            for(ui i = n; i > 0; i--) pIdx[i] = pIdx[i-1];
            pIdx[0] = 0;
    
            clique.clear();
            nxtC.clear();
    
            ui queue_n = 0, new_size = 0;
            for(ui i = 0;i < n;i ++) removed[i] = false;
            for(ui i = 0;i < n;i ++) {
                if(pIdx[i+1]-pIdx[i] < threshold) id_s[queue_n ++] = i;
                degree[i] = pIdx[i+1]-pIdx[i];
            }
            for(ui i = 0;i < queue_n;i ++) {
                ui uu = id_s[i]; degree[uu] = 0;
                removed[uu] = true;
                for(ui j = pIdx[uu];j < pIdx[uu+1];j ++) if(degree[pEdge[j]] > 0) {
                    if((degree[pEdge[j]] --) == threshold) id_s[queue_n ++] = pEdge[j];
                }
            }
            
            for(ui i = 0;i < n;i ++) {
                if(degree[i] >= threshold) id_s[queue_n + (new_size ++)] = i;
                else {
                    removed[i] = true;
                }
            }
    
            if(new_size == 0) {
                for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
                    reId[g.pEdge[i]] = g.n;
                }
                continue;
            }
            else {
#ifdef DEBUG
if(xx == xxx) {
printf("queue_n:%u\n", queue_n);//should be zero by default
printf("To peel:");
for(ui i = 0; i < new_size; i++) {
    printf("%u-%u ", id_s[queue_n + i]+1, 
        degree[id_s[queue_n + i]]);
    
}
printf("\n");fflush(stdout);
}
#endif
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
    
                        // assert(x_size == new_size);
        // printf("	last key %u, x_size %u\n", key, x_size);
                        for(ui j = i;j < new_size;j ++) {
#ifdef DEBUG
printf("findClique, u %u, v %u, pClique %u\n", 
        u, g.pEdge[g.pIdx2[u]+id_s[queue_n+j]], pClique[u]);fflush(stdout);
#endif
                            // clique.push_back(g.pEdge[g.pIdx2[u]+id_s[queue_n+j]]);
                            coreCliquesBuffer[pClique[u]++] = g.pEdge[g.pIdx2[u]+id_s[queue_n+j]];
                        }
                        
                        for(ui j = 0; j < i; j++) {
                            ui tu = g.pEdge[g.pIdx2[u]+id_s[queue_n+j]];
                            nxtCandidates[pCandidate[u]].first = u;
                            nxtCandidates[pCandidate[u]].second = tu;
                            pCandidate[u]++;
                            // nxtC.push_back(tu);
                        }

                        bool extConnectAll = false;
                        for(ui jk = g.pIdx[u]; jk < g.pIdx2[u]; jk++) {
                            ui vtmp = g.pEdge[jk];
                            bool connectAll = true;
                            for(ui j = i;j < new_size;j ++) {
                                ui c = g.pEdge[g.pIdx2[u]+id_s[queue_n+j]];
                                if(!g.connectHash(vtmp, c)) {
                                    connectAll = false; break;
                                }
                            }
                            if(connectAll) {
                                extConnectAll = true;
                                break;
                            }
                        }
                        for(ui j = 0; j < i; j++) {
                            ui tu = g.pEdge[g.pIdx2[u]+id_s[queue_n+j]];
                            bool connectAll = true;
                            for(ui l = new_size-i-1; l >= 0; l--) {
                                ui c = coreCliquesBuffer[pClique[u]-l-1];
                                if(!g.connectHash(tu, c)) {
                                    connectAll = false; break;
                                }
                            }
                            if(connectAll) {
                                extConnectAll = true;
                                break;
                            }
                        }

                        if(!extConnectAll) {

                            ans += 1;
#ifdef DEBUG
printf("add 1, ans %u \n", ans);
#endif
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
    
    // printf("findClique\n");fflush(stdout);        
    
            
            for(ui i = g.pIdx[u]; i < g.pIdx[u + 1]; i++) {
                reId[g.pEdge[i]] = g.n;
            }
    
        }

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
            ui vv = nxtCandidates[e].second;

            for(ui i = g.pIdx[u]; i < g.pIdx[u + 1]; i++) {
                reId[g.pEdge[i]] = i - g.pIdx[u];
            }
            

            clique.clear();
            for(ui i = pSt[u]; i < pClique[u]; i++) {
                ui c = coreCliquesBuffer[i];
                if(g.connectHash(vv, c)) {
                    clique.push_back(c);
                }
            }
            for(auto u : clique) inC[reId[u]] = true;
            nxtC.clear();
            for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
                ui w = g.pEdge[i];
                if(w <= vv || inC[reId[w]]) continue;
                if(g.connectHash(vv, w)) {
                    nxtC.push_back(w);
                }
            }
    
#ifdef DEBUG
// if(xx == xxx) {

printf("\nnxtC: ");for(auto v:nxtC)printf("%u ", v); printf("\n");
printf("clique: ");for(auto v:clique)printf("%u ", v); printf("\n");
fflush(stdout);
// }
#endif
            

            ui EC = 0;
            ui EX = g.maxD, SX = g.maxD;
            ui ER = 0;
            ui EP = 0;
            for(ui i = g.pIdx[u]; i < g.pIdx[u+1]; i++)  {
                ui v = g.pEdge[i];
                if(v == vv) break;
                if(!g.connectHash(vv, v) || inC[reId[v]]) continue;
                Xv3.changeTo(reId[g.pEdge[i]], --SX);
// printf("pushX %u-%u\n", g.pEdge[i], reId[g.pEdge[i]]);
            }
            for(ui v : nxtC)  {
                if(!g.connectHash(vv, v)) continue;
                Pv3.changeTo(reId[v], EP++);
            }
            for(ui v : clique) {
                if(!g.connectHash(vv, v)) continue;
                Cv3.changeTo(reId[v], EC++);
            }
#ifdef DEBUG
// if(xx == xxx) {
printf("EC %u, EP %u, EX-SX %u\n", EC, EP, EX-SX);
fflush(stdout);
// }
#endif
            //X nodes, the neighbors in P
            // for(ui i = g.pIdx[u]; i < g.pIdx2[u]; i++) {
            //     ui v = g.pEdge[i];
            for(ui i = SX; i < EX; i++) {
                ui vId = Xv3[i];
                ui v = g.pEdge[g.pIdx[u] + vId];
                // ui vId = i - g.pIdx[u];
                // ui vId = reId[v];
// printf("nn u%u-vv%u-v%u-vid%u\n ", u, vv, v, vId);
                edAdj[0][vId] = 0;
                for(ui j = g.pIdx2[v]; j < g.pIdx[v+1]; j++) {
                    ui w = g.pEdge[j]; if(reId[w] == g.n) continue;
                    if(Pv3.isIn(reId[w], 0, EP)) {
                        adj[vId][edAdj[0][vId]++] = reId[w];
// printf("addXP:%u %u %u\n", w, reId[w], edAdj[0][vId]);
                    }
                }
            }
            //P/C nodes, the neighbors in P
            // for(ui i = g.pIdx2[u]; i < g.pIdx[u+1]; i++) {
            //     ui vId = i - g.pIdx[u]; edAdj[0][vId] = 0;
            // }
            for(ui v : clique) edAdj[0][ reId[v]] = 0;
            for(ui v : nxtC) edAdj[0][ reId[v]] = 0;
            for(ui v : nxtC) {
                ui vId = reId[v];
                for(ui j = g.pIdx2[v]; j < g.pIdx[v+1]; j++) {
                    ui w = g.pEdge[j];
                    ui wId = reId[w];  if(reId[w] == g.n) continue;
                    if(Pv3.isIn(wId, 0, EP)) {
                        adj[vId][edAdj[0][vId]++] = wId;
                        adj[wId][edAdj[0][wId]++] = reId[v];
                    }
                    else if(Cv3.isIn(wId, 0, EC)) {
                        adj[wId][edAdj[0][wId]++] = reId[v];
                    }
                }
            }
            for(ui v : clique) {
                ui vId = reId[v];
                for(ui j = g.pIdx2[v]; j < g.pIdx[v+1]; j++) {
                    ui w = g.pEdge[j];
                    ui wId = reId[w];  if(reId[w] == g.n) continue;
                    if(Pv3.isIn(wId, 0, EP)) {
                        adj[vId][edAdj[0][vId]++] = wId;
                    }
                }
            }
            //X nodes, the neighbors in C
            // for(ui i = g.pIdx[u]; i < g.pIdx2[u]; i++) {
            //     ui v = g.pEdge[i];
            for(ui i = SX; i < EX; i++) {
                ui vId = Xv3[i];
                // ui vId = reId[v];
                ui v = g.pEdge[g.pIdx[u] + vId];
                stAdj[0][vId] = g.coreNumber;
    
                // for(ui j = g.pIdx2[v]; j < g.pIdx[v+1]; j++) {
                //     ui w = g.pEdge[j];  if(reId[w] == g.n) continue;
                //     if(Cv3.isIn(reId[w], 0, EC)) {
                //         adj[vId][--stAdj[0][vId]] = reId[w];
                //     }
                // }
                for(ui c : clique) {
                    if(g.connectHash(v, c)) {
                        adj[vId][--stAdj[0][vId]] = reId[c];
                    }
                }
            }
            //P/C nodes, the neighbors in C
            for(ui v : nxtC) {
// printf("u%u-vv%u-v%u-vId%u-%u\n", u, vv, v, reId[v], stAdj[0].size());
                ui vId = reId[v]; stAdj[0][vId] = g.coreNumber;
            }
            for(ui v : clique) {
                ui vId = reId[v]; stAdj[0][vId] = g.coreNumber;
            }
            for(ui v : nxtC) {
                ui vId = reId[v];
                for(ui j = g.pIdx2[v]; j < g.pIdx[v+1]; j++) {
                    ui w = g.pEdge[j];  if(reId[w] == g.n) continue;
                    if(Cv3.isIn(reId[w], 0, EC)) {
                        adj[vId][--stAdj[0][vId]] = reId[w];
                    }
                }
            }
            for(ui w : clique) {
                ui wId = reId[w];
                // inC[wId] = true;
                for(ui j = g.pIdx2[w]; j < g.pIdx[w+1]; j++) {
                    ui v = g.pEdge[j];
                    ui vId = reId[v];
                    if(vId == g.n) continue;
                    if(Pv3.isIn(vId, 0, EP)) {
                        adj[vId][--stAdj[0][vId]] = wId;
                    }
                    // else if(Cv3.isIn(vId, 0, EC)) {
                    //     adj[vId][--stAdj[0][vId]] = wId;
                    //     adj[wId][--stAdj[0][wId]] = vId;
                    // }
                }
            }       
    
    
#ifdef DEBUG
ui preA = ans;
#endif
            numOfNonCNodess[threadId] = 1;
            BKTomitaRecCCRV3(0, ER, EC, EP, SX, EX, &ans, threadId);
#ifdef DEBUG
// if(xx == xxx)
printf("local:%u-%u:%u, ans %u\n", u, vv, ans - preA, ans);fflush(stdout);
#endif
        
            for(ui w : clique) {
                ui wId = reId[w];
                inC[wId] = false;
            }
            // for(auto u : clique) vis[u] = false;

            for(ui i = g.pIdx[u]; i < g.pIdx[u + 1]; i++) {
                reId[g.pEdge[i]] = g.n;
            }
    
        }
        

        // double tE = omp_get_wtime();

        // threadTT[threadId] =  tE-tS;

        cnt += ans;
    }
    

    // for(int i = 0; i < nThreads; i++) {
    //     std::cout << i << ':' << threadTT[i] << std::endl;
    // }

    return cnt;
}

ui parallelMCE::runCoreCliqueRemovalV3() {
    printf("MCE::runCoreCliqueRemovalV3\n");
    
#ifdef DEBUG
g.print();
#endif

    double threadTT[50];
    #pragma omp parallel reduction(+:cnt)
    {
        double tS = omp_get_wtime();

        int threadId = omp_get_thread_num();
        std::vector<ui> newCliqueNei;

        ListLinearHeap heap(g.coreNumber, g.coreNumber);
        heap.reserve(g.coreNumber, g.coreNumber);

        std::vector<ui> id_s(g.coreNumber+1);
        std::vector<ui> degree(g.coreNumber);
        std::vector<bool> removed(g.coreNumber);

        std::vector<ui> pIdx(g.coreNumber+1);
        std::vector<ui> eIdx(g.coreNumber+1);
        std::vector<ui> pEdge(g.coreNumber * g.coreNumber);

        std::vector<ui> clique(g.coreNumber);
        std::vector<ui> nxtC(g.coreNumber);

        LinearSet & V = Vs[threadId];
        std::vector<std::vector<ui>> & adj = adjs[threadId];
        std::vector<std::vector<ui>> & edAdj = edAdjs[threadId];
        std::vector<std::vector<ui>> & stAdj = stAdjs[threadId];
        std::vector<ui> & reId = reIds[threadId];
        std::vector<bool> & inC = inCs[threadId];

        LinearSet & Cv3 = Cv3s[threadId];
        LinearSet & Rv3 = Rv3s[threadId];
        LinearSet & Pv3 = Pv3s[threadId];
        LinearSet & Xv3 = Xv3s[threadId];

        ui ans = 0;

        // double tS = omp_get_wtime();
        
        #pragma omp for schedule(dynamic, 1) nowait
        for(ui u = 0; u < g.n; u++) {
#ifdef DEBUG
xx=u;
if(xx == xxx) 
printf("\n    st %u\n", u);
#endif  
        
            for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
                reId[g.pEdge[i]] = i-g.pIdx2[u]+1;
            }
    
            ui n = g.pIdx[u+1] - g.pIdx2[u];
            for(ui i = 0; i <= n; i++) pIdx[i] = 0;
    
            for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
                ui v = g.pEdge[i];
    
                for(ui j = g.pIdx2[v]; j < g.pIdx[v+1]; j++) {
                    ui w = g.pEdge[j];
                    if(reId[w] < g.n) {
                        pIdx[i-g.pIdx2[u]+1]++;
                        pIdx[reId[w]]++;
                    }
                }
            }
            for(ui i = 1; i <= n; i++) pIdx[i] += pIdx[i-1];
            for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
                ui v = g.pEdge[i];
                for(ui j = g.pIdx2[v]; j  < g.pIdx[v+1]; j++) {
                    ui w = g.pEdge[j];
                    if(reId[w] < g.n) {
                        pEdge[pIdx[i-g.pIdx2[u]]++] = reId[w]-1;
                        pEdge[pIdx[reId[w]-1]++] = i-g.pIdx2[u];
                    }
                }
            }
            for(ui i = n; i > 0; i--) pIdx[i] = pIdx[i-1];
            pIdx[0] = 0;
    
            clique.clear();
            nxtC.clear();
    
            ui queue_n = 0, new_size = 0;
            for(ui i = 0;i < n;i ++) removed[i] = false;
            for(ui i = 0;i < n;i ++) {
                if(pIdx[i+1]-pIdx[i] < threshold) id_s[queue_n ++] = i;
                degree[i] = pIdx[i+1]-pIdx[i];
            }
            for(ui i = 0;i < queue_n;i ++) {
                ui uu = id_s[i]; degree[uu] = 0;
                removed[uu] = true;
                for(ui j = pIdx[uu];j < pIdx[uu+1];j ++) if(degree[pEdge[j]] > 0) {
                    if((degree[pEdge[j]] --) == threshold) id_s[queue_n ++] = pEdge[j];
                }
            }
            
            for(ui i = 0;i < n;i ++) {
                if(degree[i] >= threshold) id_s[queue_n + (new_size ++)] = i;
                else {
                    removed[i] = true;
                }
            }
    
            if(new_size == 0) {
                for(ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
                    reId[g.pEdge[i]] = g.n;
                }
                continue;
            }
            else {
    #ifdef DEBUG
    if(xx == xxx) {
    printf("queue_n:%u\n", queue_n);//should be zero by default
    printf("To peel:");
    for(ui i = 0; i < new_size; i++) {
        printf("%u-%u ", id_s[queue_n + i]+1, 
            degree[id_s[queue_n + i]]);
        
    }
    printf("\n");fflush(stdout);
    }
    #endif
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
                            clique.push_back(g.pEdge[g.pIdx2[u]+id_s[queue_n+j]]);
                        }
                        
                        for(ui j = 0; j < i; j++) {
                            //gready extend
                            ui tu = g.pEdge[g.pIdx2[u]+id_s[queue_n+j]];
                            // bool connectAll = true;
                            // for(ui cu : clique) {
                            //     if(cu < tu) {
                            //         if(!g.connectOut(cu, tu)) {
                            //             connectAll = false; break;
                            //         }
                            //     }
                            //     else if(!g.connectOut(tu, cu)) {
                            //         connectAll = false; break;
                            //     }
                            // }
                            // if(connectAll)  clique.push_back(tu);
                            // else 
                                nxtC.push_back(tu);
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
    
    // printf("findClique\n");fflush(stdout);        
    
            for(ui i = g.pIdx[u]; i < g.pIdx[u + 1]; i++) {
                reId[g.pEdge[i]] = i - g.pIdx[u];
            }
    
    #ifdef DEBUG
    if(xx == xxx) {
    
    printf("nxtC: ");for(auto v:nxtC)printf("%u ", reId[v]); printf("\n");
    printf("clique: ");for(auto v:clique)printf("%u ", reId[v]); printf("\n");
    fflush(stdout);
    }
    #endif
            ui EC = 0;
            ui EX = g.maxD, SX = g.maxD;
            ui ER = 0;
            ui EP = 0;
            for(ui i = g.pIdx[u]; i < g.pIdx2[u]; i++)  {
                Xv3.changeTo(reId[g.pEdge[i]], --SX);
            }
            for(ui v : nxtC)  {
                Pv3.changeTo(reId[v], EP++);
            }
            for(ui v : clique) {
                Cv3.changeTo(reId[v], EC++);
            }
    #ifdef DEBUG
    if(xx == xxx) {
    printf("EC %u, EP %u, EX %u\n", EC, EP, EX);
    fflush(stdout);
    }
    #endif
            //X nodes, the neighbors in P
            for(ui i = g.pIdx[u]; i < g.pIdx2[u]; i++) {
                ui v = g.pEdge[i];
                ui vId = i - g.pIdx[u];
                edAdj[0][vId] = 0;
                for(ui j = g.pIdx2[v]; j < g.pIdx[v+1]; j++) {
                    ui w = g.pEdge[j]; if(reId[w] == g.n) continue;
                    if(Pv3.isIn(reId[w], 0, EP)) {
                        adj[vId][edAdj[0][vId]++] = reId[w];
                    }
                }
            }
            //P/C nodes, the neighbors in P
            for(ui i = g.pIdx2[u]; i < g.pIdx[u+1]; i++) {
                ui vId = i - g.pIdx[u]; edAdj[0][vId] = 0;
            }
            for(ui v : nxtC) {
                ui vId = reId[v];
                for(ui j = g.pIdx2[v]; j < g.pIdx[v+1]; j++) {
                    ui w = g.pEdge[j];
                    ui wId = reId[w];  if(reId[w] == g.n) continue;
                    if(Pv3.isIn(wId, 0, EP)) {
                        adj[vId][edAdj[0][vId]++] = wId;
                        adj[wId][edAdj[0][wId]++] = reId[v];
                    }
                    else if(Cv3.isIn(wId, 0, EC)) {
                        adj[wId][edAdj[0][wId]++] = reId[v];
                    }
                }
            }
            for(ui v : clique) {
                ui vId = reId[v];
                for(ui j = g.pIdx2[v]; j < g.pIdx[v+1]; j++) {
                    ui w = g.pEdge[j];
                    ui wId = reId[w];  if(reId[w] == g.n) continue;
                    if(Pv3.isIn(wId, 0, EP)) {
                        adj[vId][edAdj[0][vId]++] = wId;
                    }
                }
            }
            //X nodes, the neighbors in C
            for(ui i = g.pIdx[u]; i < g.pIdx2[u]; i++) {
                ui v = g.pEdge[i];
                ui vId = i - g.pIdx[u];
                stAdj[0][vId] = g.coreNumber;
    
                for(ui j = g.pIdx2[v]; j < g.pIdx[v+1]; j++) {
                    ui w = g.pEdge[j];  if(reId[w] == g.n) continue;
                    if(Cv3.isIn(reId[w], 0, EC)) {
                        adj[vId][--stAdj[0][vId]] = reId[w];
                    }
                }
            }
            //P/C nodes, the neighbors in C
            for(ui v : nxtC) {
                ui vId = reId[v]; stAdj[0][vId] = g.coreNumber;
            }
            for(ui v : clique) {
                ui vId = reId[v]; stAdj[0][vId] = g.coreNumber;
            }
            for(ui v : nxtC) {
                ui vId = reId[v];
                for(ui j = g.pIdx2[v]; j < g.pIdx[v+1]; j++) {
                    ui w = g.pEdge[j];  if(reId[w] == g.n) continue;
                    if(Cv3.isIn(reId[w], 0, EC)) {
                        adj[vId][--stAdj[0][vId]] = reId[w];
                    }
                }
            }
            for(ui w : clique) {
                ui wId = reId[w];
                inC[wId] = true;
                for(ui j = g.pIdx2[w]; j < g.pIdx[w+1]; j++) {
                    ui v = g.pEdge[j];
                    ui vId = reId[v];
                    if(vId == g.n) continue;
                    if(Pv3.isIn(vId, 0, EP)) {
                        adj[vId][--stAdj[0][vId]] = wId;
                    }
                    // else if(Cv3.isIn(vId, 0, EC)) {
                    //     adj[vId][--stAdj[0][vId]] = wId;
                    //     adj[wId][--stAdj[0][wId]] = vId;
                    // }
                }
            }       
    
    
    #ifdef DEBUG
    ui preA = ans;
    #endif
            numOfNonCNodess[threadId] = 1;
            BKTomitaRecCCRV3(0, ER, EC, EP, SX, EX, &ans, threadId);
    #ifdef DEBUG
    // if(xx == xxx)
    printf("local:%u:%u, ans %u\n", u, ans - preA, ans);fflush(stdout);
    #endif
    
            for(ui w : clique) {
                ui wId = reId[w];
                inC[wId] = false;
            }
    
            for(ui i = g.pIdx[u]; i < g.pIdx[u + 1]; i++) {
                reId[g.pEdge[i]] = g.n;
            }
    
        }
        

        double tE = omp_get_wtime();

        threadTT[threadId] =  tE-tS;

        cnt += ans;
    }
    

    for(int i = 0; i < nThreads; i++) {
        std::cout << i << ':' << threadTT[i] << std::endl;
    }

    return cnt;
}

// #define DEBUG
// ui xx, xxx;

// #undef DEBUG

void parallelMCE::BKTomitaRecCCRV3(ui deep, ui ER, ui EC, 
    ui EP, ui SX, ui EX, ui * ans, int threadId) {

    bool f = true;
    std::vector<ui> & ed = edAdjs[threadId][deep];
    std::vector<ui> & st = stAdjs[threadId][deep];
    std::vector<ui> & nxtEd = edAdjs[threadId][deep+1];
    std::vector<ui> & nxtSt = stAdjs[threadId][deep+1];

    LinearSet & Cv3 = Cv3s[threadId];
    LinearSet & Rv3 = Rv3s[threadId];
    LinearSet & Pv3 = Pv3s[threadId];
    LinearSet & Xv3 = Xv3s[threadId];

    std::vector<std::vector<ui>> & adj = adjs[threadId];
    std::vector<std::vector<ui>> & edAdj = edAdjs[threadId];
    std::vector<std::vector<ui>> & stAdj = stAdjs[threadId];

    std::vector<ui> & reId = reIds[threadId];
        std::vector<bool> & inC = inCs[threadId];

#ifdef DEBUG
// if(xx == xxx) {
printf("\n    deep %u\n", deep);fflush(stdout);
printf("C:"); for(ui i = 0; i < EC; i++) printf("%u ", Cv3[i]); printf("\n");
printf("R:"); for(ui i = 0; i < ER; i++) printf("%u ", Rv3[i]); printf("\n");
printf("P:"); for(ui i = 0; i < EP; i++) printf("%u ", Pv3[i]); printf("\n");
printf("X:"); for(ui i = SX; i < EX; i++) printf("%u ", Xv3[i]); printf("\n");
printf("P: P/C\n");
for(ui i = 0; i < EP; i++) {
    ui vId = Pv3[i];
    printf("%u: ", vId);
    for(ui j = 0; j < edAdj[deep][vId]; j++) {
        printf("%u ", adj[vId][j]);
    }
    printf(" / ");
    for(ui j = stAdj[deep][vId]; j < g.coreNumber; j++) {
        printf("%u ", adj[vId][j]);
    }
    printf("\n");
}

printf("PXsg:\n");
for(ui i = 0; i < EP; i++) {
    ui vId = Pv3[i];
    printf("%u: ", vId);
    for(ui j = SX; j < EX; j++) {
        if(g.connect(g.pEdge[g.pIdx[xx]+vId], g.pEdge[g.pIdx[xx]+Xv3[j]]))
            printf("%u ", Xv3[j]);
    }
    printf("\n");
}

printf("XPsg:\n");
for(ui i = SX; i < EX; i++) {
    ui vId = Xv3[i];
    printf("%u: ", vId);
    for(ui j = 0; j < ed[vId]; j++) {
        if(g.connect(g.pEdge[g.pIdx[xx]+vId], g.pEdge[g.pIdx[xx]+Xv3[j]]))
            printf("%u ", Xv3[j]);
    }
    printf("\n");
}

printf("PPsg:\n");
for(ui i = 0; i < EP; i++) {
    ui vId = Pv3[i];
    printf("%u: ", vId);
    for(ui j = 0; j < EP; j++) {
        printf("%u ", adj[vId][j]);
    }
    printf("\n");
}
printf("CPsg:\n");
for(ui i = 0; i < EC; i++) {
    ui vId = Cv3[i];
    printf("%u: ", vId);
    for(ui j = 0; j < edAdj[deep][vId]; j++) {
        printf("%u ", adj[vId][j]);
    }
    printf("\n");
}
printf("XCsg:\n");
for(ui i = SX; i < EX; i++) {
    ui vId = Xv3[i];
    printf("%u: ", vId);
    for(ui j = stAdj[deep][vId]; j < g.coreNumber; j++) {
        printf("%u ", adj[vId][j]);
    }
    printf("\n");
}
printf("CCsg:\n");
for(ui i = 0; i < EC; i++) {
    ui vId = Cv3[i];
    printf("%u: ", vId);
    for(ui j = stAdj[deep][vId]; j < g.coreNumber; j++) {
        printf("%u ", adj[vId][j]);
    }
    printf("\n");
}
// }
#endif

    // if(numOfNonCNodes == 0 && deep > 0) f = false;
    if(numOfNonCNodess[threadId] == 0) {  
        f = false; 
    }
    if(f && EX-SX + EP > 0) {
        for(ui i = 0; i < EP; i++) {
            ui v = Pv3[i];
            if(g.coreNumber-st[v] == EC) {
                f = false;
                break; 
            }
        }
        if(f) {
            for(ui i = SX; i < EX; i++) {
                ui v = Xv3[i];
                if( inC[v] || g.coreNumber-st[v] == EC) {
                    f = false;
                    break; 
                }
            }
        }
        

    }
#ifdef DEBUG
// if(xx == xxx) {
printf("F:%d\n", (int)f);
// }
#endif
    if(f) {
        // cnt++;
        (*ans)++;
    }

    if(EP == 0) return;

    if(EP == 1) {
        ui u = Pv3[0];
#ifdef DEBUG
// if(xx == xxx) {
printf("EP==1\n");
// }
#endif
        numOfNonCNodess[threadId] = 1;
        Rv3.changeTo(u, ER);
        ui newEC = 0;
        for(ui i = st[u]; i < g.coreNumber; i++) {
            Cv3.changeTo(adj[u][i], newEC++);
        }
        ui newEX = SX;
        for(ui i = SX; i < EX; i++) {
            ui v = Xv3[i]; 
            
            bool cc = false;
            for(ui j = 0; j < ed[v]; j++) {
                ui w = adj[v][j];
                if(w == u) {
                    cc = true; break;
                }
            }
            if(cc) Xv3.changeTo(v, newEX++);
        }
        auto mainNei2 = [&](ui v) {
            ui & nst = stAdj[deep+1][v];
            ui & ned = edAdj[deep+1][v];
            nst = g.coreNumber;
            ned = 0;
            if(newEC > 0)
            for(ui j = st[v]; j < nst; ) {
                if(Cv3.isIn(adj[v][j], 0, newEC)) {
                    std::swap(adj[v][j], adj[v][--nst]);
                }
                else j++;
            } 
        };
        
        for(ui i = SX; i < newEX; i++) {
            ui v = Xv3[i]; 
            if(inC[v]) continue;
            mainNei2(v);
        }

        BKTomitaRecCCRV3(deep + 1, ER+1, newEC, 0, SX, newEX, ans, threadId);
        //Rv3.changeTo
        return;
    }

    //remove the X nodes that has no neighbors in P
    for(ui i = SX; i < EX; ) {
        ui v = Xv3[i];
        if(ed[v] == 0) {
            Xv3.changeTo(v, --EX);
        }
        else i++;
    }

    //remove the C nodes that connects to all nodes in P
    bool update = false;
    for(ui i = 0; i < EC; ) {
        ui v = Cv3[i];
        if(ed[v] == EP) {
            Cv3.changeTo(v, --EC);
            Rv3.changeTo(v, ER++);
            update = true;
#ifdef DEBUG
// if(xx == xxx)
printf("remove C %u\n", v);
#endif
            for(ui i = SX; i < EX; ) {
                ui vv = Xv3[i];

                if(inC[vv]) {
                    i++;
                    continue;
                }

                bool cc = false;
                for(ui j = st[vv]; j < g.coreNumber; j++) {
                    ui w = adj[vv][j];
                    if(w == v) {
                        cc = true; break;
                    }
                }
                if(cc) i++;
                else Xv3.changeTo(vv, --EX);
            }
        }
        else i++;
    }
    if(update) {
        auto mainNeiInC = [&](ui v) {
            ui nst = g.coreNumber;

            if(EC > 0) {
                for(ui j = st[v]; j < nst; ) {
                    if(Cv3.isIn(adj[v][j], 0, EC)) {
                        std::swap(adj[v][j], adj[v][--nst]);
                    }
                    else j++;
                } 
            }
            
            st[v] = nst;
        };
        
        for(ui i = 0; i < EP; i++) {
            ui v = Pv3[i]; mainNeiInC(v);
        }
        for(ui i = SX; i < EX; i++) {
            ui v = Xv3[i]; 
            if(inC[v]) continue;
            mainNeiInC(v);
        }
        // for(ui i = 0; i < EC; i++) {
        //     ui v = Cv3[i]; mainNeiInC(v);
        // }
    }

    //remove the P nodes that connects to all PC nodes
    update = false;
    for(ui i = 0, rawEp = EP; i < EP; ) {
        ui v = Pv3[i];
        if(ed[v] == rawEp-1 && g.coreNumber-st[v] == EC) {
            Pv3.changeTo(v, --EP);
            Rv3.changeTo(v, ER++);
            update = true;
#ifdef DEBUG
// if(xx == xxx)
printf("remove P %u\n", v);
#endif
            for(ui i = SX; i < EX; ) {
                ui vv = Xv3[i];
                bool cc = false;
                for(ui j = 0; j < ed[vv]; j++) {
                    ui w = adj[vv][j];
                    if(w == v) {
                        cc = true; break;
                    }
                }
                if(cc) i++;
                else Xv3.changeTo(vv, --EX);
            }
        }
        else {
            i++;
        }
    }
    if(update) {
        auto mainNeiInP = [&](ui v) {
            for(ui j = 0; j < ed[v]; ) {
                if(Pv3.isIn(adj[v][j], 0, EP)) {
                    j++;
                }
                else std::swap(adj[v][j], adj[v][--ed[v]]);
            }
        };
        
        for(ui i = 0; i < EP; i++) {
            ui v = Pv3[i]; mainNeiInP(v);
        }
        for(ui i = SX; i < EX; i++) {
            ui v = Xv3[i]; mainNeiInP(v);
        }
        for(ui i = 0; i < EC; i++) {
            ui v = Cv3[i]; mainNeiInP(v);
        }

        bool f = true;
        for(ui i = 0; i < EP; i++) {
            ui v = Pv3[i];
            if(g.coreNumber-st[v] == EC) {
                f = false;
                break; 
            }
        }
        if(f)
        for(ui i = SX; i < EX; i++) {
            ui v = Xv3[i];
            if(inC[v] || g.coreNumber-st[v] == EC) {
                f = false;
                break; 
            }
        }
        if(f) {
            // cnt++;
            (*ans)++;
        }

        if(EP == 0) return;
    }

    
    //find the best pivot
    ui pu = Pv3[0], maxD = 0;
    for(ui i = 0; i < EP; i++) {
        ui v = Pv3[i];
        ui d = ed[v] + g.coreNumber-st[v];
        if(d > maxD) {
            maxD = d; pu = v;
        }
    }
    for(ui i = 0; i < EC; i++) {
        ui v = Cv3[i];
        // ui d = ed[v] + g.coreNumber-st[v];
        ui d = ed[v] + EC-1;
        if(d > maxD) {
            maxD = d; pu = v;
        }
    }
    for(ui i = SX; i < EX; i++) {
        ui v = Xv3[i];
        ui neiInC = inC[v]?EC-1:g.coreNumber-st[v];
        ui d = ed[v] + neiInC;
        if(d > maxD) {
            maxD = d; pu = v;
        }
        // continue;
        if(ed[v] == EP-1 && neiInC == EC) {
            //remove the X nodes that has one-non-neighbor in P
            
            for(ui i = 0, ttt = 0; i < ed[v]; i++) {
                ui tmpv = adj[v][i];
                Pv3.changeTo(tmpv, ttt++);
            }
            ui u = Pv3[EP-1];

            ui newEP = 0;
            ui newEX = SX;
            ui newEC = 0;
    #ifdef DEBUG
    if(xx == xxx) {
    printf("XoneNonNei:accNonNei %u, deep %u\n", u, deep);fflush(stdout);
    printf("Xv3 %u:", deep);
    for(ui i = SX; i < EX; i++) printf("%u ", Xv3[i]);
    printf("\n");
    }
    #endif
            //build new P
            for(ui i = 0; i < ed[u]; i++) {
                ui v = adj[u][i];
                if(Pv3.isin(v, 0, EP)) Pv3.changeTo(v, newEP++);
            }
            //build new C
            for(ui i = st[u]; i < g.coreNumber; i++) {
                ui v = adj[u][i];
                if(Cv3.isin(v, 0, EC)) Cv3.changeTo(v, newEC++);
            }
    #ifdef DEBUG
    if(xx == xxx) {
    printf("newEP %u, newEC %u\n", newEP, newEC);fflush(stdout);
    }
    #endif
            //build new X
            numOfNonCNodess[threadId] = 1;
            for(ui i = SX; i < EX; i++) {
                ui v = Xv3[i];
                bool cc = false;
                for(ui j = 0; j < ed[v]; j++) {
                    ui w = adj[v][j];
                    if(w == u) {
                        cc = true; break;
                    }
                }
                if(cc) Xv3.changeTo(v, newEX++);
            }
    #ifdef DEBUG
    if(xx == xxx) {
    printf("newEX %u\n", newEX-SX);fflush(stdout);
    }
    #endif

            // R.push_back(u);

            //maintain the neighbors in P / C
            auto mainNei = [&](ui v) {
                
                ui & ned = nxtEd[v];
                // nst = g.coreNumber;
                // ned = 0;
                // if(EP > 0)
                // for(ui j = 0; j < ed[v]; j++) {
                //     if(Pv3.isIn(adj[v][j], 0, newEP)) {
                //         std::swap(adj[v][j], adj[v][ned++]);
                //     }
                // }
                // if(EC > 0)
                // for(ui j = st[v]; j < nst; ) {
                //     if(Cv3.isIn(adj[v][j], 0, newEC)) {
                //         std::swap(adj[v][j], adj[v][--nst]);
                //     }
                //     else j++;
                // } 
                
                ned = ed[v];
                
                for(ui j = 0; j < ned; ) {
                    if(Pv3.isIn(adj[v][j], 0, newEP)) {
                        j++;
                    }
                    else std::swap(adj[v][j], adj[v][--ned]);
                }

                if(inC[v]) return;

                ui & nst = nxtSt[v];
                nst = st[v];
                for(ui j = st[v]; j < g.coreNumber; j++) {
                    if(!Cv3.isIn(adj[v][j], 0, newEC)) {
                        std::swap(adj[v][j], adj[v][nst++]);
                    }
                } 
            };
            
            for(ui i = 0; i < newEP; i++) {
                ui v = Pv3[i]; mainNei(v);
            }
            for(ui i = SX; i < newEX; i++) {
                ui v = Xv3[i]; mainNei(v);
            }
            for(ui i = 0; i < newEC; i++) {
                ui v = Cv3[i]; mainNei(v);
            }
            
            Rv3.changeTo(u, ER);

            BKTomitaRecCCRV3(deep + 1, ER+1, newEC, newEP, SX, newEX, ans, threadId);

            
            numOfNonCNodess[threadId]--;

            return;
        }
    }
#ifdef DEBUG
if(xx == xxx) {
printf("pu:%u,maxD %u, EC %u, EP %u, ed[pu] %u\n", pu, maxD, EC, EP, ed[pu]);fflush(stdout);
}
#endif
    std::vector<ui> nonNeighbors(EC+EP - maxD + 1);
    ui d = 0;
    ui newEP = 0;
    for(ui i = 0; i < ed[pu]; i++) {
        ui v = adj[pu][i];
        Pv3.changeTo(v, newEP++);
    }
#ifdef DEBUG
if(xx == xxx) {
printf("EP-newEP:%u,EC+EP - maxD %u\n", EP-newEP, EC+EP - maxD);fflush(stdout);
}
#endif
    for(ui i = newEP; i < EP; i++) {
        nonNeighbors[d++] = Pv3[i];
    }
    ui newEC = 0;
    if(Cv3.isIn(pu, 0, EC)) {
        newEC = EC-1;
        Cv3.changeTo(pu, EC-1);
    }
    else if(inC[pu]) {
        newEC = EC;
    }
    else {
        for(ui i = st[pu]; i < g.coreNumber; i++) {
            ui v = adj[pu][i];
            Cv3.changeTo(v, newEC++);
        }
    }
    for(ui i = newEC, dc = 0; i < EC; i++) {
        ui v = Cv3[i];

        //与pu的邻居相连
        bool f = false;
        for(ui j = 0; j < ed[v]; j++) {
            ui w = adj[v][j];
            if(Pv3.isIn(w, 0, newEP)) {
                f = true;
                break;
            }
        }
        if(!f) {
            continue;
        }
        // f = false;
        // //包含前面C点的P中的非邻居
        // for(ui j = d-dc; j < d; j++) {
        //     ui preCv = nonNeighbors[j];
                
        // }

        
        nonNeighbors[d++] = v;
        dc++;
    }
    
#ifdef DEBUG
if(xx == xxx) {
printf("nonNei:");
for(ui ii = 0; ii < d; ii++) printf("%u ", nonNeighbors[ii]);
printf("\n");fflush(stdout);
}
#endif
    
    for(ui ii = 0; ii < d; ii++) {
        ui u = nonNeighbors[ii];
        ui newEP = 0;
        ui newEX = SX;
        ui newEC = 0;
#ifdef DEBUG
if(xx == xxx) {
printf("accNonNei %u, deep %u\n", u, deep);fflush(stdout);
printf("Xv3 %u:", deep);
for(ui i = SX; i < EX; i++) printf("%u ", Xv3[i]);
printf("\n");
printf("Cv3 %u:", deep);
for(ui i = 0; i < EC; i++) printf("%u ", Cv3[i]);
printf("\n");
}
#endif
        //build new P
        for(ui i = 0; i < ed[u]; i++) {
            ui v = adj[u][i];
            if(Pv3.isin(v, 0, EP)) Pv3.changeTo(v, newEP++);
        }
        //build new C
        if(Pv3.isIn(u, 0, EP))
        for(ui i = st[u]; i < g.coreNumber; i++) {
            ui v = adj[u][i];
            if(Cv3.isin(v, 0, EC)) Cv3.changeTo(v, newEC++);
        }
        else {
            Cv3.changeTo(u, EC-1);
            newEC = EC-1;
        }
#ifdef DEBUG
if(xx == xxx) {
printf("newEP %u, newEC %u\n", newEP, newEC);fflush(stdout);
}
#endif
#ifdef DEBUG
if(xx == xxx) {
printf("accNonNei %u, deep %u\n", u, deep);fflush(stdout);
printf("Xv3 %u:", deep);
for(ui i = SX; i < EX; i++) printf("%u ", Xv3[i]);
printf("\n");
printf("Cv3 %u:", deep);
for(ui i = 0; i < EC; i++) printf("%u ", Cv3[i]);
printf("\n");
}
#endif
        //build new X
        if(Pv3.isIn(u, 0, EP)) {
            numOfNonCNodess[threadId] = 1;
            for(ui i = SX; i < EX; i++) {
                ui v = Xv3[i];
                bool cc = false;
                for(ui j = 0; j < ed[v]; j++) {
                    ui w = adj[v][j];
                    if(w == u) {
                        cc = true; break;
                    }
                }
                if(cc) Xv3.changeTo(v, newEX++);
            }
        }
        else {
            numOfNonCNodess[threadId] = 0;
            for(ui i = SX; i < EX; i++) {
                ui v = Xv3[i];

                if(inC[v]) {
#ifdef DEBUG
if(xx == xxx) {
printf("xcv %u\n", v);fflush(stdout);
}
#endif
                    Xv3.changeTo(v, newEX++);
                    continue;
                }

                bool cc = false;
                for(ui j = st[v]; j < g.coreNumber; j++) {
                    ui w = adj[v][j];
                    if(w == u) {
                        cc = true; break;
                    }
                }
                if(cc) {
#ifdef DEBUG
if(xx == xxx) {
printf("xv %u\n", v);fflush(stdout);
}
#endif
                    Xv3.changeTo(v, newEX++);
                }
            }
        }
#ifdef DEBUG
if(xx == xxx) {
printf("newEX %u\n", newEX-SX);fflush(stdout);
}
#endif

        // R.push_back(u);

        //maintain the neighbors in P / C
        auto mainNei = [&](ui v) {
            ui & nst = nxtSt[v];
            ui & ned = nxtEd[v];
            nst = g.coreNumber;
            ned = 0;

            if(EP > 0)
            for(ui j = 0; j < ed[v]; j++) {
                if(Pv3.isIn(adj[v][j], 0, newEP)) {
                    std::swap(adj[v][j], adj[v][ned++]);
                }
            }

            if(inC[v]) return;

            if(EC > 0)
            for(ui j = st[v]; j < nst; ) {
                if(Cv3.isIn(adj[v][j], 0, newEC)) {
                    std::swap(adj[v][j], adj[v][--nst]);
                }
                else j++;
            } 
        };
        
        for(ui i = 0; i < newEP; i++) {
            ui v = Pv3[i]; mainNei(v);
        }
        for(ui i = SX; i < newEX; i++) {
            ui v = Xv3[i]; mainNei(v);
        }
        for(ui i = 0; i < newEC; i++) {
            ui v = Cv3[i]; //mainNei(v);
           
            ui & ned = nxtEd[v];
            ned = 0;
            if(EP > 0)
            for(ui j = 0; j < ed[v]; j++) {
                if(Pv3.isIn(adj[v][j], 0, newEP)) {
                    std::swap(adj[v][j], adj[v][ned++]);
                }
            }
        }
        
        Rv3.changeTo(u, ER);

        BKTomitaRecCCRV3(deep + 1, ER+1, newEC, newEP, SX, newEX, ans, threadId);

        Xv3.changeTo(u, --SX);
#ifdef DEBUG
if(xx == xxx) {
printf("Xv322 %u:", deep);
for(ui i = SX; i < EX; i++) printf("%u ", Xv3[i]);
printf("\n");
}
#endif
        if(Pv3.isIn(u, 0, EP)) {
            numOfNonCNodess[threadId]--;
            Pv3.changeTo(u, --EP);
        }
        else {
            --EC;

            // Cv3.changeTo(u, --EC);
        }
    }


    
    for(ui ii = 0; ii < d; ii++) {
        Xv3.changeTo(nonNeighbors[ii], SX++);
        // if(inC[nonNeighbors[ii]]) {
        //     Cv3.changeTo(nonNeighbors[ii], EC++);
        // }
    }

#ifdef DEBUG
if(xx == xxx) {
printf("Xvlast %u:", deep);
for(ui i = SX; i < EX; i++) printf("%u ", Xv3[i]);
printf("\n");
printf("Cvlast %u:", deep);
for(ui i = 0; i < EC; i++) printf("%u ", Cv3[i]);
printf("\n\n");
}
#endif
}