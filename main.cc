//
// Created by david on 2020. 12. 22..
//


#include "vrp.h"

using namespace std;
using namespace lemon;

#define IS_MAP true

int main(){
    #if IS_MAP
        VRP vrp(IS_MAP, "hun-sc-ncn.lgf");
        vrp.generateCostumersGraph(10);
        vrp.printCostumerCoordinates();
        vrp.shortestPaths();
        vrp.printShortestPathsFromDepot();
        //vrp.printToEps("graph.eps");
    #else
        VRP vrp(IS_MAP, "inp.lgf");
    #endif

    vrp.createMasterLP();
    vrp.solveMasterLP();

    vrp.checkLP();
    return 0;
}
