//
// Created by david on 2020. 12. 22..
//


#include "vrp.h"

using namespace std;
using namespace lemon;

int main(){
    //VRP vrp(true, "hun-sc-ncn.lgf");
    //vrp.generateCostumersGraph(50);
    //vrp.printCostumerCoordinates();
    //vrp.shortestPaths();
    //vrp.printShortestPathsFromDepot();
    //vrp.printToEps("graph.eps");

    VRP vrp(false, "inp.lgf");
    vrp.createMasterLP();
    vrp.printMasterLPSolution();
    return 0;
}
