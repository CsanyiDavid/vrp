//
// Created by david on 2020. 12. 22..
//


#include "vrp.h"

int mySeed;

using namespace std;
using namespace lemon;

#define IS_MAP true

int main(int argc, char** argv){
    int costumerCnt;
    ArgParser ap(argc, argv);
    ap.refOption("cnt", "Costumer count", costumerCnt, true);
    ap.refOption("seed", "Random seed", mySeed, true);
    ap.run();
    #if IS_MAP
        VRP vrp(IS_MAP, "hun-sc-ncn.lgf", costumerCnt);
        //vrp.printToEps("graph.eps");
        vrp.createMasterLP();
        vrp.solveMasterLP();
        vrp.checkMIP(false);
    #else
        VRP vrp(IS_MAP, "inp.lgf");
        vrp.createMasterLP();
        vrp.solveMasterLP();
        vrp.checkMIP();
    #endif
    return 0;
}
