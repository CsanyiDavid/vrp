//
// Created by david on 2020. 12. 22..
//


#include "vrp.h"

using namespace std;
using namespace lemon;

#define IS_MAP true

int main(){
    #if IS_MAP
        VRP vrp(IS_MAP, "hun-sc-ncn.lgf", 13);
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
