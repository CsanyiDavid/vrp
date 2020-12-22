//
// Created by david on 2020. 12. 22..
//


#include "vrp.h"

using namespace std;
using namespace lemon;

int main(){
    VRP vrp("hun-sc-ncn.lgf");
    vrp.generateCostumersGraph(50);
    return 0;
}
