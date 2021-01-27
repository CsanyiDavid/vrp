//
// Created by david on 2020. 12. 22..
//


#include "vrp.h"

int mySeed;

using namespace std;
using namespace lemon;

#define IS_MAP true

int main(){
    int costumerCnt;
    //ArgParser ap(argc, argv);
    //ap.refOption("isMap", "Is map?", isMap, true);
    //ap.refOption("cnt", "Costumer count", costumerCnt, true);
    //ap.refOption("seed", "Random seed", mySeed, true);
    //ap.run();
    #if IS_MAP
        VRP vrp(IS_MAP, "hun-sc-ncn.lgf");
        //vrp.printToEps("graph.eps");
    #else
        VRP vrp(IS_MAP, "inp.lgf");
        vrp.createMasterLP();
        vrp.solveMasterLP();
        vrp.checkMIP();
    #endif

    string input;
    while(true){
        cout << endl << " >> ";
         cin >> input;
         if(input=="exit"){
             break;
         } else if(input=="init") {
             cout << "cnt: ";
             cin >> costumerCnt;
             cout << "seed: ";
             cin >> mySeed;
             cout << endl;
             vrp.init(costumerCnt);
         } else if(input=="run") {
            vrp.branchAndBound();
         } else if(input=="check") {
             vrp.checkMIP();
         }else if(input=="print"){
             vrp.printToEps("OUT.eps");
         } else if(input=="help"){
             cout << "exit" << endl;
             cout << "init" << endl;
             cout << "run" << endl;
             cout << "check" << endl;
             cout << "print" << endl;
             cout << "help" << endl;
         } else {
             cout << "unkown command" << endl;
         };
    }
    return 0;
}
