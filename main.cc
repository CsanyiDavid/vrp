//
// Created by david on 2020. 12. 22..
//


#include "vrp.h"

int mySeed;
int myMaxWeight;

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
         if(input=="exit" || input=="e"){
             break;
         } else if(input=="init" || input=="i") {
             cout << "cnt: ";
             cin >> costumerCnt;
             cout << "seed: ";
             cin >> mySeed;
             cout << "max weight: ";
             cin >> myMaxWeight;
             cout << endl;
             vrp.init(costumerCnt);
         } else if(input=="run" || input=="r") {
             vrp.branchAndBound();
         } else if(input=="initrun" || input=="ir"){
             cout << "cnt: ";
             cin >> costumerCnt;
             cout << "seed: ";
             cin >> mySeed;
             cout << "max weight: ";
             cin >> myMaxWeight;
             cout << endl;
             vrp.init(costumerCnt);
             vrp.branchAndBound();
         } else if(input=="check" || input=="c") {
             vrp.checkMIP();
         }else if(input=="print" || input=="p"){
             vrp.printToEps("OUT.eps");
         } else if(input=="help" || input=="h"){
             cout << "exit (e)" << endl;
             cout << "init (i)" << endl;
             cout << "run (r)" << endl;
             cout << "initrun (ir)" << endl;
             cout << "check (c)" << endl;
             cout << "print (p)" << endl;
             cout << "help (h)" << endl;
         } else {
             cout << "unkown command" << endl;
         };
    }
    return 0;
}
