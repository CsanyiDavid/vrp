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
        cout << "(Type help for the available commands)" << endl;
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
             if(costumerCnt==-1 || mySeed==-1 || myMaxWeight==-1){
                 cout << "aborted" << endl;
             } else {
                 vrp.init(costumerCnt);
                 vrp.branchAndBound();
             }
         } else if(input=="check" || input=="c") {
             vrp.checkMIP();
         }else if(input=="print" || input=="p") {
             vrp.printToEps("OUT.eps");
         } else if(input=="printroutes" || input=="pr"){
             vrp.printRoutes();
         } else if(input=="help" || input=="h"){
             cout << "exit (e): \t exit from the program" << endl;
             cout << "init (i): \t initialize the problem with costumer count, random seed" << endl;
             cout << "\t\t and maximum demand of a costumer" << endl;
             cout << "run (r): \t run the column generation and branch and bound algorithm" << endl;
             cout << "\t\t to solve the problem (init must be called first!)" << endl;
             cout << "initrun (ir): \t init and run" << endl;
             cout << "check (c): \t solve the problem with the CPLEX MIP(init must be called first!)" << endl;
             cout << "print (p): \t print the map and the found solution to an eps file" << endl;
             cout << "printroutes(pr):\t print the column generated routes' nodes and cost" << endl;
             cout << "help (h): \t print the list of commands" << endl;
         } else {
             cout << "unkown command" << endl;
         };
    }
    return 0;
}
