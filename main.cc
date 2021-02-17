//
// Created by david on 2020. 12. 22..
//
/*
 * Kérdések
 * node id k nem változnak meg?
 */

#include "vrp.h"

int mySeed;
int myMaxWeight;

using namespace std;
using namespace lemon;

void readInitData(int &costumerCnt, int &seed, int &maxWeight){
    if(PRINT) cout << "costumer cnt: ";
    cin >> costumerCnt;
    if(PRINT) cout << "random seed: ";
    cin >> seed;
    if(PRINT) cout << "max weight: ";
    cin >> maxWeight;
    if(PRINT) cout << endl;
}

int main(){
    int costumerCnt;
    int seed;
    int maxWeight;
    VRP vrp(true, "hun-sc-ncn.lgf");
    if(PRINT) cout << "(Type help for the available commands)" << endl;
    string input;

    //Input and initialization check is done in the vrp class
    while(true){
        cout << endl << " >> ";
        cin >> input;
        if(input=="exit" || input=="e"){
            break;
        } else if(input=="init" || input=="i") {
            readInitData(costumerCnt, seed, maxWeight);
            vrp.init(costumerCnt, seed, maxWeight);
        } else if(input=="check" || input=="c") {
            vrp.checkMIP(false);
        }else if(input=="print" || input=="p"){
            string name;
            cout << "Eps file's name: ";
            cin >> name;
            vrp.printToEps(name);
        /*} else if(input=="printroutes" || input=="pr") {
            vrp.printRoutes();
       */} else if(input=="printcost" || input=="pc"){
            int sId;
            int tId;
            cout << "Source ID: ";
            cin >> sId;
            cin.clear();
            cin.ignore(1000, '\n');
            cout << "Target ID: ";
            cin >> tId;
            cin.clear();
            cin.ignore(1000, '\n');
            vrp.printCost(sId, tId);
        } else if (input=="branchandprice" || input=="bap") {
            vrp.callBranchAndPrice();
        } else if(input=="clarkewright"  || input=="cw"){
            vrp.callClarkeWright();
        } else if (input=="help" || input=="h"){
                cout << "exit(e): \t exit from the program" << endl;
                cout << "init(i): \t initialize the problem with costumer count, random seed" << endl;
                cout << "\t\t and maximum demand of a costumer" << endl;
                cout << "check(c): \t solve the problem with the CPLEX MIP(init must be called first!)" << endl;
                cout << "print(p): \t print the map and the found solution to an eps file" << endl;
                cout << "printroutes(pr): print the column generated routes' nodes and cost" << endl;
                cout << "printcost(pc):\t print the cost between two nodes" << endl;
                cout << "branchandprice(bap)" << endl;
                cout << "help(h): \t print the list of commands" << endl;
        } else {
            cout << "unkown command" << endl;
        };
        cin.clear();
        cin.ignore(1000, '\n');
    }
    return 0;
}
