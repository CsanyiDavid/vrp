//
// Created by david on 2020. 12. 22..
//

#include "vrp.h"

using namespace std;
using namespace lemon;

bool readInitData(int &costumerCnt, int &seed, int &maxWeight){
    if(PRINT) cout << "costumer cnt: ";
    cin >> costumerCnt;
    if(cin.fail()){
        cout << "Invalid value" << endl;
        cin.clear();
        cin.ignore(1000, '\n');
        return false;
    }
    cin.ignore(1000, '\n');
    if(PRINT) cout << "random seed: ";
    cin >> seed;
    if(cin.fail()){
        cout << "Invalid value" << endl;
        cin.clear();
        cin.ignore(1000, '\n');
        return false;
    }
    cin.ignore(1000, '\n');
    if(PRINT) cout << "max weight: ";
    cin >> maxWeight;
    if(cin.fail()){
        cout << "Invalid value" << endl;
        cin.clear();
        cin.ignore(1000, '\n');
        return false;
    }
    cin.ignore(1000, '\n');
    if(PRINT) cout << endl;
    return true;
}

int main(){
    int costumerCnt;
    int seed;
    int maxWeight;
    VRP vrp("hun-sc-ncn.lgf");
    if(PRINT) cout << "(Type help for the available commands)" << endl;
    string input;
    //Input and initialization check is done in the vrp class
    while(true){
        cout << endl << " >> ";
        cin >> input;
        if(input=="exit" || input=="e"){
            break;
        } else if(input=="init" || input=="i") {
            if(readInitData(costumerCnt, seed, maxWeight)) {
                vrp.init(costumerCnt, seed, maxWeight);
            }
        } else if(input=="check" || input=="c") {
            vrp.checkMIP(false);
        } else if(input=="checkprintsolutionarcs"){
            vrp.checkMIP(false, true);
        } else if (input=="checkwithconditions" || input=="cc"){
            cout << "Enter the conditions or invalid value to end" << endl;
            int sId, tId, value;
            vector<tuple<int, int, int>> conditions;
            bool continueRead=true;
            while(continueRead){
                cout << "source: ";
                cin >> sId;
                if(cin.fail()){
                    cout << "Invalid value" << endl;
                    cin.clear();
                    continueRead=false;
                }
                cin.ignore(1000, '\n');
                if(continueRead) {
                    cout << "target: ";
                    cin >> tId;
                    if (cin.fail()) {
                        cout << "Invalid value" << endl;
                        cin.clear();
                        continueRead = false;
                    }
                    cin.ignore(1000, '\n');
                }
                if(continueRead) {
                    cout << "value(0/1): ";
                    cin >> value;
                    if (cin.fail()) {
                        cout << "Invalid value" << endl;
                        cin.clear();
                        continueRead = false;
                    }
                    cin.ignore(1000, '\n');
                }
                if(continueRead) {
                    conditions.push_back(tuple<int, int, int>(sId, tId, value));
                }
            }
            vrp.checkMIP(false, false, conditions);
        } else if(input=="print" || input=="p") {
            string name;
            cout << "Eps file's name: ";
            cin >> name;
            vrp.printToEps(name);
        } else if(input=="printsolution" || input=="ps"){
            vrp.printSolution();
       } else if(input=="printcost" || input=="pc"){
            bool error=false;
            int sId;
            int tId;
            cout << "Source ID: ";
            cin >> sId;
            if(cin.fail()){
                cin.clear();
                error=true;
            }
            cin.ignore(1000, '\n');
            cout << "Target ID: ";
            cin >> tId;
            if(cin.fail()){
                cin.clear();
                error=true;
            }
            cin.ignore(1000, '\n');
            if(!error) {
                vrp.printCost(sId, tId);
            } else {
                cerr << "Invalid value" << endl;
            }
        } else if (input=="branchandprice" || input=="bap") {
            vrp.callBranchAndPrice();
        } else if(input=="ir") {
            if(readInitData(costumerCnt, seed, maxWeight)) {
                vrp.init(costumerCnt, seed, maxWeight);
                vrp.callBranchAndPrice();
            }
        } else if (input=="help" || input=="h"){
                cout << "exit(e): \t exit from the program" << endl;
                cout << "init(i): \t initialize the problem with costumer count, random seed" << endl;
                cout << "\t\t and maximum demand of a costumer" << endl;
                cout << "check(c): \t solve the problem with the CPLEX MIP" << endl;
                cout << "checkprintsolutionarcs" << endl;
                cout << "checkwithconditions(cc)" << endl;
                cout << "print(p): \t print the map and the found solution to an eps file" << endl;
                cout << "printsolution(ps)" << endl;
                cout << "printcost(pc):\t print the cost between two nodes" << endl;
                cout << "branchandprice(bap)" << endl;
                cout << "help(h): \t print the list of commands" << endl;
        } else {
            cout << "unkown command" << endl;
        };
    }
    return 0;
}
