//
// Created by david on 2020. 12. 22..
//
#include <iostream>
#include <lemon/list_graph.h>

using namespace std;
using namespace lemon;

int main(){
    cout << "Hello" << endl;
    ListDigraph g;
    g.addNode();
    cout << countNodes(g) << endl;
    return 0;
}
