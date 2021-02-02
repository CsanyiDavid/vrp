//
// Created by david on 2021. 01. 13..
//

#include <lemon/lp.h>
#include <lemon/list_graph.h>
#include <vector>

using namespace std;
using namespace lemon;

void myAssert(bool bo, const string& errorType)
{
    if(!bo) {
        cout << "Assertion error: " << errorType << endl;
        exit(-1);
    }
}

int main()
{
    ListDigraph g;
    vector<ListDigraph::Node> nodes;
    for(int i=0; i<5; ++i){
        nodes.push_back(g.addNode());
    }
    ListDigraph::Arc arc;
    g.addArc(nodes[2], nodes[3]);
    arc=g.addArc(nodes[0], nodes[1]);
    ListDigraph copy;
    digraphCopy(g, copy).run();
    cout << countNodes(g) << endl;
    cout << countNodes(copy) << endl;
    cout << countArcs(g) << endl;
    cout << countArcs(copy) << endl;
    g.erase(arc);
    g.erase(nodes[1]);
    cout << countNodes(g) << endl;
    cout << countNodes(copy) << endl;
    cout << countArcs(g) << endl;
    cout << countArcs(copy) << endl;
    digraphCopy(copy, g).run();
    cout << countNodes(g) << endl;
    cout << countNodes(copy) << endl;
    cout << countArcs(g) << endl;
    cout << countArcs(copy) << endl;
    g.addArc(nodes[3], nodes[4]);
    digraphCopy(copy, g).run();
    cout << g.id(g.source(arc)) << " " << g.id(g.target(arc)) << endl;
    cout << copy.id(copy.source(arc)) << " " << copy.id(copy.target(arc)) << endl;
}
