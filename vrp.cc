//
// Created by david on 2020. 12. 22..
//

#include "vrp.h"

using namespace std;
using namespace lemon;

VRP::VRP(string inputMapName)   :
    maxspeed{map},
    length{map},
    lat{map},
    lon{map}
{
    Timer timer(1);
    cout << "Reading the map... " << flush;
    digraphReader(map, string(inputMapName))
            .arcMap("maxspeed",maxspeed)
            .arcMap("length",length)
            .nodeMap("lat",lat)
            .nodeMap("lon",lon)
            .run();
    cout << "Elapsed: " << timer.realTime() << "s" << endl;
    cout << "Number of nodes: " << countNodes(map) << endl;
    cout << "Number of arcs: " << countArcs(map) << endl;
}
