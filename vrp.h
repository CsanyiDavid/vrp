//
// Created by david on 2020. 12. 22..
//

#include <iostream>
#include <lemon/list_graph.h>
#include <lemon/lgf_reader.h>
#include <lemon/time_measure.h>
#include <string>

#ifndef VRP_H
#define VRP_H

using namespace std;
using namespace lemon;

class VRP{
private:
    //The map
    ListDigraph map;
    ListDigraph::ArcMap<int> maxspeed;
    ListDigraph::ArcMap<int> length;
    ListDigraph::NodeMap<double> lat;
    ListDigraph::NodeMap<double> lon;

    //The graph
    ListDigraph g;

public:
    VRP(string inputMapName);
};


#endif /* VRP_H */
