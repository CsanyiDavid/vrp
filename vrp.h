//
// Created by david on 2020. 12. 22..
//

#include <iostream>
#include <lemon/dim2.h>
#include <lemon/list_graph.h>
#include <lemon/lgf_reader.h>
#include <lemon/random.h>
#include <lemon/time_measure.h>
#include <string>

#ifndef VRP_H
#define VRP_H

using namespace std;
using namespace lemon;

class CoordMap{
    ListDigraph::NodeMap<double>& lon;
    ListDigraph::NodeMap<double>& lat;
public:
    typedef ListDigraph::Node Key;
    typedef dim2::Point<double> Value;

    CoordMap(ListDigraph::NodeMap<double>& _lon,
             ListDigraph::NodeMap<double>& _lat) : lon(_lon), lat(_lat){
    }

    Value operator[](const Key& node) const;
};

class VRP{
private:
    //The map
    ListDigraph map;
    ListDigraph::ArcMap<int> maxspeed;
    ListDigraph::ArcMap<int> length;
    ListDigraph::NodeMap<double> lat;
    ListDigraph::NodeMap<double> lon;
    int mapNodesNumber;
    int mapArcsNumber;
    CoordMap coords;
    vector<int> depotAndCostumers;

    //The graph
    ListDigraph g;
    int n;
    vector<ListDigraph::Node> nodes;
    //ListDigraph::NodeMap<int> ids;
    vector<vector<ListDigraph::Arc>> arcs;

public:
    VRP(string inputMapName);

    void generateCostumersGraph(int in_n);
};


#endif /* VRP_H */
