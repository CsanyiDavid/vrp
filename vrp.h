//
// Created by david on 2020. 12. 22..
//

#include <iostream>
#include <lemon/color.h>
#include <lemon/dijkstra.h>
#include <lemon/dim2.h>
#include <lemon/graph_to_eps.h>
#include <lemon/list_graph.h>
#include <lemon/lgf_reader.h>
#include <lemon/lp.h>
#include <lemon/random.h>
#include <lemon/time_measure.h>
#include <string>
#include <vector>

#ifndef VRP_H
#define VRP_H

#define BIG_VALUE 1000000000
#define EPSILON 0.000001

using namespace std;
using namespace lemon;

void myAssert(bool bo, const string& errorType);

class CoordMap{
    ListDigraph::NodeMap<double>& lon;
    ListDigraph::NodeMap<double>& lat;
public:
    typedef ListDigraph::Node Key;
    typedef dim2::Point<double> Value;

    CoordMap(ListDigraph::NodeMap<double>& _lon,
             ListDigraph::NodeMap<double>& _lat);

    Value operator[](const Key& node) const;
};

double haversineDist(double lat1, double lon1, double lat2, double lon2);

class MarginalCost;

class Label;

class NodeLabels;

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
    vector<vector<vector<ListDigraph::Arc>>> paths;

    ListDigraph::ArcMap<int> c;     //travel distance (meters)
    ListDigraph::ArcMap<double> t;     //travel time (minutes)
    //ListDigraph::NodeMap<int> a;
    //ListDigraph::NodeMap<int> b;
    ListDigraph::NodeMap<int> q;
    int Q;

    // The Master Problem
    Lp masterLP;
    vector<Lp::Col> cols;
    vector<vector<ListDigraph::Arc>> routes;
    ListDigraph::NodeMap<Lp::Col> startCols;
    ListDigraph::NodeMap<Lp::Row> nodeRows;
    Lp::Col vehicleNumberCol;
    Lp::Col totalCostCol;
    Lp::Row vehicleNumberRow;
    Lp::Row totalCostRow;

    //Check Mip
    Mip lp;

public:
    VRP(bool isMap, const string& inputName);

    void generateCostumersGraph(int costumerCnt);

    void printToEps(const string& filename);

    void shortestPaths();

    void printCostumerCoordinates();

    ListDigraph::Node nodeFromLatLon(double latitude, double longitude);

    void printShortestPathsFromDepot();

    void createMasterLP();

    void printMasterLPSolution();

    void solveMasterLP();

    bool extendLabel( ListDigraph::NodeMap<NodeLabels>& nodeLabels, const ListDigraph::Node& node,
                     MarginalCost& mc);

    bool generateColumn();

    void addGeneratedColumn(const Label& l);

    void printMasterLPMatrix();

    void checkLP();
    void printToEpsCheckLp();

    friend class MarginalCost;
    friend class Label;
};

#endif /* VRP_H */
