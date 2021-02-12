//
// Created by david on 2020. 12. 22..
//

#include <iostream>
#include <lemon/arg_parser.h>
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
#include <tuple>
#include <utility>
#include <vector>

#define PRINT true

#ifndef VRP_H
#define VRP_H

#define BIG_VALUE 1000000000.0
#define EPSILON 0.0001

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
    //Init data
    int seed;
    int maxWeight;
    int Q;  //vehicle capacity
    bool isInitialized=false;
    bool createdMasterLP=false;

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
    vector<vector<ListDigraph::Arc>> arcs;
    vector<vector<vector<ListDigraph::Arc>>> paths;

    ListDigraph::ArcMap<int> c;     //travel distance (meters)
    ListDigraph::ArcMap<double> t;     //travel time (minutes)
    ListDigraph::NodeMap<int> q;

    // The Master Problem
    Lp masterLP;
    vector<Lp::Col> cols;
    vector<vector<ListDigraph::Arc>> routes;
    vector<vector<ListDigraph::Node>> routeNodes;
    ListDigraph::NodeMap<Lp::Col> startCols;
    ListDigraph::NodeMap<Lp::Row> nodeRows;
    Lp::Col vehicleNumberCol;
    Lp::Row vehicleNumberRow;

    //Best solution
    ListDigraph::ArcMap<double> arcUse;
    double bestCost=BIG_VALUE;
    int bestSolutionVehicle;
    vector<double> bestSolutionColIndexs;
    ListDigraph::NodeMap<double> bestSolutionStartCols;

public:
    VRP(bool isMap, const string& inputName);
    void init(int initCostumerCnt, int initSeed, int initMaxWeight, int initVehicleCapacity=100);
    void printToEps(const string& filename);
    void printCostumerCoordinates();
    ListDigraph::Node nodeFromLatLon(double latitude, double longitude);
    void printShortestPathsFromDepot();
    //void printMasterLPSolution();
    void printRoutes(int index=-1);
    void printMasterLPMatrix();
    void checkMIP(bool printEps=false);
    void printToEpsCheckMIP(const string &, const Mip& mip,
                            const ListDigraph::ArcMap<Mip::Col>& cols);
    void branchAndBound();
    void printCost(int sourceId, int targetId);
    friend class MarginalCost;
    friend class Label;

private:
    void generateCostumersGraph();
    void shortestPaths();
    void createMasterLP();
    bool solveMasterLP();
    bool extendLabel( ListDigraph::NodeMap<NodeLabels>& nodeLabels, const ListDigraph::Node& node,
                      MarginalCost& mc);
    bool generateColumn();
    void addGeneratedColumn(const Label& l);
    void recursiveBranch(int&  branchedNodes);
    void changeObjCoeffs(ListDigraph::Arc arc, vector<pair<int, int>>& changedCosts,
                         vector<pair<ListDigraph::Node, int>>& changedStartCosts);
    void calculateArcUse();
};

#endif /* VRP_H */
