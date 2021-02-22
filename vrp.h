//
// Created by david on 2020. 12. 22..
//

#include <algorithm>
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
#define COLUMN_PRINT false

#ifndef VRP_H
#define VRP_H

#define BIG_VALUE 1000000000.0
#define EPSILON 0.000001

using namespace std;
using namespace lemon;

///Prints errorType and aborts program running if bo is false
void myAssert(bool bo, const string& errorType);

///Creates plane coordinates from longitude and latitude
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

///Calculates haversine distance between two map nodes from longitudes and langitudes
double haversineDist(double lat1, double lon1, double lat2, double lon2);

class MarginalCost;
class Label;
class NodeLabels;

/// Solves the VRP with the branch and price algorithm
class BranchAndPrice{
private:
    bool createdMasterLP=false;

    //The problem
    ListDigraph& g;
    int& n;     //costumer+depo count
    int& Q;     //vehicle capacity
    vector<ListDigraph::Node>& nodes;
    vector<vector<ListDigraph::Arc>>& arcs;
    ListDigraph::ArcMap<int>& c;     //travel distance (meters)
    ListDigraph::NodeMap<int>& q;   //costumer demands

    // The Master Problem
    Lp masterLP;
    vector<Lp::Col> cols;
    vector<vector<ListDigraph::Node>> routeNodes;   //added column/route nodes
    ListDigraph::NodeMap<Lp::Col> startCols;
    ListDigraph::NodeMap<Lp::Row> nodeRows;
    Lp::Col vehicleNumberCol;
    Lp::Row vehicleNumberRow;

    vector<ListDigraph::Arc> branchedArcs;

    //Best solution
    ListDigraph::ArcMap<double> arcUse;
    double bestCost=BIG_VALUE;
    int bestSolutionVehicle;
    vector<int> bestSolutionColIndexs;
    ListDigraph::NodeMap<double> bestSolutionStartCols;
public:
    BranchAndPrice(ListDigraph& in_g, int& n, int& in_Q,
                   vector<ListDigraph::Node>& in_nodes,
                   vector<vector<ListDigraph::Arc>>& in_arcs,
                   ListDigraph::ArcMap<int>& in_c, ListDigraph::NodeMap<int>& in_q);
    void branchAndBound();
    void printMasterLPSolution();
    void printMasterLPMatrix();
    void printBranchedArcs();
    void saveSolution(vector<vector<ListDigraph::Node>>& solution, int& solutionCost);

    friend class MarginalCost;
    friend class Label;
private:
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

/*
typedef tuple<ListDigraph::Node, ListDigraph::Node, int> savingType;

///Solves the VRP with the Clarke-Wright algorithm
class ClarkeWright{
private:
    //The problem
    ListDigraph& g;
    int& n;     //costumer+depo count
    int& Q;     //vehicle capacity
    vector<ListDigraph::Node>& nodes;
    vector<vector<ListDigraph::Arc>>& arcs;
    ListDigraph::ArcMap<int>& c;     //travel distance (meters)
    ListDigraph::NodeMap<int>& q;   //costumer demands

    //
    vector<savingType> savings;
public:
    ClarkeWright(ListDigraph& in_g, int& n, int& in_Q,
                 vector<ListDigraph::Node>& in_nodes,
                 vector<vector<ListDigraph::Arc>>& in_arcs,
                 ListDigraph::ArcMap<int>& in_c, ListDigraph::NodeMap<int>& in_q);
    void calculateSavings();
    void printSavings();
    void run();
};
*/

/// Manages and stores the VRP: reads map, finds shortest paths, calls solvers, prints solution
class VRP{
private:
    //Init data
    int seed;
    int maxWeight;
    int Q;  //vehicle capacity
    bool isInitialized=false;

    //The map
    ListDigraph map;
    ListDigraph::ArcMap<int> maxspeed;  //(km/h)
    ListDigraph::ArcMap<int> length;    //(m)
    ListDigraph::NodeMap<double> lat;
    ListDigraph::NodeMap<double> lon;
    int mapNodesNumber;
    int mapArcsNumber;
    CoordMap coords;
    vector<int> depotAndCostumers;  //Which map nodes are costumers or depo

    //The graph
    ListDigraph g;
    int n;      //costumer+depo count
    vector<ListDigraph::Node> nodes;
    vector<vector<ListDigraph::Arc>> arcs;
    vector<vector<vector<ListDigraph::Arc>>> paths;     //arcs of the paths between two important nodes on the map
    ListDigraph::ArcMap<int> c;     //travel distance (meters)
    ListDigraph::ArcMap<double> t;     //travel time (minutes)
    ListDigraph::NodeMap<int> q;    //costumer demands

    //The solution
    vector<vector<ListDigraph::Node>> solution;
    int solutionCost;

public:
    VRP(const string& inputName);
    void init(int initCostumerCnt, int initSeed, int initMaxWeight, int initVehicleCapacity=100);
    void printCostumerCoordinates();
    void printCost(int sourceId, int targetId);
    ListDigraph::Node nodeFromLatLon(double latitude, double longitude);
    void checkMIP(bool printEps=false,
                  bool printSolutionArcs=false,
                  vector<tuple<int, int, int>> conditions=vector<tuple<int, int, int>>(0)
                  );
    void printToEpsCheckMIP(const string &, const Mip& mip,
                            const ListDigraph::ArcMap<Mip::Col>& cols);
    void callBranchAndPrice();
    //void printRoutes(int index=-1);
    void printToEps(const string& filename);
    void printSolution();
    //void callClarkeWright();

private:
    void generateCostumersGraph();
    void shortestPaths();
};

#endif /* VRP_H */
