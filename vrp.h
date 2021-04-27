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
#define TIME_LIMIT 300

#ifndef VRP_H
#define VRP_H

#define BIG_VALUE 1000000000.0
#define EPSILON 0.000001

using namespace std;
using namespace lemon;

void myAssert(bool condition, const string& errorType);

///Creates plane coordinates from longitude and latitude
class CoordMap{
    ListDigraph::NodeMap<double>& lon;
    ListDigraph::NodeMap<double>& lat;
public:
    typedef ListDigraph::Node Key;
    typedef dim2::Point<double> Value;

    CoordMap(ListDigraph::NodeMap<double>& in_lon,
             ListDigraph::NodeMap<double>& in_lat)
            : lon(in_lon), lat(in_lat){}

    Value operator[](const Key& node) const;
};

double haversineDist(double lat1, double lon1, double lat2, double lon2);

class MarginalCost;
class Label;
class NodeLabels;

/// Solves the VRP with the branch and price algorithm
class BranchAndPrice{
private:
    bool createdMasterLP=false;

    //Init data
    ListDigraph& g;
    const int& n;     //costumer+depo count
    const int& Q;     //vehicle capacity
    vector<vector<ListDigraph::Arc>>& arcs;
    ListDigraph::ArcMap<int>& c;     //travel distance (meters)
    ListDigraph::NodeMap<int>& q;   //costumer demands
    bool earlyStop;
    double smoothingParameter;

    // The Master Problem
    Lp masterLP;
    vector<Lp::Col> cols;
    vector<vector<ListDigraph::Node>> routeNodes;   //The nodes of the added columns/routes
    ListDigraph::NodeMap<Lp::Col> startCols;
    ListDigraph::NodeMap<Lp::Row> nodeRows;
    Lp::Col vehicleNumberCol;
    Lp::Row vehicleNumberRow;

    //Smoothing
    ListDigraph::NodeMap<double> ySep;
    ListDigraph::NodeMap<double> yIn;
    ListDigraph::NodeMap<double> subgradient;
    double ySepVehicle;
    double yInVehicle;
    double subgradientVehicle;
    double alpha;
    vector<double> alphas;

    ListDigraph::ArcMap<double> arcUse;

    //Best solution
    double bestCost=BIG_VALUE;
    int bestSolutionVehicle;
    vector<int> bestSolutionColIndexs;
    ListDigraph::NodeMap<double> bestSolutionStartCols;
public:
    BranchAndPrice(ListDigraph& in_g, const int& n, const int& in_Q,
                   vector<vector<ListDigraph::Arc>>& in_arcs,
                   ListDigraph::ArcMap<int>& in_c, ListDigraph::NodeMap<int>& in_q,
                   double in_smoothingParameter, bool in_earlyStop);
    void branchAndBound();
    void printMasterLPSolution();
    void saveSolution(vector<vector<ListDigraph::Node>>& solution, int& solutionCost);

    void printAlphas();

    friend class MarginalCost;
    friend class Label;
private:
    void createMasterLP();
    bool solveMasterLPwithSmoothing();
    bool extendLabel( ListDigraph::NodeMap<NodeLabels>& nodeLabels, const ListDigraph::Node& node,
                      MarginalCost& mc);
    bool generateColumn();
    void calculateSubgradient(const Label& l);
    bool misPricingSchedule();
    void getRouteFromLabel(const Label& l, vector<ListDigraph::Node>& currRouteNodes,
                            double& cost);
    void addGeneratedColumn(vector<ListDigraph::Node>& currRouteNodes,
                            double& cost);
    bool recursiveBranch(int&  branchedNodes, Timer& timer);
    void changeObjCoeffs(ListDigraph::Arc arc, vector<pair<int, int>>& changedCosts,
                         vector<pair<ListDigraph::Node, int>>& changedStartCosts);
    void calculateArcUse();
    void findArcToBranch(ListDigraph::Arc& arcToBranch);
    void dealWithFoundWholeSolution();
};

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
    vector<int> depotAndCostumers;  //Which map nodes are the costumers or the depo

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
    void callBranchAndPrice(double smoothingParameter, bool earlyStop);
    void printToEps(const string& filename);
    void printSolution();

private:
    void generateCostumersGraph();
    void shortestPaths();
};

#endif /* VRP_H */
