//
// Created by david on 2020. 12. 22..
//

#include "vrp.h"

using namespace std;
using namespace lemon;

///Prints errorType and aborts program running if bo is false
void myAssert(bool bo, const string& errorType)
{
    if(!bo) {
        cout << "Assertion error: " << errorType << endl;
        exit(-1);
    }
}

CoordMap::CoordMap(ListDigraph::NodeMap<double>& _lon,
    ListDigraph::NodeMap<double>& _lat) : lon(_lon), lat(_lat)
{
}

CoordMap::Value CoordMap::operator[](const Key& node) const {
    const double Phi = (lat[node] /* + 0.9430/60 */) / 180 * lemon::PI;
    const double Lam = (lon[node] /* + 4.0495/60 */) / 180 * lemon::PI;

    const double r = 6379743.001;
    const double k = 1.003110007693;
    const double n = 1.000719704936;
    const double eps = 0.0818205679407;
    const double Lam0 = (19 + 2.0 / 60 + 54.8584 / 3600) / 180 * lemon::PI;
    const double m0 = 0.9996;
    const double phi0 = 47.1 / 180 * lemon::PI;

    const double sinPhi = std::sin(Phi);

    const double phi
            = 2 * atan(k * std::pow(tan(lemon::PI_4 + Phi / 2), n)
                       * std::pow((1 - eps * sinPhi) / (1 + eps * sinPhi), n * eps / 2)) - lemon::PI_2;

    const double lam = n * (Lam - Lam0);

    const double phi2 = std::asin(std::cos(phi0) * std::sin(phi)
                                  - std::sin(phi0) * std::cos(phi) * std::cos(lam));
    const double lam2 = std::asin(std::cos(phi) * std::sin(lam) / std::cos(phi2));

    const double y = r * m0 * lam2 + 650000;
    const double x = r * m0 * std::log(std::tan(lemon::PI_4 + phi2 / 2)) + 200000;

    return lemon::dim2::Point<double>(y, x);
}

///Calculates haversine distance between two map nodes from longitudes and langitudes
double haversineDist(double lat1, double lon1, double lat2, double lon2){
    lat1 = lat1 / 180 * lemon::PI;
    lon1 = lon1 / 180 * lemon::PI;
    lat2 = lat2 / 180 * lemon::PI;
    lon2 = lon2 / 180 * lemon::PI;
    const double dlon = lon2 - lon1;
    const double dlat = lat2 - lat1;
    return 6371000.0*2 * asin(sqrt(pow(sin(dlat/2), 2) + cos(lat1)*cos(lat2) *
                                                         pow(sin(dlon/2), 2)));
}

///Reads the map from an lgf file and creates the map graph
VRP::VRP(const string& inputName)   :
    maxspeed{map},
    length{map},
    lat{map},
    lon{map},
    coords{lon, lat},
    c{g},
    t{g},
    q{g}
{
    Timer timer(true);
    if (PRINT) cout << "Reading the map... " << flush;
    digraphReader(map, string(inputName))
            .arcMap("maxspeed", maxspeed)
            .arcMap("length", length)
            .nodeMap("lat", lat)
            .nodeMap("lon", lon)
            .run();
    if (PRINT) cout << "Elapsed: " << timer.realTime() << "s" << endl;
    mapNodesNumber = countNodes(map);
    mapArcsNumber = countArcs(map);
    if (PRINT) cout << "Number of nodes: " << mapNodesNumber << endl;
    if (PRINT) cout << "Number of arcs: " << mapArcsNumber << endl << endl;
}

///Initializes the problem: sets init data, generates costumer graph, finds shortest paths
void VRP::init(int initCostumerCnt, int initSeed, int initMaxWeight,
               int initVehicleCapacity)
{
    if(0<initCostumerCnt
        && 0<initMaxWeight && initMaxWeight<=initVehicleCapacity
        && 0<initVehicleCapacity)
    {
        n = initCostumerCnt + 1;
        seed = initSeed;
        maxWeight = initMaxWeight;
        Q = initVehicleCapacity;
        isInitialized=true;

        depotAndCostumers.resize(0);
        g.clear();
        nodes.resize(0);
        arcs.resize(0);
        paths.resize(0);


        generateCostumersGraph();
        if (PRINT) printCostumerCoordinates();
        shortestPaths();
    } else {
        cerr << "Initialization failed!" << endl;
        if(!(0<initCostumerCnt)){
            cerr << "Nonpozitive costumer cnt" << endl;
        }
        if(!(0<initMaxWeight)) {
            cerr << "Nonpozitive max costumer demand" << endl;
        }
        if(!(initMaxWeight<=initVehicleCapacity)){
            cerr << "Max costumer demand is bigger then vehicle capacity (";
            cerr << initVehicleCapacity << ")" << endl;
        }
        if(!(0<initVehicleCapacity)){
            cerr << "Nonpozitive vehicle capacity" << endl;
        }
    }
}

///Generates random costumer points on the map and their demands and creates the costumer graph (g)
void VRP::generateCostumersGraph()
{
    Timer timer(true);
    if(PRINT) cout << "Generate costumer graph.. " << flush;
    depotAndCostumers.reserve(n);
    nodes.reserve(n);
    depotAndCostumers.push_back(0);
    nodes.push_back(g.addNode());

    Random random(seed);
    for(int i = 1; i < n; ++i){
        depotAndCostumers.push_back(random[mapNodesNumber]);
        nodes.push_back(g.addNode());
    }
    arcs.resize(n);
    for(int i = 0; i < n; ++i){
        arcs[i].resize(n);
        for(int j = 0; j < n; ++j){
            if(i != j) {
                arcs[i][j]=g.addArc(nodes[i], nodes[j]);
            }
        }
    }
    if(PRINT) cout << "Generated demands: " << endl;
    for(ListDigraph::NodeIt node(g); node != INVALID; ++node){
        if(g.id(node)==0){
            q[node]=0;
        } else {
            q[node] = random[maxWeight]+1;
        }
        if(PRINT) cout << q[node] << " ";
    }
    if(PRINT) cout << endl;
    if(PRINT) cout << "Elapsed: " << timer.realTime() << "s" << endl << endl;
}

///An arcmap class, calculates the travel time (min)
class ArcTravelTime{
private:
    const ListDigraph::ArcMap<int>& maxspeed;
    const ListDigraph::ArcMap<int>& length;
public:
    typedef double Value;
    typedef const ListDigraph::Arc Key;

    ArcTravelTime(const ListDigraph::ArcMap<int>& in_maxspeed,
                  const ListDigraph::ArcMap<int>& in_length)   :
        maxspeed{in_maxspeed},
        length{in_length}
    {
    }

    Value operator[](Key arc) const {
        return length[arc]/(maxspeed[arc]*1000.0/60.0);
    }
};

///Runs Dijkstras between all costumer nodes and depo and calculates travel times and costs for g's arcs
void VRP::shortestPaths()
{
    Timer timer(true);
    if(PRINT) cout << "Shortest paths..." << flush;
    ArcTravelTime travelTime(maxspeed, length);
    paths.resize(n);
    for(int i=0; i<n; ++i){
        paths[i].resize(n);

        ListDigraph::Node startNode=map.nodeFromId(depotAndCostumers[i]);
        Dijkstra<ListDigraph, ArcTravelTime> dijkstra(map, travelTime);
        dijkstra.run(startNode);

        ListDigraph::Node currNode=INVALID;
        ListDigraph::Arc currArc=INVALID;

        for(int j = 0; j < n; ++j){
            if(j != i){
                currNode=map.nodeFromId(depotAndCostumers[j]);
                c[arcs[i][j]]=0;
                t[arcs[i][j]]=0;
                while(currNode != startNode) {
                    currArc = dijkstra.predArc(currNode);
                    c[arcs[i][j]] += length[currArc];
                    t[arcs[i][j]] += travelTime[currArc];
                    paths[i][j].push_back(currArc);
                    currNode = dijkstra.predNode(currNode);
                }
            }
        }
    }
    if(PRINT) cout << " Elapsed: " << timer.realTime() << "s" << endl << endl;
}

///Prints all costumers, their latitudes, longitudes and demands
void VRP::printCostumerCoordinates()
{
    if(isInitialized) {
        cout << "Costumer coordinates: " << endl;
        for (int i = 0; i < n; ++i) {
            ListDigraph::Node node;
            node = map.nodeFromId(depotAndCostumers[i]);
            cout << i << " : (" << lat[node] << ", " << lon[node] << ")";
            cout << "\t demand: " << q[nodes[i]] << endl;
        }
        cout << endl;
    } else {
        cerr << "Not initialized!" << endl;
    }
}

///Prints the travelling cost between two nodes of g
void VRP::printCost(int sourceId, int targetId){
    if(isInitialized) {
        if (0 <= sourceId && sourceId < n && 0 <= targetId && targetId < n) {
            if (sourceId != targetId) {
                ListDigraph::Node node;
                node = map.nodeFromId(depotAndCostumers[sourceId]);
                cout << "Source: ";
                cout << sourceId << "  (" << lat[node] << ", " << lon[node] << ")" << endl;
                node = map.nodeFromId(depotAndCostumers[targetId]);
                cout << "Target: ";
                cout << targetId << "  (" << lat[node] << ", " << lon[node] << ")" << endl;
                cout << "Cost: " << c[arcs[sourceId][targetId]] << endl;
            } else {
                cerr << "Identical source and target nodes." << endl;
            }
        } else {
            cerr << "Invalid node ID." << endl;
        }
    } else {
        cerr << "Not initialized!" << endl;
    }
}

///Finds the node in the map, which is the closest to the given latitude and longitude
ListDigraph::Node VRP::nodeFromLatLon(double latitude, double longitude)
{
    double minDist=BIG_VALUE;
    ListDigraph::Node closestNode=INVALID;
    double currDist;
    for(ListDigraph::NodeIt node(map); node !=INVALID; ++node){
        currDist=haversineDist(lat[node], lon[node], latitude, longitude);
        if(currDist < minDist){
            minDist=currDist;
            closestNode=node;
        }
    }
    if(PRINT) cout << "Closest node: " << g.id(closestNode) << endl;
    if(PRINT) cout << "Distance: " << minDist << endl;
    return closestNode;
}

///Solves the VRP with the CPLEX MIP, extra conditions on arcs can be added
void VRP::checkMIP(bool printEps, bool printSolutionArcs,
                   vector<tuple<int, int, int>> conditions)
{
    if(isInitialized) {
        Timer timer(true);
        cout << "Check MIP started: " << endl;

        Mip mip;
        //Add cols
        ListDigraph::ArcMap<Mip::Col> checkCols(g);
        for (ListDigraph::ArcIt arc(g); arc != INVALID; ++arc) {
            checkCols[arc] = mip.addCol();
            mip.colLowerBound(checkCols[arc], 0);
            mip.colType(checkCols[arc], Mip::INTEGER);
        }
        ListDigraph::NodeMap<Lp::Col> checkWeights(g);
        for (ListDigraph::NodeIt node(g); node != INVALID; ++node) {
            checkWeights[node] = mip.addCol();
            mip.colLowerBound(checkWeights[node], 0);
            mip.colUpperBound(checkWeights[node], Q);
        }

        //Add rows
        for (ListDigraph::NodeIt node(g); node != INVALID; ++node) {
            if (g.id(node) != 0) {
                Mip::Expr eOut, eIn;
                for (ListDigraph::OutArcIt arc(g, node); arc != INVALID; ++arc) {
                    eOut += checkCols[arc];
                }
                mip.addRow(eOut == 1);
                for (ListDigraph::InArcIt arc(g, node); arc != INVALID; ++arc) {
                    eIn += checkCols[arc];
                }
                mip.addRow(eIn == 1);
            }
        }
        mip.colUpperBound(checkWeights[g.nodeFromId(0)], 0);
        for (ListDigraph::NodeIt node(g); node != INVALID; ++node) {
            if (g.id(node) != 0) {
                for (ListDigraph::InArcIt arc(g, node); arc != INVALID; ++arc) {
                    Lp::Expr e;
                    e = (checkWeights[g.source(arc)] + q[node]);
                    mip.addRow(checkWeights[node] <= BIG_VALUE * (1 - checkCols[arc]) + e);
                    mip.addRow(checkWeights[node] + BIG_VALUE * (1 - checkCols[arc]) >= e);
                }
            }
        }
        for (ListDigraph::ArcIt arc(g); arc != INVALID; ++arc) {
            mip.objCoeff(checkCols[arc], c[arc]);
        }
        for(unsigned int i=0; i<conditions.size(); ++i){
            int sId=get<0>(conditions[i]);
            int tId=get<1>(conditions[i]);
            int value=get<2>(conditions[i]);
            if(0<=sId && sId<n && 0<=tId && tId<=n && sId!=tId && (value==1 || value==0)) {
                mip.addRow(checkCols[arcs[sId][tId]] == value);
            }
        }
        mip.min();
        mip.solve();
        switch (mip.type()) {
            case 0:
                cout << "undefined" << endl;
                break;
            case 1:
                cout << "infeasible" << endl;
                break;
            case 2:
                cout << "feasible" << endl;
                break;
            case 3:
                cout << "optimal" << endl;
                break;
            case 4:
                cout << "unbounded" << endl;
                break;
        }
        cout << "CHECK MIP primal: " << mip.solValue() << endl;
        cout << "Elapsed: " << timer.realTime() << "s" << endl;
        cout << "User time: " << timer.userTime() << "s" << endl;
        if(printSolutionArcs){
            for(ListDigraph::ArcIt arc(g); arc !=INVALID; ++arc){
                if(mip.sol(checkCols[arc])>EPSILON){
                    myAssert(mip.sol(checkCols[arc])>1.0-EPSILON, "check sol error");
                    myAssert(mip.sol(checkCols[arc])<1.0+EPSILON, "check sol error");
                    cout << g.id(g.source(arc)) << "->" << g.id(g.target(arc)) << endl;
                }
            }
            cout << endl;
        }
        if (printEps) {
            printToEpsCheckMIP("checkMIP.eps", mip, checkCols);
        }
    } else {
        cerr << "Not initialized!" << endl;
    }
}

void VRP::printToEpsCheckMIP(const string& filename, const Mip& mip,
                             const ListDigraph::ArcMap<Mip::Col>& checkCols)
{
    Timer timer(true);
    cout << "Printing to eps..." << flush;
    ListDigraph::ArcMap<Color> arcColor(map, Color(0, 0, 0));
    ListDigraph::ArcMap<double> arcWidth(map, 0.0);
    for(ListDigraph::ArcIt arc(g); arc!=INVALID; ++arc){
        if(1-EPSILON<=mip.sol(checkCols[arc])
           && mip.sol(checkCols[arc])<=1+EPSILON)
        {
            //Print this path
            int s=g.id(g.source(arc));
            int t=g.id(g.target(arc));
            for(unsigned int k=0; k<paths[s][t].size(); ++k){
                arcWidth[paths[s][t][k]]=0.1;
                arcColor[paths[s][t][k]]=Color(255, 0, 0);
            }
        }
    }
    ListDigraph::NodeMap<double> nodeSize(map, 0.0);
    ListDigraph::NodeMap<int> nodeShape(map, 0);
    ListDigraph::NodeMap<Color> nodeColor(map, Color(0, 0, 255));
    ListDigraph::Node node=INVALID;
    for(int i = 1; i < n; ++i) {
        node=map.nodeFromId(depotAndCostumers[i]);
        nodeSize[node] = 0.2;
    }
    nodeShape[map.nodeFromId(depotAndCostumers[0])]=1;
    nodeSize[map.nodeFromId(depotAndCostumers[0])]=0.5;

    graphToEps(map, filename)
            .coords(coords)
            .arcWidths(arcWidth)
            .arcColors(arcColor)
            .scale(200)
            .nodeSizes(nodeSize)
            .nodeShapes(nodeShape)
            .nodeColors(nodeColor)
            .run();
    cout << " Elapsed: " << timer.realTime() << "s" << endl << endl;
}

///Creates a BranchAndPrice class, runs the algorithm and saves the solution
void VRP::callBranchAndPrice(double smoothingParameter){
    if(isInitialized) {
        BranchAndPrice bap(g, n, Q, arcs, c, q, smoothingParameter);
        bap.branchAndBound();
        bap.saveSolution(solution, solutionCost);
        if(smoothingParameter==-1) {
            bap.printAlphas();
        }
        cerr << "Save done" << endl;
    } else {
        cerr << "Not initialized!" << endl;
    }
    cerr << "After if" << endl;
}

void VRP::printToEps(const string& filename){
    if(isInitialized) {
        Timer timer(true);
        cout << "Printing to eps..." << flush;
        ListDigraph::ArcMap<Color> arcColor(map, Color(
                200, 200, 200));
        ListDigraph::ArcMap<double> arcWidth(map,
                                             0.030);
        for (unsigned int i = 0; i < solution.size(); ++i) {
            for (unsigned int j = 1; j < solution[i].size(); ++j) {
                int sId = g.id(solution[i][j - 1]);
                int tId = g.id(solution[i][j]);
                for (unsigned int a = 0; a < paths[sId][tId].size(); ++a) {
                    ListDigraph::Arc currArc;
                    currArc = paths[sId][tId][a];
                    arcWidth[currArc] = 0.45;
                    arcColor[currArc] = Color(255, 0, 0);
                }
            }
        }
        ListDigraph::NodeMap<double> nodeSize(map, 0.0);
        ListDigraph::NodeMap<int> nodeShape(map, 0);
        ListDigraph::NodeMap<Color> nodeColor(map, Color(0, 0, 255));
        ListDigraph::Node node = INVALID;
        for (int i = 1; i < n; ++i) {
            node = map.nodeFromId(depotAndCostumers[i]);
            nodeSize[node] = 0.4;
        }
        nodeShape[map.nodeFromId(depotAndCostumers[0])] = 1;
        nodeSize[map.nodeFromId(depotAndCostumers[0])] = 0.6;

        graphToEps(map, filename)
                .coords(coords)
                .arcWidths(arcWidth)
                .arcColors(arcColor)
                .scale(200)
                .nodeSizes(nodeSize)
                .nodeShapes(nodeShape)
                .nodeColors(nodeColor)
                .run();
        cout << " Elapsed: " << timer.realTime() << "s" << endl << endl;
    } else {
        cerr << "Not initialized'" << endl;
    }
}

void VRP::printSolution()
{
    if(solutionCost<BIG_VALUE-1){
        cout << "Solution cost: " << solutionCost << endl;
        for(unsigned int i=0; i<solution.size(); ++i) {
            for (unsigned int j = 0; j < solution[i].size(); ++j) {
                cout << g.id(solution[i][j]) << " ";
            }
            cout << endl;
        }
    } else {
        cerr << "The solution is empty!" << endl;
    }
}

BranchAndPrice::BranchAndPrice(
        ListDigraph& in_g,
        const int& in_n,
        const int& in_Q,
        vector<vector<ListDigraph::Arc>>& in_arcs,
        ListDigraph::ArcMap<int>& in_c,
        ListDigraph::NodeMap<int>& in_q,
        double in_smoothingParameter
)
        : g{in_g},
          n{in_n},
          Q{in_Q},
          arcs{in_arcs},
          c{in_c},
          q{in_q},
          startCols{g},
          nodeRows{g},
          smoothingParameter{in_smoothingParameter},
          ySep{g},
          yIn{g},
          subgradient{g},
          arcUse{g},
          bestSolutionStartCols{g}
{
}

void BranchAndPrice::createMasterLP()
{
    createdMasterLP = true;
    //masterLP.clear();
    cols.resize(0);
    routeNodes.resize(0);

    //Add cols
    for (ListDigraph::NodeIt node(g); node != INVALID; ++node) {
        if (g.id(node) != 0) {
            startCols[node] = masterLP.addCol();
            masterLP.colLowerBound(startCols[node], 0);
            masterLP.colUpperBound(startCols[node], 1);
            //Calculate startCol costs
            masterLP.objCoeff(
                    startCols[node],
                    c[arcs[0][g.id(node)]] + c[arcs[g.id(node)][0]]);
        }
    }
    vehicleNumberCol = masterLP.addCol();
    masterLP.objCoeff(vehicleNumberCol, 0);
    masterLP.colLowerBound(vehicleNumberCol, 0);

    //Add rows
    for (ListDigraph::NodeIt node(g); node != INVALID; ++node) {
        if (g.id(node) != 0) {
            nodeRows[node] = masterLP.addRow(startCols[node] == 1);
        }
    }
    vehicleNumberRow = masterLP.addRow();

    //Set vehicleNumber and totalCost rows coeffs
    for (ListDigraph::NodeIt node(g); node != INVALID; ++node) {
        if (g.id(node) != 0) {
            masterLP.coeff(vehicleNumberRow, startCols[node], 1);
        }
    }
    masterLP.coeff(vehicleNumberRow, vehicleNumberCol, -1);
    masterLP.rowLowerBound(vehicleNumberRow, 0);
    masterLP.rowUpperBound(vehicleNumberRow, 0);
    masterLP.min();
}

class MarginalCost {
private:
    const BranchAndPrice &bap;
public:
    typedef const ListDigraph::Arc &Key;
    typedef double Value;

    MarginalCost(const BranchAndPrice &in_bap) :
            bap{in_bap}
    {}

    Value operator[](MarginalCost::Key arc) {
        ListDigraph::Node sNode, tNode;
        sNode = bap.g.source(arc);
        tNode = bap.g.target(arc);
        double pi_s;
        if (bap.g.id(sNode) == 0) {
            pi_s = bap.ySepVehicle;
        } else {
            pi_s = bap.ySep[sNode];
        }
        double cost = 0;
        cost = bap.c[
                bap.arcs[bap.g.id(sNode)][bap.g.id(tNode)]
        ] - pi_s;
        return cost;
    }
};

bool BranchAndPrice::solveMasterLPwithSmoothing(){
    Timer timer(true);
    if (PRINT) cout << "Solving the master problem.. " << endl;
    int itCnt = 0;
    alpha=0;
    bool continueIteration=true;
    do {
        ++itCnt;
        if (cols.size() % 1000 == 0 && cols.size() > 0) {
            if (PRINT) cout << "Added cols number: " << cols.size() << endl;
        }
        masterLP.solve();
        if(COLUMN_PRINT) printMasterLPSolution();
        switch (masterLP.primalType()) {
            case 0:
                cout << "undefined" << endl;
                return false;
            case 1:
                cout << "infeasible" << endl;
                return false;
            case 2:
                cout << "feasible" << endl;
                break;
            case 3:
                //cout << "optimal" << endl;
                break;
            case 4:
                cout << "unbounded" << endl;
                return false;
        }
        for(ListDigraph::NodeIt node(g); node!=INVALID; ++node){
            ySep[node]=masterLP.dual(nodeRows[node])*(1-alpha)+yIn[node]*alpha;
        }
        ySepVehicle = masterLP.dual(vehicleNumberRow)*(1-alpha)+yInVehicle*alpha;

        continueIteration = generateColumn();
        if(!continueIteration){
            continueIteration = misPricingSchedule();
        }

        //calculate szubgradient dotproduct yOut-yIn
        double dotproduct=0;
        for(ListDigraph::NodeIt node(g); node!=INVALID; ++node){
            dotproduct+=subgradient[node]*(masterLP.dual(nodeRows[node])-yIn[node]);
        }
        dotproduct+=subgradientVehicle*(masterLP.dual(vehicleNumberRow)-yInVehicle);


        if(smoothingParameter==-1) {
            if (itCnt == 1) {
                alpha = 0.1;
                cout << "AUTO SMOOTHING own" << endl;
            } else if (dotproduct < 0) {
                alpha += (1 - alpha) * 0.1;
            } else {
                alpha = alpha - 0.1;
                if (alpha < EPSILON) {
                    alpha = 0;
                }
            }
            alphas.push_back(alpha);
        } else if(smoothingParameter==-2) {
            if (itCnt == 1) {
                alpha = 0.5;
                cout << "AUTO SMOOTHING article" << endl;
            } else if (dotproduct > 0) {
                alpha += (1 - alpha) * 0.1;
            } else {
                alpha = alpha - 0.1;
                if (alpha < EPSILON) {
                    alpha = 0;
                }
            }
        } else if(itCnt==1) {
            alpha=smoothingParameter;
            cout << "CONSTANT ALPHA = " << alpha << endl;
        }
        for(ListDigraph::NodeIt node(g); node!=INVALID; ++node){
            yIn[node]=ySep[node];
        }
        yInVehicle = ySepVehicle;
    } while (continueIteration);
    if (PRINT) {
        cout << "Master LP solved, elapsed: " << timer.realTime() << "s" << endl;
        cout << "User time: " << timer.userTime() << "s" << endl;
        cout << "Number of added columns: " << cols.size() << endl;
        cout.precision(12);
        cout << "Vehicle number: " << masterLP.primal(vehicleNumberCol) << endl;
        cout << "Total cost: " << masterLP.primal() << endl;
        cout << endl;
    }
    if(COLUMN_PRINT) cout << endl << endl << endl << endl << endl << endl << endl << endl;
    return true;
}

class Label {
public:
    const ListDigraph &g;

    ListDigraph::NodeMap<int> nodeUse;  //0 if not used, i if i-th node in path
    double cost;
    int weight;
    int nodeCnt;

    Label(const ListDigraph &inG) :
            g{inG},
            nodeUse{inG, 0},
            cost{0.0},
            weight{0},
            nodeCnt{0}
    {
    }

    Label(const Label &l) :
            g{l.g},
            nodeUse{g, 0},
            cost{l.cost},
            weight{l.weight},
            nodeCnt{l.nodeCnt}
    {
        for (ListDigraph::NodeIt node(g); node != INVALID; ++node) {
            nodeUse[node] = l.nodeUse[node];
        }
    }

    Label &operator=(const Label &inL) {
        cost = inL.cost;
        weight = inL.weight;
        nodeCnt = inL.nodeCnt;
        for (ListDigraph::NodeIt node(g); node != INVALID; ++node) {
            nodeUse[node] = inL.nodeUse[node];
        }
        return *this;
    }
};

class NodeLabels {
public:
    vector<Label> labelList;
    unsigned int nextIndex{0};

    bool next(Label &l) {
        if (nextIndex == labelList.size()) {
            return false;
        } else {
            l = labelList[nextIndex];
            ++nextIndex;
            return true;
        }
    }

    bool dominated(const Label &l) {
        bool isDominator;
        for(int i=0; i<static_cast<int>(labelList.size()); ++i){
            isDominator=true;
            if(l.cost>=labelList[i].cost
                && l.weight>=labelList[i].weight)
            {
                for(ListDigraph::NodeIt node(l.g); node != INVALID; ++node){
                    if(l.nodeUse[node]==0 && labelList[i].nodeUse[node]>0){
                        isDominator=false;
                        break;
                    }
                }
            } else {
                isDominator=false;
            }
            if(isDominator){
                return true;
            }
        }
        return false;
    }
};


bool BranchAndPrice::extendLabel( ListDigraph::NodeMap<NodeLabels>& nodeLabels,
         const ListDigraph::Node& node,
         MarginalCost& mc)
{
    ListDigraph::Node otherNode=INVALID;
    Label l(g);
    Label extendedL(g);
    bool extend=false;
    while (nodeLabels[node].next(l)) {
        extend=true;
        for (ListDigraph::OutArcIt arc(g, node); arc != INVALID; ++arc) {
            otherNode = g.target(arc);
            extendedL=l;
            if (extendedL.nodeUse[otherNode] == 0 &&
                extendedL.weight + q[otherNode] <= Q)
            {
                ++extendedL.nodeCnt;
                extendedL.nodeUse[otherNode] = extendedL.nodeCnt;
                extendedL.weight += q[otherNode];
                extendedL.cost += mc[arc];
                if (!nodeLabels[otherNode].dominated(extendedL)) {
                    //cerr << g.id(node) << "  -> " << g.id(otherNode);
                    //cerr << "  : " << extendedL.cost << endl;
                    nodeLabels[otherNode].labelList.push_back(extendedL);
                }
            }
        }
    }
    return extend;
}

bool BranchAndPrice::generateColumn()
{
    if(COLUMN_PRINT) cout << "Generating column" << endl;
    Timer timer(true);
    MarginalCost mc{*this};
    ListDigraph::NodeMap<NodeLabels> nodeLabels(g);
    ListDigraph::Node depot=g.nodeFromId(0);
    nodeLabels[depot].labelList.emplace_back(g);
    int notExtendedCnt=0;
    extendLabel(nodeLabels, depot, mc);
    int nextCheckDepotIndex=1;
    while(notExtendedCnt<n-1) {
        for (ListDigraph::NodeIt node(g); node != INVALID; ++node) {
            if(g.id(node)==0){
                continue;
            }
            if(notExtendedCnt>=n-1){
                break;
            }
            if(extendLabel(nodeLabels, node, mc)){
                notExtendedCnt=0;
                while (nextCheckDepotIndex < static_cast<int>
                (nodeLabels[depot].labelList.size())) {
                    if (nodeLabels[depot].labelList[nextCheckDepotIndex].
                            cost < -EPSILON) {
                        if (COLUMN_PRINT) {
                            cout << "Early exit!" << endl;
                            cout << "Found min cost: "
                                 << nodeLabels[depot].labelList[nextCheckDepotIndex].cost << endl;
                            cout << "Elapsed: " << timer.realTime() << "s" << endl << endl;
                        }

                        Label l = nodeLabels[depot].labelList[nextCheckDepotIndex];

                        //calculate if mis-pricing
                        vector<ListDigraph::Node> currRouteNodes;
                        double c_i;
                        getRouteFromLabel(l, currRouteNodes, c_i);
                        double dotproduct=0; //y_out*a_{i^t}
                        for(unsigned int i=0; i<currRouteNodes.size(); ++i) {
                            if(currRouteNodes[i]!=g.nodeFromId(0)){
                                dotproduct+=masterLP.dual(nodeRows[currRouteNodes[i]]);
                            }
                        }
                        calculateSubgradient(l);

                        if(c_i-dotproduct<-EPSILON){
                            addGeneratedColumn(currRouteNodes, c_i);
                            return true;
                        } else {
                            cout << "Early mis-pricing" << endl;
                            return false;
                        }
                    } else {
                        ++nextCheckDepotIndex;
                    }
                }
            } else {
                ++notExtendedCnt;
            }
        }
    }
    //search min cost path
    int minIndex=-1;
    double minCost=BIG_VALUE;
    for(int i=0; i<static_cast<int>(nodeLabels[depot].labelList.size()); ++i){
        if(nodeLabels[depot].labelList[i].cost<minCost){
            minCost=nodeLabels[depot].labelList[i].cost;
            minIndex=i;
        }
    }
    if(COLUMN_PRINT) cout << "Found min cost: " << minCost << endl;
    if(COLUMN_PRINT) cout << "Elapsed: " << timer.realTime() << "s" << endl << endl;


    Label l = nodeLabels[depot].labelList[minIndex];

    //calculate if mis-pricing
    vector<ListDigraph::Node> currRouteNodes;
    double c_i;
    getRouteFromLabel(l, currRouteNodes, c_i);
    double dotproduct=0; //y_out*a_{i^t}
    for(unsigned int i=0; i<currRouteNodes.size(); ++i) {
        if(currRouteNodes[i]!=g.nodeFromId(0)){
            dotproduct+=masterLP.dual(nodeRows[currRouteNodes[i]]);
        }
    }
    calculateSubgradient(l);

    if(c_i-dotproduct<-EPSILON){
        addGeneratedColumn(currRouteNodes, c_i);
        return true;
    } else {
        cout << "Late mis-pricing" << endl;
        return false;
    }
}

void BranchAndPrice::calculateSubgradient(const Label& l){
    for(ListDigraph::NodeIt node(g); node!=INVALID; ++node){
        subgradient[node]=1;
        if(l.nodeUse[node]>1){
            subgradient[node]-=n;
        }
    }
    subgradientVehicle=1-n;
}

bool BranchAndPrice::misPricingSchedule(){
    cout << "Mis-pricing sequence started!" << endl;
    cout << "Alpha = " << alpha << endl;
    int k=2;
    double tildeAlpha;
    ListDigraph::NodeMap<double> yIIn(g);
    for(ListDigraph::NodeIt node(g); node!=INVALID; ++node){
        yIIn[node] = yIn[node];
    }
    double yIInVehicle=yInVehicle;
    bool generatedValidCuttingPlane=false;
    do{
        tildeAlpha=1-k*(1-alpha);
        if(tildeAlpha<EPSILON){
            tildeAlpha = 0;
        }
        cout << "TildeAlpha: " << tildeAlpha << endl;
        for(ListDigraph::NodeIt node(g); node!=INVALID; ++node){
            ySep[node]=tildeAlpha*yIIn[node]+(1-tildeAlpha)*masterLP.dual(nodeRows[node]);
        }
        ySepVehicle=tildeAlpha*yIInVehicle+(1-tildeAlpha)*masterLP.dual(vehicleNumberRow);

        generatedValidCuttingPlane = generateColumn();
        ++k;
    } while(!generatedValidCuttingPlane && tildeAlpha>EPSILON);

    return generatedValidCuttingPlane;
}

///Calculates route nodes in order and route cost from a label
void BranchAndPrice::getRouteFromLabel(const Label&l,
                                        vector<ListDigraph::Node>& currRouteNodes,
                                        double& cost){
    int nodeCount=l.nodeCnt+1;
    currRouteNodes.resize(nodeCount, INVALID);
    currRouteNodes[0]=g.nodeFromId(0);
    for(ListDigraph::NodeIt node(g); node != INVALID; ++node){
        if(l.nodeUse[node]!=0) {
            currRouteNodes[l.nodeUse[node]] = node;
        }
    }
    cost=0;
    for(int i=0; i<nodeCount; ++i){
        if(i>0){
            cost+=c[arcs[g.id(currRouteNodes[i-1])][g.id(currRouteNodes[i])]];
        }
    }
}

///Adds a new column to the master problem from routeNodes vector and route cost
void BranchAndPrice::addGeneratedColumn(vector<ListDigraph::Node>& currRouteNodes,
                                        double& cost)
{
    Lp::Col col;
    col=masterLP.addCol();
    cols.push_back(col);
    routeNodes.push_back(currRouteNodes);
    masterLP.colLowerBound(col, 0.0);
    if(COLUMN_PRINT) cout << "Added column's nodes:  ";
    for(unsigned int i=0; i<currRouteNodes.size(); ++i){
        if(COLUMN_PRINT) cout << g.id(currRouteNodes[i]) << " ";
        if(g.id(currRouteNodes[i]) != 0) {
            masterLP.coeff(nodeRows[currRouteNodes[i]], col, 1.0);
        }
    }
    if(COLUMN_PRINT) cout << endl;
    masterLP.coeff(vehicleNumberRow, col, 1);
    masterLP.objCoeff(col, cost);
}

void BranchAndPrice::printMasterLPSolution()
{
    cout << "Master LP solution:" << endl;
    cout << "Number of added columns: " << cols.size() << endl;
    for(unsigned int i=0; i<cols.size(); ++i){
        if(masterLP.primal(cols[i])>EPSILON) {
            cout << "value: " << masterLP.primal(cols[i]) << ", cost: ";
            cout << masterLP.objCoeff(cols[i]) << endl;
            for(unsigned int j=0; j<routeNodes[i].size(); ++j){
                cout << g.id(routeNodes[i][j]) << " ";
            }
            cout << endl;
        }
    }
    for(ListDigraph::NodeIt node(g); node !=INVALID; ++node) {
        if(g.id(node)!=0) {
            cout << g.id(node) << " . node: ";
            cout << masterLP.primal(startCols[node]) << endl;
        }
    }
    cout.precision(12);
    cout << "Vehicle number: " << masterLP.primal(vehicleNumberCol) << endl;
    cout << "Total cost: " << masterLP.primal() << endl;


    //Calculate dual cost
    double dualCost=0;
    for(Lp::RowIt row(masterLP); row!=INVALID; ++row){
        dualCost+=masterLP.dual(row);
    }
    dualCost-=masterLP.dual(vehicleNumberRow);
    //dualCost-=masterLP.dual(totalCostRow);
    cout << "Dual cost: " << dualCost << endl;
    cout << endl;
}

void BranchAndPrice::branchAndBound()
{
    Timer timer(true);
    bestCost = BIG_VALUE;
    if (PRINT) cout << "Branch and bound started: " << endl;
    int branchedNodes = 0;
    createMasterLP();
    recursiveBranch(branchedNodes);

    if (PRINT) {
        cout << endl;
        cout << "Finished" << endl;
        //printBranchedArcs();
    }
    //cout << n - 1 << " " << seed << " " << maxWeight << endl;
    myAssert(countArcs(g) == n * (n - 1), "Wrong arc count at the end of branchAndBound!");
    cout << "Elapsed real time: " << timer.realTime() << "s" << endl;
    cout << "Elapsed user time: " << timer.userTime() << "s" << endl;
    cout << "Solution: " << bestCost << endl;
    cout << "Vehicle count: " << bestSolutionVehicle << endl;
    cout << "Number of added cols: " << cols.size() << endl;
    cout << "Visited branching nodes: " << branchedNodes << endl;
}

bool isWhole(double x){
    if(x-floor(x)<EPSILON || ceil(x)-x<EPSILON){
        return true;
    } else {
        return false;
    }
}

void BranchAndPrice::recursiveBranch(int& branchedNodes)
{
    ++branchedNodes;
    if(PRINT) cout << "New branching node: " << branchedNodes << endl;
    if(!solveMasterLPwithSmoothing()){
        return;
    } else if(masterLP.primal()>bestCost) {
        if(PRINT) cout << "BOUND" << endl;
        return;
    }

    double vehicleNumber=masterLP.primal(vehicleNumberCol);
    Lp::Row tempRow=INVALID;
    if(!isWhole(vehicleNumber)){
        //branch on vehicle number
        tempRow=masterLP.addRow(vehicleNumberCol<=floor(vehicleNumber));
        if(PRINT) cout << "Add: vehicleNumberCol <= " << floor(vehicleNumber) << endl;
        recursiveBranch(branchedNodes);
        if(PRINT) cout << "Remove vehicleNumberCol <= " << floor(vehicleNumber) << endl;
        masterLP.erase(tempRow);

        tempRow=masterLP.addRow(vehicleNumberCol>=ceil(vehicleNumber));
        if(PRINT) cout << "Add: vehicleNumberCol >= " << ceil(vehicleNumber) << endl;
        recursiveBranch(branchedNodes);
        if(PRINT) cout << "Remove: vehicleNumberCol >= " << ceil(vehicleNumber) << endl;
        masterLP.erase(tempRow);
    } else {
        //branch on arcs
        calculateArcUse();
        double maxFractionalUse=0.0;
        ListDigraph::Arc arcToBranch=INVALID;
        for(int i=0; i<n; ++i){
            for(int j=0; j<n; ++j) {
                if (i != j && arcs[i][j] != INVALID) {
                    ListDigraph::Arc arc;
                    arc = arcs[i][j];
                    if (!isWhole(arcUse[arc])) {
                        if (COLUMN_PRINT) {
                            cout << g.id(g.source(arc)) << " " << g.id(g.target(arc)) << " " << arcUse[arc] << "  "
                                 << endl;
                        } else if (PRINT) {
                            cout << arcUse[arc] << "  ";
                        }
                        if (arcUse[arc] > maxFractionalUse + EPSILON && arcUse[arc] < 1) {
                            maxFractionalUse = arcUse[arc];
                            arcToBranch = arc;
                        }
                    }
                }
            }
        }
/*
        for(ListDigraph::ArcIt arc(g); arc!=INVALID; ++arc){
            if(!isWhole(arcUse[arc])){
                if(COLUMN_PRINT){
                    cout << g.id(g.source(arc)) << " " << g.id(g.target(arc)) << " " << arcUse[arc] << "  " << endl;
                } else {
                    cout << arcUse[arc] << "  ";
                }
                if(arcUse[arc]>maxFractionalUse+EPSILON && arcUse[arc]<1){
                    maxFractionalUse=arcUse[arc];
                    arcToBranch=arc;
                }
                //If arcUse is the same
                //Choose according to alphabetical order to avoid random branching arc
                if(abs(arcUse[arc]-maxFractionalUse)<EPSILON){
                    if(g.id(g.source(arcToBranch))>g.id(g.source(arc))){
                        arcToBranch=arc;
                    } else if(g.id(g.source(arcToBranch))==g.id(g.source(arc))
                        && g.id(g.target(arcToBranch))>g.id(g.target(arc))){
                        arcToBranch=arc;
                    }
                }
            }
        }
*/
        cout  << endl;
        if(arcToBranch==INVALID){
            cout.precision(12);
            if(PRINT) cout << "WHOLE solution found with cost: ";
            if(PRINT) cout << masterLP.primal();

            if(masterLP.primal()<bestCost){
                //This solution is better than the current best
                //Save solution
                if(PRINT) cout << " Better" << endl;
                bestCost=masterLP.primal();
                bestSolutionVehicle=masterLP.primal(vehicleNumberCol);
                bestSolutionColIndexs.resize(0);
                for(unsigned int i=0; i<cols.size(); ++i){
                    if(masterLP.primal(cols[i])>EPSILON) {
                        if(PRINT) {
                            cout << i << " .col,";
                            cout << " value: ";
                            cout << masterLP.primal(cols[i]) << endl;
                            cout << "Route: ";
                            //printRoutes(i);
                            cout << "Route cost: " << masterLP.objCoeff(cols[i]) << endl;
                            cout << endl;
                        }
                    }
                }
                for(ListDigraph::NodeIt node(g); node!=INVALID; ++node){
                    if(g.id(node)!=0)// &&
                        //masterLP.primal(startCols[node])>EPSILON)
                    {
                        if(PRINT) {
                            cout << g.id(node) << " . node col,";
                            cout << " value: ";
                            cout << masterLP.primal(startCols[node]) << endl;
                            cout << "Route cost: ";
                            cout << masterLP.objCoeff(startCols[node]) << endl;
                            cout << endl;
                        }
                    }
                }
                //printRoutes(-1);
                for(unsigned int i=0; i<cols.size(); ++i){
                    if(masterLP.primal(cols[i])>EPSILON){
                        myAssert(masterLP.primal(cols[i])<1.0+EPSILON, ">1+eps col");
                        myAssert(masterLP.primal(cols[i])>1.0-EPSILON, "<1-eps nonzero col");
                        bestSolutionColIndexs.push_back(i);
                    }
                }
                for(ListDigraph::NodeIt node(g); node!=INVALID; ++node){
                    myAssert(masterLP.primal(startCols[node])<1.0+EPSILON, ">1+eps col");
                    myAssert(masterLP.primal(startCols[node])>1.0-EPSILON
                        || masterLP.primal(startCols[node])<EPSILON, "<1-eps nonzero col");
                    bestSolutionStartCols[node]=masterLP.primal(startCols[node]);
                }
            } else {
                if(PRINT) cout << " Worse" << endl;
            }
        } else {
            //Branch on this arc
            if(PRINT) cout << "Branch on arc with use: ";
            if(PRINT) cout << arcUse[arcToBranch] << endl;
            //branchedArcs.push_back(arcToBranch);
            if(PRINT) cout << "Source: " << g.id(g.source(arcToBranch));
            if(PRINT) cout << ", target: " << g.id(g.target(arcToBranch)) << endl;

            int originalArcCount=countArcs(g); //save original arc count
            ListDigraph::Node sourceNode=g.source(arcToBranch);
            ListDigraph::Node targetNode=g.target(arcToBranch);
            vector<pair<int, int>> changedCosts;    // <index, oldCost>
            vector<pair<ListDigraph::Node, int>> changedStartCosts;

            //arc=1
            if(PRINT) cout << "Add: arc =1" << endl;
            changedCosts.resize(0);
            changedStartCosts.resize(0);
            vector<tuple<ListDigraph::Node, ListDigraph::Node, int>> erasedArcs;
            if(sourceNode!=g.nodeFromId(0)) {
                for (ListDigraph::OutArcIt arc(g, sourceNode); arc != INVALID;) {
                    ListDigraph::Arc a = arc;
                    ++arc;
                    if (g.target(a) != targetNode) {
                        changeObjCoeffs(a, changedCosts, changedStartCosts);
                        erasedArcs.push_back(tuple<ListDigraph::Node, ListDigraph::Node, int>
                                                     (g.source(a), g.target(a), c[a]));
                        arcs[g.id(g.source(a))][g.id(g.target(a))]=INVALID;
                        g.erase(a);
                    }
                }
            }
            if(targetNode!=g.nodeFromId(0)) {
                for (ListDigraph::InArcIt arc(g, targetNode); arc != INVALID;) {
                    ListDigraph::Arc a = arc;
                    ++arc;
                    if (g.source(a) != sourceNode) {
                        changeObjCoeffs(a, changedCosts, changedStartCosts);
                        erasedArcs.push_back(tuple<ListDigraph::Node, ListDigraph::Node, int>
                                                     (g.source(a), g.target(a), c[a]));
                        arcs[g.id(g.source(a))][g.id(g.target(a))]=INVALID;
                        g.erase(a);
                    }
                }
            }
            recursiveBranch(branchedNodes);
            for(unsigned int i=0; i<changedCosts.size(); ++i){
                masterLP.objCoeff(cols[changedCosts[i].first], changedCosts[i].second);
            }
            for(unsigned int i=0; i<changedStartCosts.size(); ++i){
                masterLP.objCoeff(startCols[changedStartCosts[i].first], changedStartCosts[i].second);
            }
            for(unsigned int i=0; i<erasedArcs.size(); ++i){
                ListDigraph::Node sNode=std::get<0>(erasedArcs[i]);
                ListDigraph::Node tNode=std::get<1>(erasedArcs[i]);
                arcs[g.id(sNode)][g.id(tNode)]
                        =g.addArc(sNode, tNode);
                c[arcs[g.id(sNode)][g.id(tNode)]]=std::get<2>(erasedArcs[i]);
            }

            //0
            if(PRINT) cout << "Add: arc =0" << endl;
            changedCosts.resize(0);
            changedStartCosts.resize(0);
            int originalArcCost=c[arcToBranch];
            changeObjCoeffs(arcToBranch, changedCosts, changedStartCosts);
            arcs[g.id(g.source(arcToBranch))][g.id(g.target(arcToBranch))]=INVALID;
            g.erase(arcToBranch);
            recursiveBranch(branchedNodes);
            //Set the original costs
            for(unsigned int i=0; i<changedCosts.size(); ++i){
                masterLP.objCoeff(cols[changedCosts[i].first], changedCosts[i].second);
            }
            for(unsigned int i=0; i<changedStartCosts.size(); ++i){
                masterLP.objCoeff(startCols[changedStartCosts[i].first], changedStartCosts[i].second);
            }
            arcs[g.id(sourceNode)][g.id(targetNode)]=g.addArc(sourceNode, targetNode);
            c[arcs[g.id(sourceNode)][g.id(targetNode)]]=originalArcCost;
            myAssert(countArcs(g)==originalArcCount, "Arc(s) vanished!");
            if(PRINT) cout << "Remove: arc =0" << endl;

            myAssert(countArcs(g)==originalArcCount, "Arc(s) vanished!");
            if(PRINT) cout << "Remove: arc =1" << endl;
        }
    }
}

void BranchAndPrice::changeObjCoeffs(ListDigraph::Arc arc, vector<pair<int, int>>& changedCosts,
    vector<pair<ListDigraph::Node, int>>& changedStartCosts)
{
    if(COLUMN_PRINT) cout << "CHANGE OBJ COEFF on arc: " << g.id(g.source(arc)) << "->" << g.id(g.target(arc)) << endl;
    for(unsigned int i=0; i<cols.size(); ++i){
        bool containsArc=false;
        for(unsigned int j=1; j<routeNodes[i].size(); ++j){
            if(routeNodes[i][j-1]==g.source(arc)
                && routeNodes[i][j]==g.target(arc))
            {
                containsArc=true;
                break;
            }
        }
        if(containsArc && masterLP.objCoeff(cols[i])<BIG_VALUE-1)
        {
            changedCosts.push_back(pair<int, int>(i, masterLP.objCoeff(cols[i])));
            masterLP.objCoeff(cols[i], BIG_VALUE);
        }
    }
    if(g.id(g.source(arc))==0){
        ListDigraph::Node node=g.target(arc);
        if(masterLP.objCoeff(startCols[node])<BIG_VALUE-1) {
            changedStartCosts.push_back(
                    pair<ListDigraph::Node, int>(node, masterLP.objCoeff(startCols[node]))
            );
            masterLP.objCoeff(startCols[node], BIG_VALUE);
        }
    }
    if(g.id(g.target(arc))==0){
        ListDigraph::Node node=g.source(arc);
        if(masterLP.objCoeff(startCols[node])<BIG_VALUE-1) {
            changedStartCosts.push_back(
                    pair<ListDigraph::Node, int>(node, masterLP.objCoeff(startCols[node]))
            );
            masterLP.objCoeff(startCols[node], BIG_VALUE);
        }
    }
}

void BranchAndPrice::calculateArcUse(){
    double usage;
    for(ListDigraph::ArcIt arc(g); arc!=INVALID; ++arc){
        arcUse[arc]=0.0;
    }
    for(unsigned int i=0; i<cols.size(); ++i){
        usage=masterLP.primal(cols[i]);
        for(unsigned int nodeIndex=1; nodeIndex<routeNodes[i].size(); ++nodeIndex){
            int sId=g.id(routeNodes[i][nodeIndex-1]);
            int tId=g.id(routeNodes[i][nodeIndex]);
            arcUse[arcs[sId][tId]]+=usage;
        }
    }
    for(ListDigraph::NodeIt node(g); node!=INVALID; ++node) {
        if (g.id(node) != 0) {
            usage = masterLP.primal(startCols[node]);
            arcUse[arcs[0][g.id(node)]] += usage;
            arcUse[arcs[g.id(node)][0]] += usage;
        }
    }
}

void BranchAndPrice::saveSolution(vector<vector<ListDigraph::Node>> &solution,
                                  int &solutionCost)
{
    if(bestCost<BIG_VALUE-1){
        cout << "Save solution!" << endl;
        solutionCost=bestCost;
        solution.resize(0);
        for(unsigned int i=0; i<bestSolutionColIndexs.size(); ++i){
            solution.push_back(routeNodes[bestSolutionColIndexs[i]]);
        }
        for(ListDigraph::NodeIt node(g); node!=INVALID; ++node){
            if(g.id(node)!=0 && bestSolutionStartCols[node]>EPSILON){
                myAssert(1-EPSILON<bestSolutionStartCols[node], "Save solution: fractional route error!");
                myAssert((bestSolutionStartCols[node]<1+EPSILON), "Save solution: fractional route error!");
                vector<ListDigraph::Node> temp;
                temp.push_back(g.nodeFromId(0));
                temp.push_back(node);
                temp.push_back(g.nodeFromId(0));
                solution.push_back(temp);
            }
        }
        cout << "Solution size: " << solution.size() << endl;
        cout << "Solution cost: " << solutionCost << endl << endl;
    } else {
        cerr << "The solution is empty!" << endl;
    }
}

void BranchAndPrice::printAlphas(){
    cout << "ALPHAAS:" << endl;
    for(unsigned int i=0; i<alphas.size(); ++i){
        cout << alphas[i] << "  ";
    }
    cout << endl;
    cout << "END OF ALPHAS!" << endl;
}