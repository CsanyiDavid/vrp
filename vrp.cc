//
// Created by david on 2020. 12. 22..
//

#include "vrp.h"

extern int mySeed;
extern int myMaxWeight;

using namespace std;
using namespace lemon;

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

//Read a Map and create graph if isMap, else read graph
VRP::VRP(bool isMap, const string& inputName)   :
    maxspeed{map},
    length{map},
    lat{map},
    lon{map},
    coords{lon, lat},
    c{g},
    t{g},
    q{g},
    startCols{g},
    nodeRows{g},
    arcUse{g},
    bestSolutionStartCols{g}
{
    if(isMap) {
        Timer timer(true);
        cout << "Reading the map... " << flush;
        digraphReader(map, string(inputName))
                .arcMap("maxspeed", maxspeed)
                .arcMap("length", length)
                .nodeMap("lat", lat)
                .nodeMap("lon", lon)
                .run();
        cout << "Elapsed: " << timer.realTime() << "s" << endl;
        mapNodesNumber = countNodes(map);
        mapArcsNumber = countArcs(map);

        Q=100;
        cout << "Number of nodes: " << mapNodesNumber << endl;
        cout << "Number of arcs: " << mapArcsNumber << endl << endl;
    } else {
        Timer timer(true);
        cout << "Reading the graph.. " << flush;
        digraphReader(g, inputName)
                .arcMap("c", c).arcMap("t", t)
                .nodeMap("q", q)
                .attribute("Q", Q)
                .run();
        n = countNodes(g);
        arcs.resize(n, vector<ListDigraph::Arc>(n));
        for (ListDigraph::ArcIt arc(g); arc != INVALID; ++arc) {
            arcs[g.id(g.source(arc))][g.id(g.target(arc))] = arc;
        }
        cout << "Elapsed: " << timer.realTime() << "s" << endl;
        cout << "Number of nodes: " << n << endl << endl;
    }
    //Check
    for(ListDigraph::NodeIt node(g); node!=INVALID; ++node){
        if(g.id(node)==0){
            q[node]=0;
        } else {
            myAssert(q[node]<=Q, "Too big demand(q) value!");
            myAssert(q[node]>0, "Too low demand(q) value'");
        }
    }
}

void VRP::init(int costumerCnt)
{
    g.clear();
    generateCostumersGraph(costumerCnt);
    printCostumerCoordinates();
    shortestPaths();
    //printShortestPathsFromDepot();
}

void VRP::generateCostumersGraph(int costumerCnt)
{
    Timer timer(true);
    cout << "Generate costumer graph.. " << flush;
    n=costumerCnt+1;
    depotAndCostumers.resize(0);
    depotAndCostumers.reserve(n);
    nodes.resize(0);
    nodes.reserve(n);
    depotAndCostumers.push_back(0);
    nodes.push_back(g.addNode());

    Random random(mySeed);
    for(int i = 1; i < n; ++i){
        depotAndCostumers.push_back(random[mapNodesNumber]);
        nodes.push_back(g.addNode());
    }
    arcs.resize(0);
    arcs.resize(n);
    for(int i = 0; i < n; ++i){
        arcs[i].resize(n);
        for(int j = 0; j < n; ++j){
            if(i != j) {
                arcs[i][j]=g.addArc(nodes[i], nodes[j]);
            }
        }
    }
    cout << "Generated demands: " << endl;
    for(ListDigraph::NodeIt node(g); node != INVALID; ++node){
        if(g.id(node)==0){
            q[node]=0;
        } else {
            q[node] = random[myMaxWeight]+1;
        }
        cout << q[node] << " ";
    }
    cout << endl;
    cout << "Elapsed: " << timer.realTime() << "s" << endl << endl;
}

void VRP::printToEps(const string& filename){
    Timer timer(true);
    cout << "Printing to eps..." << flush;
    ListDigraph::ArcMap<Color> arcColor(map, Color(170, 255, 0));
    ListDigraph::ArcMap<double> arcWidth(map, 0.035);
    for(unsigned int i=0; i<bestSolutionColIndexs.size(); ++i){
        int colIndex=bestSolutionColIndexs[i];
        for(unsigned int arcIndex=0; arcIndex<routes[colIndex].size(); ++arcIndex){
            ListDigraph::Arc arc=routes[colIndex][arcIndex];
            for(unsigned int a=0; a<paths[g.id(g.source(arc))][g.id(g.target(arc))].size(); ++a){
                ListDigraph::Arc currArc;
                currArc=paths[g.id(g.source(arc))][g.id(g.target(arc))][a];
                arcWidth[currArc]=0.15;
                arcColor[currArc]=Color(255, 0, 0);
            }
        }
    }
    for(ListDigraph::NodeIt node(g); node!=INVALID; ++node){
        if(g.id(node)!=0 && bestSolutionStartCols[node]>EPSILON){
            ListDigraph::Arc arc;
            arc=arcs[0][g.id(node)];
            for(unsigned int a=0; a<paths[g.id(g.source(arc))][g.id(g.target(arc))].size(); ++a){
                ListDigraph::Arc currArc;
                currArc=paths[g.id(g.source(arc))][g.id(g.target(arc))][a];
                arcWidth[currArc]=0.15;
                arcColor[currArc]=Color(255, 0, 0);
            }
            arc=arcs[g.id(node)][0];
            for(unsigned int a=0; a<paths[g.id(g.source(arc))][g.id(g.target(arc))].size(); ++a){
                ListDigraph::Arc currArc;
                currArc=paths[g.id(g.source(arc))][g.id(g.target(arc))][a];
                arcWidth[currArc]=0.15;
                arcColor[currArc]=Color(255, 0, 0);
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

void VRP::shortestPaths()
{
    Timer timer(true);
    cout << "Shortest paths..." << flush;
    ArcTravelTime travelTime(maxspeed, length);

    paths.resize(0);
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
    cout << " Elapsed: " << timer.realTime() << "s" << endl << endl;
}

void VRP::printCostumerCoordinates()
{
    cout << "Costumer coordinates: " << endl;
    for(int i = 0; i < n; ++i){
        ListDigraph::Node node;
        node=map.nodeFromId(depotAndCostumers[i]);
        cout << i << " : (" << lat[node] << ", " << lon[node] << ")";
        cout << "\t demand: " << q[nodes[i]] << endl;
    }
    cout << endl;
}

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
    return closestNode;
}

void VRP::printShortestPathsFromDepot()
{
    cout << "Shortest paths from depot: " << endl;
    for(int i = 1; i < n; ++i){
        ListDigraph::Node node;
        node=map.nodeFromId(depotAndCostumers[i]);
        cout << i << " : (" << lat[node] << ", " << lon[node] << ")" << endl;
        cout << "     time: " << t[arcs[0][i]] << " minutes" << endl;
        cout << "     distance: " << c[arcs[0][i]] << " meters" << endl;
    }
    cout << endl;
}

void VRP::createMasterLP()
{
    masterLP.clear();
    cols.resize(0);
    routes.resize(0);
    routeNodes.resize(0);

    //Add cols
    for(ListDigraph::NodeIt node(g); node !=INVALID; ++node){
        if(g.id(node)!=0){
            startCols[node]=masterLP.addCol();
            masterLP.colLowerBound(startCols[node], 0);
            masterLP.colUpperBound(startCols[node], 1);
            //Calculate startCol costs
            masterLP.objCoeff(
                    startCols[node],
                    c[arcs[0][g.id(node)]]+c[arcs[g.id(node)][0]]);
        }
    }
    vehicleNumberCol=masterLP.addCol();
    totalCostCol=masterLP.addCol();
    masterLP.objCoeff(vehicleNumberCol, 0);
    masterLP.objCoeff(totalCostCol, 0);
    masterLP.colLowerBound(vehicleNumberCol, 0);
    masterLP.colLowerBound(totalCostCol, 0);

    //Add rows
    for(ListDigraph::NodeIt node(g); node!=INVALID; ++node) {
        if (g.id(node) != 0) {
            nodeRows[node] = masterLP.addRow(startCols[node] >= 1);
        }
    }
    vehicleNumberRow=masterLP.addRow();
    totalCostRow=masterLP.addRow();

    //Set vehicleNumber and totalCost rows coeffs
    for(ListDigraph::NodeIt node(g); node !=INVALID; ++node){
        if(g.id(node) != 0){
            masterLP.coeff(vehicleNumberRow, startCols[node], 1);
            masterLP.coeff(totalCostRow, startCols[node],
                    masterLP.objCoeff(startCols[node]));
        }
    }
    masterLP.coeff(vehicleNumberRow, vehicleNumberCol, -1);
    masterLP.coeff(totalCostRow, totalCostCol, -1);
    masterLP.rowLowerBound(vehicleNumberRow, 0);
    masterLP.rowUpperBound(vehicleNumberRow, 0);
    masterLP.rowLowerBound(totalCostRow, 0);
    masterLP.rowUpperBound(totalCostRow, 0);

    masterLP.min();
}

void VRP::printMasterLPSolution()
{
    cout << "Master LP solution:" << endl;
    cout << "Number of added columns: " << cols.size() << endl;
    /*for(Lp::Col col : cols){
        if(masterLP.primal(col)>EPSILON) {
            cout << masterLP.primal(col) << endl;
        }
    }
    for(ListDigraph::NodeIt node(g); node !=INVALID; ++node) {
        if(g.id(node)!=0) {
            cout << g.id(node) << " . node: ";
            cout << masterLP.primal(startCols[node]) << endl;
        }
    }*/
    cout.precision(12);
    cout << "Vehicle number: " << masterLP.primal(vehicleNumberCol) << endl;
    cout << "Total cost: " << masterLP.primal(totalCostCol) << endl;


    //Calculate dual cost
    double dualCost=0;
    for(Lp::RowIt row(masterLP); row!=INVALID; ++row){
        dualCost+=masterLP.dual(row);
    }
    dualCost-=masterLP.dual(vehicleNumberRow);
    dualCost-=masterLP.dual(totalCostRow);
    //cout << "Dual cost: " << dualCost << endl;
    cout << endl;
}

class MarginalCost {
private:
    const VRP &vrp;
public:
    typedef const ListDigraph::Arc &Key;
    typedef double Value;

    MarginalCost(const VRP &in_vrp) :
            vrp{in_vrp}
    {}

    Value operator[](MarginalCost::Key arc) {
        ListDigraph::Node sNode, tNode;
        sNode = vrp.g.source(arc);
        tNode = vrp.g.target(arc);
        double pi_c = vrp.masterLP.dual(vrp.totalCostRow);
        double pi_s;
        if (vrp.g.id(sNode) == 0) {
            pi_s = vrp.masterLP.dual(vrp.vehicleNumberRow);
        } else {
            pi_s = vrp.masterLP.dual(vrp.nodeRows[sNode]);
        }
        double cost = 0;
        cost = (1 - pi_c) * vrp.c[
                vrp.arcs[vrp.g.id(sNode)][vrp.g.id(tNode)]
        ] - pi_s;
        return cost;
    }
};

bool VRP::solveMasterLP()
{
    Timer timer(true);
    cout << "Solving the master problem.. " << endl;
    int itCnt=0;
    do{
        ++itCnt;
        if(cols.size()%1000==0 && cols.size()>0){
            cout << "Added cols number: " << cols.size() << endl;
        }
        masterLP.solve();
        //myAssert(masterLP.primalType()==Lp::OPTIMAL, "Not optimal");
        switch(masterLP.primalType()) {
            case 0:
                cout << "undefined" << endl; return false;
            case 1:
                cout << "infeasible" << endl; return false;
            case 2:
                cout << "feasible" << endl; break;
            case 3:
                //cout << "optimal" << endl;
                break;
            case 4:
                cout << "unbounded" << endl; return false;
        }
        //printMasterLPMatrix();
        //printMasterLPSolution();
    } while(generateColumn() && itCnt<4000);

    cout << "Master LP solved, elapsed: " << timer.realTime() << "s" << endl;
    cout << "User time: " << timer.userTime() << "s" << endl;
    cout << "Number of added columns: " << cols.size() << endl;
    cout.precision(12);
    cout << "Vehicle number: " << masterLP.primal(vehicleNumberCol) << endl;
    cout << "Total cost: " << masterLP.primal(totalCostCol) << endl;
    cout << endl;

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


bool VRP::extendLabel( ListDigraph::NodeMap<NodeLabels>& nodeLabels,
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
            if (extendedL.nodeUse[otherNode] == 0
                && extendedL.weight + q[otherNode] <= Q)
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

bool VRP::generateColumn()
{
    //cout << "Generating column" << endl;
    static int generationCnt=0;
    ++generationCnt;
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
                if(generationCnt<5000){
                    while(nextCheckDepotIndex<static_cast<int>
                        (nodeLabels[depot].labelList.size()))
                    {
                        if(nodeLabels[depot].labelList[nextCheckDepotIndex].
                                cost<-EPSILON)
                        {
                            //cout << "Early exit!" << endl;
                            //cout << "Found min cost: "
                            //    << nodeLabels[depot]
                            //        .labelList[nextCheckDepotIndex]
                            //        .cost << endl;
                            //cout << "Elapsed: " << timer.realTime() << "s" << endl << endl;
                            addGeneratedColumn(nodeLabels[depot].
                                    labelList[nextCheckDepotIndex]);
                            return true;
                        } else {
                            ++nextCheckDepotIndex;
                        }
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
            //cout << minCost << endl;
            minIndex=i;
        } else {
            //cout << "Nope" << endl;
        }
    }
    //cout << "Found min cost: " << minCost << endl;
    //cout << "Elapsed: " << timer.realTime() << "s" << endl << endl;
    if(minCost<-EPSILON){
        addGeneratedColumn(nodeLabels[depot].labelList[minIndex]);
        return true;
    } else {
        return false;
    }
}

void VRP::addGeneratedColumn(const Label& l)
{
    int nodeCount=l.nodeCnt+1;

    Lp::Col col;
    col=masterLP.addCol();
    cols.push_back(col);
    routes.push_back(vector<ListDigraph::Arc> (nodeCount-1, INVALID));
    vector<ListDigraph::Node> currRouteNodes(nodeCount, INVALID);
    currRouteNodes[0]=g.nodeFromId(0);
    for(ListDigraph::NodeIt node(g); node != INVALID; ++node){
        if(l.nodeUse[node]!=0) {
            currRouteNodes[l.nodeUse[node]] = node;
        }
    }
    routeNodes.push_back(currRouteNodes);
    masterLP.colLowerBound(col, 0.0);
    int currCost=0;
    //cout << "Added column's nodes:  ";
    for(int i=0; i<nodeCount; ++i){
        //cout << g.id(currRouteNodes[i]) << " ";
        if(g.id(currRouteNodes[i]) != 0) {
            masterLP.coeff(nodeRows[currRouteNodes[i]], col, 1.0);
        }
        if(i>0){
            routes[routes.size()-1][i-1]=
                    arcs[g.id(currRouteNodes[i-1])][g.id(currRouteNodes[i])];
            currCost+=c[routes[routes.size()-1][i-1]];
        }
    }
    //cout << endl;
    masterLP.coeff(vehicleNumberRow, col, 1);
    masterLP.coeff(totalCostRow, col, currCost);
    masterLP.objCoeff(col, currCost);
}

void VRP::printRoutes(int index){
    if(0<=index && index<static_cast<int>(routeNodes.size())){
        cout << "Print route with index: " << index << endl;
        for(unsigned int j=0; j<routeNodes[index].size(); ++j){
            cout << g.id(routeNodes[index][j]) << " ";
        }
        cout << endl;
        cout << "Cost: " << masterLP.objCoeff(cols[index]) << endl;
    } else if(index==-1) {
        cout << "Print routes: " << endl;
        for (unsigned int i = 0; i < routeNodes.size(); ++i) {
            for(unsigned int j=0; j<routeNodes[i].size(); ++j){
                cout << g.id(routeNodes[i][j]) << " ";
            }
            cout << "    Cost: " << masterLP.objCoeff(cols[i]) << endl;
        }
        cout << endl;
    }
}

void VRP::printMasterLPMatrix(){
    cout << "PRINT MASTERLP MATRIX: " << endl;
    cout << "\t";
    for(Lp::ColIt col(masterLP); col!=INVALID; ++col){
        cout << masterLP.primal(col) << "\t  ";
    }
    cout << "|  " << endl;


    cout << "\t";
    for(Lp::ColIt col(masterLP); col!=INVALID; ++col){
        cout << "-" << "\t  ";
    }
    cout << "|  " << endl;


    for(Lp::RowIt row(masterLP); row !=INVALID; ++row){
        cout << masterLP.dual(row) << " |\t";
        for(Lp::ColIt col(masterLP); col!=INVALID; ++col){
            cout << masterLP.coeff(row, col) << "\t  ";
        }
        cout << "|  ";
        cout << masterLP.rowLowerBound(row) << "  ";
        cout << masterLP.rowUpperBound(row) << endl;
    }


    cout << "\t";
    for(Lp::ColIt col(masterLP); col!=INVALID; ++col){
        cout << "-" << "\t  ";
    }
    cout << "|  " << endl;


    cout << "\t";
    for(Lp::ColIt col(masterLP); col!=INVALID; ++col){
        cout << masterLP.objCoeff(col) << "\t  ";
    }
    cout << "|  " << endl;

    cout << endl;
}


void VRP::checkMIP(bool printEps)
{
    Timer timer(true);
    cout << "Check MIP started: " << endl;

    Mip mip;
    //Add cols
    ListDigraph::ArcMap<Mip::Col> checkCols(g);
    for(ListDigraph::ArcIt arc(g); arc!=INVALID; ++arc){
        checkCols[arc]=mip.addCol();
        mip.colLowerBound(checkCols[arc], 0);
        mip.colType(checkCols[arc], Mip::INTEGER);
    }
    ListDigraph::NodeMap<Lp::Col> checkWeights(g);
    for(ListDigraph::NodeIt node(g); node !=INVALID; ++node){
        checkWeights[node]=mip.addCol();
        mip.colLowerBound(checkWeights[node], 0);
        mip.colUpperBound(checkWeights[node], Q);
    }

    //Add rows
    for(ListDigraph::NodeIt node(g); node !=INVALID; ++node){
        if(g.id(node)!=0) {
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
    for(ListDigraph::NodeIt node(g); node !=INVALID; ++node){
        if(g.id(node)!=0) {
            for (ListDigraph::InArcIt arc(g, node); arc != INVALID; ++arc) {
                Lp::Expr e;
                e = (checkWeights[g.source(arc)] + q[node]);
                mip.addRow(checkWeights[node]<=BIG_VALUE*(1-checkCols[arc])+e);
                mip.addRow(checkWeights[node]+BIG_VALUE*(1-checkCols[arc])>=e);
            }
       }
    }
    for(ListDigraph::ArcIt arc(g); arc!=INVALID; ++arc){
        mip.objCoeff(checkCols[arc], c[arc]);
    }
    mip.min();
    mip.solve();
    switch(mip.type()) {
        case 0:
            cout << "undefined" << endl; break;
        case 1:
            cout << "infeasible" << endl; break;
        case 2:
            cout << "feasible" << endl; break;
        case 3:
            cout << "optimal" << endl; break;
        case 4:
            cout << "unbounded" << endl; break;
    }
    cout << "CHECK MIP primal: " << mip.solValue() << endl;
    cout << "Elapsed: " << timer.realTime() << "s" << endl;
    cout << "User time: " << timer.userTime() << "s" << endl;

    if(printEps){
        printToEpsCheckMIP("checkMIP.eps", mip, checkCols);
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

void VRP::branchAndBound()
{
    Timer timer(true);
    bestCost=BIG_VALUE;
    cout << "Branch and bound started: " << endl;
    int branchedNodes=0;
    createMasterLP();
    recursiveBranch(branchedNodes);

    cout << endl;
    cout << "Finished" << endl;
    cout << "Arc count of g: " << countArcs(g) << endl;
    cout << "Elapsed real time: " << timer.realTime() << "s" << endl;
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

void VRP::recursiveBranch(int& branchedNodes)
{
    ++branchedNodes;
    cout << "New branching node: " << branchedNodes << endl;
    if(!solveMasterLP()){
        return;
    } else if(masterLP.primal()>bestCost) {
        cout << "BOUND" << endl;
        return;
    }

    double vehicleNumber=masterLP.primal(vehicleNumberCol);
    Lp::Row tempRow=INVALID;
    if(!isWhole(vehicleNumber)){
        //branch on vehicle number
        tempRow=masterLP.addRow(vehicleNumberCol<=floor(vehicleNumber));
        cout << "Add: vehicleNumberCol <= " << floor(vehicleNumber) << endl;
        recursiveBranch(branchedNodes);
        cout << "Remove vehicleNumberCol <= " << floor(vehicleNumber) << endl;
        masterLP.erase(tempRow);

        tempRow=masterLP.addRow(vehicleNumberCol>=ceil(vehicleNumber));
        cout << "Add: vehicleNumberCol >= " << ceil(vehicleNumber) << endl;
        recursiveBranch(branchedNodes);
        cout << "Remove: vehicleNumberCol >= " << ceil(vehicleNumber) << endl;
        masterLP.erase(tempRow);
    } else {
        //branch on arcs
        calculateArcUse();
        double maxFractionalUse=0.0;
        ListDigraph::Arc arcToBranch=INVALID;
        for(ListDigraph::ArcIt arc(g); arc!=INVALID; ++arc){
            if(!isWhole(arcUse[arc])){
                cout << arcUse[arc] << "  ";
                if(arcUse[arc]>maxFractionalUse){
                    maxFractionalUse=arcUse[arc];
                    arcToBranch=arc;
                }
            }
        }
        //myAssert(arcUse[arcToBranch]<1, "Arc with >1 use");
        cout  << endl;
        if(arcToBranch==INVALID){
            cout.precision(12);
            cout << "WHOLE solution found with cost: ";
            cout << masterLP.primal();

            if(masterLP.primal()<bestCost){
                cout << " Better" << endl;
                bestCost=masterLP.primal();
                bestSolutionVehicle=masterLP.primal(vehicleNumberCol);
                bestSolutionColIndexs.resize(0);
                for(unsigned int i=0; i<cols.size(); ++i){
                    if(masterLP.primal(cols[i])>EPSILON) {
                        cout << i << " .col,";
                        cout << " value: ";
                        cout << masterLP.primal(cols[i]) << endl;
                        cout << "Route: ";
                        printRoutes(i);
                        cout << "Route cost: " << masterLP.objCoeff(cols[i]) << endl;
                        cout << endl;
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
                    bestSolutionStartCols[node]=masterLP.primal(startCols[node]);
                }
            } else {
                cout << " Worse" << endl;
            }
        } else {
            //Branch on this arc
            cout << "Branch on arc with use: ";
            cout << arcUse[arcToBranch] << endl;

            //0
            cout << "Add: arc =0" << endl;
            int originalArcCount=countArcs(g); //save original arc count
            ListDigraph::Node sourceNode=g.source(arcToBranch);
            ListDigraph::Node targetNode=g.target(arcToBranch);
            int originalArcCost=c[arcToBranch];
            vector<pair<int, int>> changedCosts;    // <index, oldCost>
            changeObjCoeffs(arcToBranch, changedCosts);
            g.erase(arcToBranch);
            int startOldCost=-1;
            if(g.id(sourceNode)==0){
                startOldCost=masterLP.objCoeff(startCols[targetNode]);
                masterLP.objCoeff(startCols[targetNode], BIG_VALUE);
            } else if(g.id(targetNode)==0){
                startOldCost=masterLP.objCoeff(startCols[sourceNode]);
                masterLP.objCoeff(startCols[sourceNode], BIG_VALUE);
            }
            recursiveBranch(branchedNodes);
            for(unsigned int i=0; i<changedCosts.size(); ++i){
                masterLP.objCoeff(cols[changedCosts[i].first], changedCosts[i].second);
            }
            if(g.id(sourceNode)==0){
                masterLP.objCoeff(startCols[targetNode], startOldCost);
            } else if(g.id(targetNode)==0){
                masterLP.objCoeff(startCols[sourceNode], startOldCost);
            }
            arcs[g.id(sourceNode)][g.id(targetNode)]=g.addArc(sourceNode, targetNode);
            c[arcs[g.id(sourceNode)][g.id(targetNode)]]=originalArcCost;
            myAssert(countArcs(g)==originalArcCount, "Arc(s) vanished!");
            cout << "Remove: arc =0" << endl;


            //arc=1
            cout << "Add: arc =1" << endl;
            changedCosts.resize(0);
            vector<tuple<ListDigraph::Node, ListDigraph::Node, int>> erasedArcs;
            for(ListDigraph::OutArcIt arc(g, sourceNode); arc!=INVALID;){
                ListDigraph::Arc a=arc;
                ++arc;
                changeObjCoeffs(a, changedCosts);
                erasedArcs.push_back(tuple<ListDigraph::Node, ListDigraph::Node, int>
                        (g.source(a), g.target(a), c[a]));
                g.erase(a);
            }
            for(ListDigraph::InArcIt arc(g, targetNode); arc!=INVALID;){
                ListDigraph::Arc a=arc;
                ++arc;
                changeObjCoeffs(a, changedCosts);
                erasedArcs.push_back(tuple<ListDigraph::Node, ListDigraph::Node, int>
                        (g.source(a), g.target(a), c[a]));
                g.erase(a);
            }
            int startSourceOldCost=-1;
            int startTargetOldCost=-1;
            if(g.id(sourceNode)!=0) {
                startSourceOldCost = masterLP.objCoeff(startCols[sourceNode]);
                masterLP.objCoeff(startCols[sourceNode], BIG_VALUE);
            }
            if(g.id(targetNode)!=0){
                startTargetOldCost = masterLP.objCoeff(startCols[targetNode]);
                masterLP.objCoeff(startCols[targetNode], BIG_VALUE);
            }
            recursiveBranch(branchedNodes);
            for(unsigned int i=0; i<changedCosts.size(); ++i){
                masterLP.objCoeff(cols[changedCosts[i].first], changedCosts[i].second);
            }
            if(g.id(sourceNode)!=0){
                masterLP.objCoeff(startCols[targetNode], startSourceOldCost);
            } else if(g.id(targetNode)!=0){
                masterLP.objCoeff(startCols[sourceNode], startTargetOldCost);
            }
            for(unsigned int i=0; i<erasedArcs.size(); ++i){
                ListDigraph::Node sNode=std::get<0>(erasedArcs[i]);
                ListDigraph::Node tNode=std::get<1>(erasedArcs[i]);
                arcs[g.id(sNode)][g.id(tNode)]
                    =g.addArc(sNode, tNode);
                c[arcs[g.id(sNode)][g.id(tNode)]]=std::get<2>(erasedArcs[i]);
            }
            myAssert(countArcs(g)==originalArcCount, "Arc(s) vanished!");
            cout << "Remove: arc =1" << endl;
        }
    }
}

void VRP::changeObjCoeffs(ListDigraph::Arc arc, vector<pair<int, int>>& changedCosts)
{
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
}

void VRP::calculateArcUse(){
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

void VRP::printCost(int sourceId, int targetId){
    //TODO: check is valid id?
    cout << c[arcs[sourceId][targetId]] << endl;
}