//
// Created by david on 2020. 12. 22..
//

#include "vrp.h"

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
    nodeRows{g}
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
        n=countNodes(g);
        arcs.resize(n, vector<ListDigraph::Arc> (n));
        for(ListDigraph::ArcIt arc(g); arc!=INVALID; ++arc){
            arcs[g.id(g.source(arc))][g.id(g.target(arc))]=arc;
        }

        //Check
        for(ListDigraph::NodeIt node(g); node!=INVALID; ++node){
            if(g.id(node)==0){
                q[node]=0;
            } else {
                myAssert(q[node]<=Q, "Too big demand(q) value!");
            }
        }
        cout << "Elapsed: " << timer.realTime() << "s" << endl;
        cout << "Number of nodes: " << n << endl << endl;
    }
}

void VRP::generateCostumersGraph(int costumerCnt)
{
    Timer timer(true);
    cout << "Generate costumer graph.. " << flush;
    n=costumerCnt+1;
    depotAndCostumers.reserve(n);
    nodes.reserve(n);
    depotAndCostumers.push_back(0);
    nodes.push_back(g.addNode());

    Random random(42);
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
    for(ListDigraph::NodeIt node(g); node != INVALID; ++node){
        q[node]=random[20]+1;
    }
    cout << "Elapsed: " << timer.realTime() << "s" << endl << endl;
}

void VRP::printToEps(const string& filename){
    Timer timer(true);
    cout << "Printing to eps..." << flush;
    ListDigraph::ArcMap<Color> arcColor(map, Color(0, 0, 0));
    ListDigraph::ArcMap<double> arcWidth(map, 0.0);
    for(int j=0; j<n; ++j) {
        for(unsigned int k=0; k<paths[0][j].size(); ++k){
            arcWidth[paths[0][j][k]]=0.1;
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
        cout << i << " : (" << lat[node] << ", " << lon[node] << ")" << endl;
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
    for(Lp::Col col : cols){
        if(masterLP.primal(col)>EPSILON) {
            cout << masterLP.primal(col) << endl;
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
    cout << "Total cost: " << masterLP.primal(totalCostCol) << endl;


    //Calculate dual cost
    double dualCost=0;
    for(Lp::RowIt row(masterLP); row!=INVALID; ++row){
        dualCost+=masterLP.dual(row);
    }
    dualCost-=masterLP.dual(vehicleNumberRow);
    dualCost-=masterLP.dual(totalCostRow);
    cout << "Dual cost: " << dualCost << endl;
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

void VRP::solveMasterLP()
{
    Timer timer(true);
    cout << "Solving the master problem.. " << endl;
    int itCnt=0;
    do{
        ++itCnt;
        masterLP.solve();
        myAssert(masterLP.primalType()==Lp::OPTIMAL, "Not optimal");
        //printMasterLPMatrix();
        printMasterLPSolution();
    } while(generateColumn() && itCnt<1000);

    cout << "Master LP solved, elapsed: " << timer.realTime() << "s" << endl;
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
    cout << "Generating column" << endl;
    Timer timer(true);
    MarginalCost mc{*this};
    ListDigraph::NodeMap<NodeLabels> nodeLabels(g);
    ListDigraph::Node depot=g.nodeFromId(0);
    nodeLabels[depot].labelList.emplace_back(g);

    int notExtendedCnt=0;
    extendLabel(nodeLabels, depot, mc);

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
    cout << "Found min cost: " << minCost << endl;
    cout << "Elapsed: " << timer.realTime() << "s" << endl << endl;
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
    masterLP.colLowerBound(col, 0.0);
    masterLP.colUpperBound(col, 1.0);
    int currCost=0;
    cout << "Added column's nodes:  ";
    for(int i=0; i<nodeCount; ++i){
        cout << g.id(currRouteNodes[i]) << " ";
        if(g.id(currRouteNodes[i]) != 0) {
            masterLP.coeff(nodeRows[currRouteNodes[i]], col, 1.0);
        }
        if(i>0){
            routes[routes.size()-1][i-1]=
                    arcs[g.id(currRouteNodes[i-1])][g.id(currRouteNodes[i])];
            currCost+=c[routes[routes.size()-1][i-1]];
        }
    }
    cout << endl;
    masterLP.coeff(vehicleNumberRow, col, 1);
    masterLP.coeff(totalCostRow, col, currCost);
    masterLP.objCoeff(col, currCost);
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


void VRP::checkLP() {
    Timer timer(true);
    cout << "Check LP started: " << endl;

    //Add cols
    ListDigraph::ArcMap<Lp::Col> cols(g);
    for(ListDigraph::ArcIt arc(g); arc!=INVALID; ++arc){
        cols[arc]=lp.addCol();
        lp.colLowerBound(cols[arc], 0);
        lp.colUpperBound(cols[arc], 1);
        lp.colType(cols[arc], Mip::INTEGER);
    }
    ListDigraph::NodeMap<Lp::Col> weights(g);
    for(ListDigraph::NodeIt node(g); node !=INVALID; ++node){
        weights[node]=lp.addCol();
        lp.colLowerBound(weights[node], 0);
        lp.colUpperBound(weights[node], Q);
    }

    //Add rows
    for(ListDigraph::NodeIt node(g); node !=INVALID; ++node){
        if(g.id(node)!=0) {
            Lp::Expr eOut, eIn;
            for (ListDigraph::OutArcIt arc(g, node); arc != INVALID; ++arc) {
                eOut += cols[arc];
            }
            lp.addRow(eOut == 1);
            for (ListDigraph::InArcIt arc(g, node); arc != INVALID; ++arc) {
                eIn += cols[arc];
            }
            lp.addRow(eIn == 1);
        }
    }
    lp.colUpperBound(weights[g.nodeFromId(0)], 0);
    for(ListDigraph::NodeIt node(g); node !=INVALID; ++node){
        if(g.id(node)!=0) {
            for (ListDigraph::InArcIt arc(g, node); arc != INVALID; ++arc) {
                Lp::Expr e;
                e = (weights[g.source(arc)] + q[node]);
                lp.addRow(weights[node]<=BIG_VALUE*(1-cols[arc])+e);
                lp.addRow(weights[node]+BIG_VALUE*(1-cols[arc])>=e);
            }
       }
    }
    for(ListDigraph::ArcIt arc(g); arc!=INVALID; ++arc){
        lp.objCoeff(cols[arc], c[arc]);
    }
    lp.min();
    lp.solve();
    cout << "CHECK LP primal: " << lp.solValue() << endl;
    cout << "Elapsed: " << timer.realTime() << "s" << endl;
}
/*
void VRP::printToEpsCheckLp(const string& filename){
    Timer timer(true);
    cout << "Printing to eps..." << flush;
    ListDigraph::ArcMap<Color> arcColor(map, Color(0, 0, 0));
    ListDigraph::ArcMap<double> arcWidth(map, 0.0);
    for(ListDigraph::ArcIt arc(g); arc!=INVALID; ++arc){
        if()
    }


    for(int j=0; j<n; ++j) {
        for(unsigned int k=0; k<paths[0][j].size(); ++k){
            arcWidth[paths[0][j][k]]=0.1;
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
}*/