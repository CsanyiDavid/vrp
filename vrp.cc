//
// Created by david on 2020. 12. 22..
//

#include "vrp.h"

using namespace std;
using namespace lemon;

void myAssert(bool bo, string errorType)
{
    if(!bo) {
        cout << "Assertion error: " << errorType << endl;
        exit(-1);
    }
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
VRP::VRP(bool isMap, string inputName)   :
    maxspeed{map},
    length{map},
    lat{map},
    lon{map},
    coords{lon, lat},
    ids{g},
    c{g},
    t{g},
    a{g},
    b{g},
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
        cout << "Number of nodes: " << mapNodesNumber << endl;
        cout << "Number of arcs: " << mapArcsNumber << endl << endl;
    } else {
        digraphReader(g, inputName)
                .arcMap("c", c).arcMap("t", t)
                .nodeMap("a",a).nodeMap("b",b)
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
                myAssert(a[node]<=b[node], "Incorrect time window");
            }
        }
    }
}

void VRP::generateCostumersGraph(int costumerCnt)
{
    Timer timer(true);
    cout << "Generate costumer graph..." << flush;
    n=costumerCnt+1;
    depotAndCostumers.reserve(n);
    nodes.reserve(n);
    depotAndCostumers.push_back(0);
    nodes.push_back(g.addNode());
    ids[nodes[0]]=0;

    Random random(42);
    for(int i = 1; i < n; ++i){
        depotAndCostumers.push_back(random[mapNodesNumber]);
        nodes.push_back(g.addNode());
        ids[nodes[i]]=i;
    }
    arcs.resize(n);
    for(int i = 0; i < n; ++i){
        arcs[i].reserve(n);
        for(int j = 0; j < n; ++j){
            arcs[i].push_back(g.addArc(nodes[i], nodes[j]));
        }
    }
    cout << " Elapsed: " << timer.realTime() << "s" << endl << endl;
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

    //ListDigraph::NodeMap<string> nodeText(map, "");
    //ListDigraph::NodeMap<Color> nodeTextColor(map, Color(255, 255, 0));
    ListDigraph::Node node=INVALID;
    for(int i = 1; i < n; ++i) {
        node=map.nodeFromId(depotAndCostumers[i]);
        nodeSize[node] = 0.2;
        //nodeText[node] = to_string(i);
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
            //.nodeTexts(nodeText)
            //.nodeTextColors(nodeTextColor)
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
        //Dijkstra<ListDigraph> dijkstra(map, length);
        dijkstra.run(startNode);

        ListDigraph::Node currNode=INVALID;
        ListDigraph::Arc currArc=INVALID;

        for(int j = 0; j < n; ++j){
            if(j != i){
                currNode=map.nodeFromId(depotAndCostumers[j]);

                c[arcs[i][j]]=0;
                t[arcs[i][j]]=0;
                while(currNode != startNode){
                    currArc=dijkstra.predArc(currNode);
                    c[arcs[i][j]]+=length[currArc];
                    t[arcs[i][j]]+=travelTime[currArc];
                    paths[i][j].push_back(currArc);
                    currNode = dijkstra.predNode(currNode);
                }
            } else {
                //c[arcs[i][j]]=0; ????????
                //t[arcs[i][j]]=0; ????????
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

void VRP::updateTimeWindows()
{
    int T=0;
    int curr_t;
    for(ListDigraph::NodeIt node(g); node !=INVALID; ++node){
        curr_t=b[node]+t[arcs[g.id(node)][0]];
        if(curr_t>T){
            T=curr_t;
        }
    }
    a[g.nodeFromId(0)]=0;
    b[g.nodeFromId(0)]=T;
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
            myAssert(c[arcs[0][g.id(node)]]<=b[node],
                     "Costumer can't be served in time!");
            masterLP.objCoeff(
                    startCols[node],
                    c[arcs[0][g.id(node)]]+c[arcs[g.id(node)][0]]);
        }
    }
    vehicleNumberCol=masterLP.addCol();
    totalCostCol=masterLP.addCol();
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

    //Set vehicleNUmber and totalCost rows coeffs
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
    for(Lp::Col col : cols){
        cout << masterLP.primal(col) << endl;
    }
    for(ListDigraph::NodeIt node(g); node !=INVALID; ++node) {
        if(g.id(node)!=0) {
            cout << g.id(node) << " . node: ";
            cout << masterLP.primal(startCols[node]) << endl;
        }
    }
    cout << "Vehicle number: " << masterLP.primal(vehicleNumberCol) << endl;
    cout << "Total cost: " << masterLP.primal(totalCostCol) << endl;
}

void VRP::solveMasterLP()
{
    do{
        masterLP.solve();
        printMasterLPSolution();
    } while(generateColumn());
}

class MarginalCost{
private:
    const VRP& vrp;
public:
    typedef const ListDigraph::Arc& Key;
    typedef double Value;

    MarginalCost(const VRP& in_vrp) :
        vrp{in_vrp}
    {}

    Value operator[](Key arc){
        ListDigraph::Node sNode, tNode;
        sNode=vrp.g.source(arc);
        tNode=vrp.g.target(arc);
        double pi_c=vrp.masterLP.dual(vrp.totalCostRow);
        double pi_s;
        if(vrp.g.id(sNode)==0){
            pi_s=vrp.masterLP.dual(vrp.vehicleNumberRow);
        } else {
            pi_s=vrp.masterLP.dual(vrp.nodeRows[sNode]);
        }
        int cost=0;
        cost=(1-pi_c)*vrp.c[
                vrp.arcs[vrp.g.id(sNode)][vrp.g.id(tNode)]
                ]-pi_s;
        return cost;
    }
};

class Labels{
public:
    int lower;
    int upper;
    int value;
    bool equal; //True is lower==upper

    Labels()    :
        lower(-BIG_VALUE),
        upper(BIG_VALUE),
        equal(false)
    {
    }
};

bool VRP::generateColumn()
{
    //Step 1 : Initialization
    MarginalCost mc(*this);
    ListDigraph::Node node=INVALID;
    vector<vector<vector<Labels>>> G(n, vector<vector<Labels>> (Q+1));
    int T=b[g.nodeFromId(0)];
    for(int jj=0; jj<n; ++jj) {
        for(int qq=0; qq<=Q; ++qq){
            node=g.nodeFromId(jj);
            G[jj][qq].resize(b[node]);
        }
    }
    G[0][0][0].upper=0;
    G[0][0][0].lower=0;
    G[0][0][0].value=0;
    G[0][0][0].equal=true;

    //Step 2 : Search for (q,t) to be treated
    ListDigraph::Node jNode, iNode;

    for(int qq=0; qq<=Q; ++qq){
        for(int tt=0; tt<=T; ++tt){
            for(int jj=0; jj<n; ++jj){
                jNode=g.nodeFromId(jj);


                if(a[jNode]<=tt && tt<=b[jNode]){
                    for(int ii=0; ii<n; ++ii){
                        iNode=g.nodeFromId(ii);
                        for(int qq2=0; qq2<=Q; ++qq2){
                            for(int tt2=a[iNode]; tt2<=b[iNode]; ++tt2){

                                if(tt2+t[arcs[ii][jj]]<=tt && qq2+q[jNode]<=qq) {
                                    int currUpper
                                            = G[ii][qq2][tt2].upper + mc[arcs[ii][jj]];
                                    int currLower
                                            = G[ii][qq2][tt2].upper + mc[arcs[ii][jj]];
                                    if (currUpper < G[jj][qq][tt].upper) {
                                        G[jj][qq][tt].upper = currUpper;
                                    }
                                    if (currLower < G[jj][qq][tt].lower) {
                                        G[jj][qq][tt].lower = currLower;
                                    }
                                }
                            }
                        }

                    }

                }

                //myAssert(G[jj][qq][tt].upper==G[jj][qq][tt].lower, "upper != lower");
            }
        }
    }
    cout << G[0][Q][T].lower << "   " << G[0][Q][T].upper << endl;
    return false;
}
