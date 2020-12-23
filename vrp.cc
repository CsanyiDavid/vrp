//
// Created by david on 2020. 12. 22..
//

#include "vrp.h"

using namespace std;
using namespace lemon;

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


VRP::VRP(string inputMapName)   :
    maxspeed{map},
    length{map},
    lat{map},
    lon{map},
    coords{lon, lat},
    ids{g},
    c{g},
    t{g}
{
    Timer timer(true);
    cout << "Reading the map... " << flush;
    digraphReader(map, string(inputMapName))
            .arcMap("maxspeed",maxspeed)
            .arcMap("length",length)
            .nodeMap("lat",lat)
            .nodeMap("lon",lon)
            .run();
    cout << "Elapsed: " << timer.realTime() << "s" << endl;
    mapNodesNumber=countNodes(map);
    mapArcsNumber=countArcs(map);
    cout << "Number of nodes: " << mapNodesNumber << endl;
    cout << "Number of arcs: " << mapArcsNumber << endl;
}

void VRP::generateCostumersGraph(int in_n)
{
    Timer timer(true);
    cout << "Generate costumer graph..." << flush;
    n=in_n+1;
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
    cout << " Elapsed: " << timer.realTime() << "s" << endl;
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
    cout << " Elapsed: " << timer.realTime() << "s" << endl;
}

class ArcTravelTime{
private:
    const ListDigraph::ArcMap<int>& maxspeed;
    const ListDigraph::ArcMap<int>& length;
public:
    typedef int Value;
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
            }
        }
    }
    cout << " Elapsed: " << timer.realTime() << "s" << endl;
}

void VRP::printCostumerCoordinates()
{
    cout << "Costumer coordinates: " << endl;
    for(int i=0; i<n ; ++i){
        ListDigraph::Node node;
        node=map.nodeFromId(depotAndCostumers[i]);
        cout << i << " : (" << lat[node] << ", " << lon[node] << ")" << endl;
    }
    cout << endl;
}

ListDigraph::Node VRP::nodeIdFromLatLon(double latitude, double longitude)
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