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


VRP::VRP(string inputMapName)   :
    maxspeed{map},
    length{map},
    lat{map},
    lon{map},
    coords{lon, lat}
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
    n=in_n;
    depotAndCostumers.reserve(n+1);
    nodes.reserve(n+1);
    depotAndCostumers.push_back(0);
    nodes.push_back(g.addNode());

    Random random(42);
    for(int i = 1; i <= n; ++i){
        depotAndCostumers.push_back(random[mapNodesNumber]);
        nodes.push_back(g.addNode());
    }
    arcs.resize(n+1);
    for(int i = 0; i <=n; ++i){
        arcs[i].reserve(n+1);
        for(int j = 0; j <= n; ++j){
            arcs[i].push_back(g.addArc(nodes[i], nodes[j]));
        }
    }
}
