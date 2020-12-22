//
// Created by david on 2020. 12. 22..
//

#include <iostream>
#include <lemon/dim2.h>
#include <lemon/list_graph.h>
#include <lemon/lgf_reader.h>
#include <lemon/time_measure.h>
#include <string>

#ifndef VRP_H
#define VRP_H

using namespace std;
using namespace lemon;

class CoordMap{
    ListDigraph::NodeMap<double>& lon;
    ListDigraph::NodeMap<double>& lat;
public:
    typedef ListDigraph::Node Key;
    typedef dim2::Point<double> Value;

    CoordMap(ListDigraph::NodeMap<double>& _lon,
             ListDigraph::NodeMap<double>& _lat) : lon(_lon), lat(_lat){
    }

    Value operator[](const Key& node) const{
        const double Phi = (lat[node] /* + 0.9430/60 */)/180*lemon::PI;
        const double Lam = (lon[node] /* + 4.0495/60 */)/180*lemon::PI;

        const double r = 6379743.001;
        const double k = 1.003110007693;
        const double n = 1.000719704936;
        const double eps = 0.0818205679407;
        const double Lam0 = (19+2.0/60+54.8584/3600)/180*lemon::PI;
        const double m0 = 0.9996;
        const double phi0 = 47.1/180*lemon::PI;

        const double sinPhi = std::sin(Phi);

        const double phi
                = 2*atan(k*std::pow(tan(lemon::PI_4+Phi/2),n)
                         *std::pow((1-eps*sinPhi)/(1+eps*sinPhi),n*eps/2))-lemon::PI_2;

        const double lam = n*(Lam-Lam0);

        const double phi2 = std::asin(std::cos(phi0)*std::sin(phi)
                                      -std::sin(phi0)*std::cos(phi)*std::cos(lam));
        const double lam2 = std::asin(std::cos(phi)*std::sin(lam)/std::cos(phi2));

        const double y = r*m0*lam2 + 650000;
        const double x = r*m0*std::log(std::tan(lemon::PI_4+phi2/2)) + 200000;

        return lemon::dim2::Point<double>(y,x);
    }
};

class VRP{
private:
    //The map
    ListDigraph map;
    ListDigraph::ArcMap<int> maxspeed;
    ListDigraph::ArcMap<int> length;
    ListDigraph::NodeMap<double> lat;
    ListDigraph::NodeMap<double> lon;

    //The graph
    ListDigraph g;

public:
    VRP(string inputMapName);

};


#endif /* VRP_H */
