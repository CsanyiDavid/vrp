//
// Created by david on 2021. 01. 13..
//

#include <lemon/lp.h>
#include <vector>

using namespace std;
using namespace lemon;

int main(){
    Lp lp;
    Lp::Col x = lp.addCol();
    Lp::Col y = lp.addCol();
    Lp::Row r;
    r= lp.addRow(4*x+2*y <= 3);
    Lp::Row r2;
    r2=lp.addRow(y>=0);
    lp.obj(x);
    lp.max();
    lp.solve();
    cout << lp.primal() << endl;

    cout << lp.coeff(r, x) << " " << lp.coeff(r, y) << endl;
    cout << lp.coeff(r2, x) << " " << lp.coeff(r2, y) << endl;
    return 0;
}