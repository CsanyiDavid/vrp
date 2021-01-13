//
// Created by david on 2021. 01. 13..
//

#include <lemon/lp.h>
#include <vector>

using namespace std;
using namespace lemon;

void myAssert(bool bo, const string& errorType)
{
    if(!bo) {
        cout << "Assertion error: " << errorType << endl;
        exit(-1);
    }
}

void printLPMatrix(const Lp& lp);

int main(){
    Lp lp;
    vector<Lp::Col> cols;
    for(int i=0; i<9; ++i){
        cols.push_back(lp.addCol());
        lp.colLowerBound(cols[i], 0);
        lp.colUpperBound(cols[i], 1);
    }

    vector<Lp::Row> rows;
    for(int i=0; i<9; ++i){
        rows.push_back(lp.addRow(cols[i]>=1));
    }
    Lp::Expr o;
    o=14*cols[0]
      +6*cols[1]
      +8*cols[2]
      +12*cols[3]
      +14*cols[4]
      +26*cols[5]
      +24*cols[6]
      +16*cols[7]
      +12*cols[8];
    lp.obj(o);
    lp.min();
    cout << lp.solve() << endl;
    printLPMatrix(lp);

    //add new col
    cols.push_back(lp.addCol());
    lp.colLowerBound(cols[9], 0);
    lp.colUpperBound(cols[9], 1);

    lp.coeff(rows[4], cols[9], 1);
    lp.coeff(rows[5], cols[9], 1);
    lp.coeff(rows[6], cols[9], 1);
    lp.coeff(rows[7], cols[9], 1);
    lp.coeff(rows[8], cols[9], 1);
    lp.objCoeff(cols[9], 36);

    cout << lp.solve() << endl;
    lp.solve();
    printLPMatrix(lp);

    //add new col
    cols.push_back(lp.addCol());
    lp.colLowerBound(cols[10], 0);
    lp.colUpperBound(cols[10], 1);

    lp.coeff(rows[4], cols[10], 1);
    lp.coeff(rows[5], cols[10], 1);
    lp.coeff(rows[6], cols[10], 1);
    lp.coeff(rows[7], cols[10], 1);
    lp.coeff(rows[8], cols[10], 1);
    lp.objCoeff(cols[10], 36);

    cout << lp.solve() << endl;
    printLPMatrix(lp);
}

void printLPMatrix(const Lp& lp){
    cout << "PRINT LP MATRIX: " << endl;
    cout << "Primal: " << lp.primal() << endl;
    if(lp.primalType()==Lp::OPTIMAL){
        cout << "Optimal" << endl;
    } else {
        cout << "Not optimal" << endl;
    }

    cout << "\t";
    for(Lp::ColIt col(lp); col!=INVALID; ++col){
        cout << lp.primal(col) << "\t  ";
    }
    cout << "|  " << endl;


    cout << "\t";
    for(Lp::ColIt col(lp); col!=INVALID; ++col){
        cout << "-" << "\t  ";
    }
    cout << "|  " << endl;


    for(Lp::RowIt row(lp); row !=INVALID; ++row){
        cout << lp.dual(row) << " |\t";
        for(Lp::ColIt col(lp); col!=INVALID; ++col){
            cout << lp.coeff(row, col) << "\t  ";
        }
        cout << "|  ";
        cout << lp.rowLowerBound(row) << "  ";
        cout << lp.rowUpperBound(row) << endl;
    }


    cout << "\t";
    for(Lp::ColIt col(lp); col!=INVALID; ++col){
        cout << "-" << "\t  ";
    }
    cout << "|  " << endl;


    cout << "\t";
    for(Lp::ColIt col(lp); col!=INVALID; ++col){
        cout << lp.objCoeff(col) << "\t  ";
    }
    cout << "|  " << endl;

    cout << endl;
}
