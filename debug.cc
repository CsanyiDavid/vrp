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
	for(int i=0; i<11; ++i){
		cols.push_back(lp.addCol());
		lp.colLowerBound(cols[i], 0);
		if(i<9) lp.colUpperBound(cols[i], 1);
	}

	vector<Lp::Row> rows;
	for(int i=0; i<9; ++i){
		rows.push_back(lp.addRow(cols[i]>=1));
	}
	Lp::Expr e;
	for(int i=0; i<9; ++i){
		e+=cols[i];
	}
	e-=cols[9];
	rows.push_back(lp.addRow(e==0));
	Lp::Expr e2;
	e2=14*cols[0]
		+6*cols[1]
		+8*cols[2]
		+12*cols[3]
		+14*cols[4]
		+26*cols[5]
		+24*cols[6]
		+16*cols[7]
		+12*cols[8]
		-cols[10];
	rows.push_back(lp.addRow(e2==0));
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

	//lp.coeff(rows[0], cols[1], 3);

    cout << lp.solve() << endl;
	printLPMatrix(lp);

	//add new col
	cols.push_back(lp.addCol());
    lp.colLowerBound(cols[11], 0);
    lp.colUpperBound(cols[11], 1);

    lp.coeff(rows[4], cols[11], 1);
    lp.coeff(rows[5], cols[11], 1);
    lp.coeff(rows[6], cols[11], 1);
    lp.coeff(rows[7], cols[11], 1);
    lp.coeff(rows[8], cols[11], 1);
    lp.coeff(rows[9], cols[11], 1);
    lp.coeff(rows[10], cols[11], 36);
    lp.objCoeff(cols[11], 36);

    cout << lp.solve() << endl;
    lp.solve();

    if(lp.coeff(rows[3], cols[9])==0){
        cout << "OK" << endl;
    } else {
        cout << "ERROR" << endl;
    }

    printLPMatrix(lp);



    //add new col
    cols.push_back(lp.addCol());
    lp.colLowerBound(cols[12], 0);
    lp.colUpperBound(cols[12], 1);

    lp.coeff(rows[4], cols[12], 1);
    lp.coeff(rows[5], cols[12], 1);
    lp.coeff(rows[6], cols[12], 1);
    lp.coeff(rows[7], cols[12], 1);
    lp.coeff(rows[8], cols[12], 1);
    lp.coeff(rows[9], cols[12], 1);
    lp.coeff(rows[10], cols[12], 36);
    lp.objCoeff(cols[12], 36);

    cout << lp.solve() << endl;

    if(lp.coeff(rows[3], cols[9])==0){
        cout << "OK" << endl;
    } else {
        cout << "ERROR" << endl;
    }

    printLPMatrix(lp);
}

void printLPMatrix(const Lp& lp){
    cout << "PRINT LP MATRIX: " << endl;
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
