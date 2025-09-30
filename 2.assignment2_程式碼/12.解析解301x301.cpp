
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>
#include <iomanip>
using namespace std;
double T[301][301] ;

void output() {
    for(int j = 1; j <= 301; j++){
        for(int i = j; i <= 301; i++){
            T[i-1][j-1] = 0.0 ;
        }
    }
    for(int j = 1; j <= 301; j++){
        for(int i = 1; i < j; i++){
            T[i-1][j-1] = 1.0 ;
        }
    }
    
    ostringstream name;
    name << "解析解_" << 301 << "x" << 301 << ".vtk";
    ofstream out(name.str().c_str());
    out << "# vtk DataFile Version 3.0\n";
    out << "解析解_\n";
    out << "ASCII\n";
    out << "DATASET STRUCTURED_POINTS\n";
    out << "DIMENSIONS " << 301 << " " << 301 << " 1\n";
    out << "ORIGIN 0.0 0.0 0.0\n";
    out << "SPACING " << 1/300.0 << " " << 1/300.0 << " 1.0\n";
    out << "POINT_DATA " << 301 * 301 << "\n";
    
    out << "SCALARS Temperature double 1\n";
    out << "LOOKUP_TABLE default\n";
    for(int j = 0; j < 301; j++) {
        for(int i = 0; i < 301; i++) {
            out << scientific << setprecision(6) << T[i][j] << "\n";
        }
    }
    
    out.close();
    cout << "VTK document output: " << name.str() << endl;
}

int main() {
	output() ;
    return 0;
}