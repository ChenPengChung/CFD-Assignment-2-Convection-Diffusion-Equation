//利用有限體積法求解不可壓縮理想氣體的二維穩態熱擴散對流方程式 
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>
#include <iomanip>
using namespace std;
//會用到的參數
//會用到的參數
int Nx[] = {40,80};
//有限體積法採用Wet node處理，故計算點數目 = 網格數目, 節點數目 = 網格數目+1
int NX; //NX = NY => dx = dy
int n;
vector<vector<double> > a; //a二維矩陣為等式左側係數矩陣
vector<double> b; //b一維矩陣為等式右側之外源項，方程組之非齊性項
vector<double> x;
vector<double> x_old;
vector<double> x_old2;
vector<double> x_nuse ;
vector<vector<double> > T;
const double Pelect = 1000.0 ;
const double  Gamma = (sqrt(1.0)/Pelect) ;//熱擴散係數 alpha
double dx ,dy  ;  // 修正：網格間距應為1/(NX-1)
const double lamda = 0.01 ;
int G, max_G = 10000000000;
double maxerror; 
const float tolerance = 1e-10;
bool steadystate;


//初始化矩陣
void initial(vector<vector<double> >& a, vector<double>& b , vector<double>& x) {
    // 初始化所有元素為0
    for(int i = 1; i <= n; i++) {
        for(int j = 1; j <= n; j++) {
            a[i][j] = 0.0;
        }
        b[i] = 0.0;
        x[i] = 0.0;
        x_old[i] = 0.0;
    }
    
    // 設定邊界條件和係數矩陣
    // 四個角點
    a[1][1] = (6.0*Gamma) + (dx/2.0) + (dy/2.0) ; 
    a[1][2] = -Gamma + (dy/2.0);
    a[1][NX+1] = -Gamma + (dx/2.0);
    b[1] = 2.0*Gamma + dy ; // 邊界溫度
    
    a[NX][NX] = (4.0*Gamma) + (dx/2.0) + (dy/2.0);
    a[NX][NX-1] = -Gamma - (dy/2.0);
    a[NX][2*NX] = -Gamma + (dy/2.0);
    b[NX] = 0.0;// 邊界溫度
    
    a[n-NX+1][n-NX+1] = (4.0*Gamma) + (dx/2.0) + (dy/2.0);
    a[n-NX+1][n-NX+2] =  -Gamma + (dy/2.0);
    a[n-NX+1][n-2*NX+1] = -Gamma - (dx/2.0);
    b[n-NX+1] = 2.0*Gamma + dy ; //邊界效應引起之外源項
    
    a[n][n] =(2.0*Gamma) + (dx/2.0) + (dy/2.0);
    a[n][n-1] = -Gamma - (dy/2.0) ;
    a[n][n-NX] = -Gamma - (dx/2.0);
    b[n] = 0.0 ;// 邊界溫度
    
    // 下邊界 (除角點外)
    for(int i = 2; i <= NX-1; i++) {
        a[i][i] = (5*Gamma)+(dx/2.0);
        a[i][i-1] =-(Gamma + (dy/2.0));
        a[i][i+1] = -(Gamma - (dy/2.0));
        a[i][i+NX] = - (Gamma - (dx/2.0));
        b[i] = 0.0; // 邊界溫度
    }
    
    // 上邊界 (除角點外)
    for(int i = n-NX+2; i < n; i++) {
        a[i][i] = (3.0*Gamma)+(dx/2.0);
        a[i][i+1] = -Gamma + (dy/2.0);
        a[i][i-1] = -Gamma - (dy/2.0);
        a[i][i-NX] = -Gamma - (dx/2.0);
        b[i] =0.0 ;
    }
    
    // 左邊界 (除角點外)
    for(int i = 1; i <= NX-2; i++) {
        int idx = NX*i+1;
        a[idx][idx] = (5*Gamma)+(dy/2.0);
        a[idx][idx+1] = -(-(dy/2.0) + Gamma) ;
        a[idx][idx-NX] = -((dx/2.0) + Gamma);
        a[idx][idx+NX] = -(-(dx/2.0) + Gamma);
        b[idx] = dy + 2.0 * Gamma ; //邊界溫度
    }
    
    // 右邊界 (除角點外)
    for(int i = 2; i <= NX-1; i++) {
        int idx = NX*i;
        a[idx][idx] = (3*Gamma)+(dy/2.0);
        a[idx][idx-1] = -((dy/2.0)+Gamma);
        a[idx][idx-NX] = -((dx/2.0)+Gamma);
        a[idx][idx+NX] = (dx/2.0)-Gamma;
        b[idx] = 0.0;
    }
    
    // 內點
    for(int i = 2; i <= NX-1; i++) {
        for(int j = 1; j <= NX-2; j++) {
            int idx = j*NX+i;
            a[idx][idx] = (4*Gamma);
            a[idx][idx+1] = -(Gamma-(dy/2.0));
            a[idx][idx-1] = -(Gamma+(dy/2.0));
            a[idx][idx-NX] = -(Gamma+(dx/2.0));
            a[idx][idx+NX] = -(Gamma-(dx/2.0));
            b[idx] = 0.0; // 內點無熱源
        }
    }
}
void Jacobi(vector<vector<double> >& a, vector<double>& b , vector<double>& x , int n) {
	
	if(G == 0){//經過0次迭代,要進行第一次迭代 
	// 先複製當前解到x_old
	//以及賦值解到前二解 //112
        for(int k = 1; k <= n; k++) {
            x_old2[k] = x[k];
            x_old[k] = x[k];
        }
	}else if (G == 1){//經過一次迭代，要進行第二次迭代 
	    // 先複製當前解到x_old
        for(int k = 1; k <= n; k++) {
            //前前解不動。 //123
            x_old[k] = x[k];
        }
	}else{//234//345//456.....
		for(int k = 1; k <= n; k++) {
            x_old2[k] = x_old[k] ;
            x_old[k] = x[k] ;
        }
	}
    // 計算新的解
    for(int k = 1; k <= n; k++) {
        double sum = 0;
        for(int p = 1; p <= n; p++) {
            if(p != k) {
                sum += a[k][p] * x_old[p];
            }
        }
        x[k] = ((b[k] - sum) / a[k][k]);
        //產生新的溫度場 
	    x_nuse[k] = x[k] ;//未使用鬆弛技術之新迭代溫度場 
        //under relaxtion(欠鬆弛迭代) && over relation(超鬆弛迭代)
        // 計算最大誤差
        if(fabs(x_nuse[k] - x_old[k]) > fabs (x_old[k] - x_old2[k]) ){
    	    //才能採用欠鬆弛技術
    	    x[k] = x_old[k] + lamda * (x_nuse[k] - x_old[k]) ;
	    }
    }
    maxerror = 0;
    for(int k = 1; k <= n; k++) {
        double error = fabs(x[k] - x_old[k]);
        if(maxerror < error) {
            maxerror = error;
        }
    }
}

void output(int m) {
    // 將一維解轉換為二維溫度場
    for(int j = 0; j < NX; j++) {
        for(int i = 0; i < NX; i++) {
            T[i][j] = x[j*NX + i + 1]; //輸入x[] : 1~n 
        }
    }
    
    int nodes_i = NX + 1;  // x方向節點數
    int nodes_j = NX + 1;  // y方向節點數
    int cells_i = NX;      // x方向網格數
    int cells_j = NX;      // y方向網格數
    
    ostringstream name;
    name << "CD_1JC_Pe = " << Pelect <<","<<"Grid "<< nodes_i << "x" << nodes_j << ","<< setfill('0') << setw(6) << m << ".vtk";
    ofstream out(name.str().c_str());
    
    // VTK 文件頭：使用 STRUCTURED_GRID
    out << "# vtk DataFile Version 3.0\n";
    out << "CD_1JC\n";
    out << "ASCII\n";
    out << "DATASET STRUCTURED_GRID\n";  // 改為 STRUCTURED_GRID
    out << "DIMENSIONS " << nodes_i << " " << nodes_j << " 1\n";
    
    // 輸出節點座標
    out << "POINTS " << nodes_i * nodes_j << " double\n";
    for(int j = 0; j < nodes_j; j++) {
        for(int i = 0; i < nodes_i; i++) {
            double x_coord = i * dx;
            double y_coord = j * dy;
            out << scientific << setprecision(6) << x_coord << " " << y_coord << " 0.0\n";
        }
    }
    
    // 輸出網格中心數據（CELL_DATA）
    out << "CELL_DATA " << cells_i * cells_j << "\n";  // 使用 CELL_DATA
    out << "SCALARS Temperature double 1\n";
    out << "LOOKUP_TABLE default\n";
    
    // 直接輸出網格中心點的溫度數據（你的 x[] 陣列）
    for(int j = 0; j < cells_j; j++) {
        for(int i = 0; i < cells_i; i++) {
            int idx = j * cells_i + i + 1;  // 對應你的 x[] 索引
            out << scientific << setprecision(6) << x[idx] << "\n";
        }
    }
    
    out.close();
    cout << "VTK document has been output: " << name.str() << endl;
}

int main() {
	for(int grid_idx = 0; grid_idx < 2; grid_idx++){
        NX = Nx[grid_idx];
        dx = 1.0/double(NX);
        dy = 1.0/double(NX);
        n = (NX)*(NX);
        // 重新調整向量大小
        a.assign(n+2, vector<double>(n+2, 0.0));
        b.assign(n+2, 0.0);
        x.assign(n+2, 0.0);
        x_old.assign(n+2, 0.0); //向量大小必須唯一 
        x_old2.assign(n+2, 0.0); //向量大小必須唯一 
        x_nuse.assign(n+2, 0.0); //向量大小必須唯一 
        T.assign(NX+2, vector<double>(NX+2, 0.0));
        cout << "C++ code is running quiuck 2025 ...." << endl;
        cout << "grid size is : " << NX+1 << " x " << NX+1 << endl;
        cout << "the mesh size : dx=" << dx << ", dy=" << dx << endl;
        cout << "Boundary condition : up&right : Neumann n dot Gradient of T = 0 " << endl;
        cout << "Boundary condition : left = 1.0, bottom = 0.0" << endl;


        steadystate = false;
        initial(a, b, x);
    
        for(G = 0; G < max_G; G++) {
            if(G > 0) {
                Jacobi(a, b, x, n);
            }
        if(G % 500 == 0) {
            cout << "itteration = " << G << endl;
            if(G > 0) {
                cout << "the max_change = " << scientific << maxerror << endl;
            }
            output(G);
        }
        
        if(G > 100 && maxerror < tolerance) {
            steadystate = true;
            cout << "will be steady state and temperature field converged!!" << endl;
            break;
        }
    }
    
        if(!steadystate) {
            cout << "\nWarning: Reached maximum iteration count, but not converged to steady state!" << endl;
        }
        
        output(G);  //輸出最終結果
        cout << "grid " << NX+1 << "x" << NX+1 << " computation completed\n" << endl;
    }
    cout << "All computations completed!" << endl;
    return 0;
}