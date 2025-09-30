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
const int NX = 81;
const int NY = 81;
const int n = NX * NY;
double a[n+1][n+1], b[n+1]; //a二維矩陣為等式左側係數矩陣;b一維矩陣為等式右側之外源項，方程組之非齊性項
double x[n+1], x_old[n+1], T[NX][NY];
const double Pelect = 10 ;
const double Gamma = (sqrt(2) /Pelect);//熱擴散係數
const double dx = 1.0/double(NX-1) ;
const double dy = 1.0/double(NY-1) ;
const double a_W = (Gamma + dy/2.0);
const double a_E = (Gamma - dy/2.0);
const double a_S = (Gamma + dx/2.0);
const double a_N = (Gamma - dx/2.0);
int G, max_G = 10000; //限制最大迭代次數
double maxerror; 
const double tolerance = 1e-6;
bool steadystate;

//初始化矩陣
void initial(double a[][n+1], double b[], int n) {
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
    a[1][2] = -a_E;
    a[1][NX+1] = -a_N;
    b[1] = 2.0*a_W; // 邊界溫度
    
    a[NX][NX] = (6.0*Gamma) - (dx/2.0) + (dy/2.0);
    a[NX][NX-1] = -a_W;
    a[NX][2*NX] = -a_N;
    b[NX] = 0.0;// 邊界溫度
    
    a[n-NX+1][n-NX+1] = (6.0*Gamma) + (dx/2.0) - (dy/2.0);
    a[n-NX+1][n-NX+2] =  -a_E;
    a[n-NX+1][n-2*NX+1] = -a_S;
    b[n-NX+1] = 2.0*(a_W+a_N); //邊界效應引起之外源項
    
    a[n][n] =(6.0*Gamma) - (dx/2.0) - (dy/2.0);
    a[n][n-1] = -a_W;
    a[n][n-NX] = -a_S;
    b[n] = 2.0*a_N;// 邊界溫度
    
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
        a[i][i] = (5*Gamma)-(dx/2.0);
        a[i][i+1] = -(Gamma + (dy/2.0));
        a[i][i-1] = -(Gamma - (dy/2.0));
        a[i][i-NX] = -(Gamma + (dx/2.0));
        b[i] =2.0*Gamma- dx;
    }
    
    // 左邊界 (除角點外)
    for(int i = 1; i <= NY-2; i++) {
        int idx = NX*i+1;
        a[idx][idx] = (5*Gamma)+(dy/2.0);
        a[idx][idx+1] = -(-(dy/2.0) + Gamma) ;
        a[idx][idx-NX] = -((dx/2.0) + Gamma);
        a[idx][idx+NX] = -(-(dx/2.0) + Gamma);
        b[idx] = dy + 2.0 * Gamma ; //邊界溫度
    }
    
    // 右邊界 (除角點外)
    for(int i = 2; i <= NY-1; i++) {
        int idx = NX*i;
        a[idx][idx] = (5*Gamma)-(dy/2.0);
        a[idx][idx-1] = -((dy/2.0)+Gamma);
        a[idx][idx-NX] = -((dx/2.0)+Gamma);
        a[idx][idx+NX] = -(-(dx/2.0)+Gamma);
        b[idx] = 0.0;
    }
    
    // 內點
    for(int i = 2; i <= NX-1; i++) {
        for(int j = 1; j <= NY-2; j++) {
            int idx = j*NX+i;
            a[idx][idx] = (4*Gamma);
            a[idx][idx+1] = -(Gamma+(dy/2.0));
            a[idx][idx-1] = -(Gamma-(dy/2.0));
            a[idx][idx-NX] = -(Gamma+(dx/2.0));
            a[idx][idx+NX] = -(Gamma-(dx/2.0));
            b[idx] = 0.0; // 內點無熱源
        }
    }
}

void GaussSeidel(double a[][n+1], double b[], double x[], int n) {
    // 先複製當前解到x_old用於計算誤差
    for(int k = 1; k <= n; k++) {
        x_old[k] = x[k];
    }
    
    // 計算新的解 - 關鍵差異：立即使用新計算的值
    for(int k = 1; k <= n; k++) {
    	double sum2 = 0.0 ;
    	for(int i = 1 ; i < k ; i++){ //迭加以更新的x[k] 值 
    		sum2 = sum2 + a[k][i] * x[i] ;
		}
    	
        double sum = 0;
        
        // 對於尚未更新的值(j > k)，使用舊的 x[p]
        for(int p = k + 1; p <= n; p++) {
            sum += a[k][p] * x[p];  // 使用舊值
        }
        
        // 計算新的 x[k]
        x[k] = (b[k] - sum2 - sum) / a[k][k];
        
        
    }
    
    // 計算最大誤差
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
    name << "CD_2GS_Pe = " << Pelect <<","<<"Grid "<< nodes_i << "x" << nodes_j << ","<< setfill('0') << setw(6) << m << ".vtk";
    ofstream out(name.str().c_str());
    
    // VTK 文件頭：使用 STRUCTURED_GRID
    out << "# vtk DataFile Version 3.0\n";
    out << "CD_2GS\n";
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
    cout << "程式碼開始執行...." << endl;
    cout << "網格大小: " << NX << " x " << NY << endl;
    cout << "網格間距: dx=" << dx << ", dy=" << dy << endl;
    cout << "邊界條件: 左邊界和上邊界=1.0, 右邊界和下邊界=0.0" << endl;
    cout << "Pelect number = " << Pelect << ",Gamma = " << Gamma << endl ;
    
    steadystate = false;
    initial(a, b, n);
    
    for(G = 0; G < max_G; G++) {
        if(G > 0) {
            GaussSeidel(a, b, x, n);
        }

        if(G % 1000 == 0) {
            cout << "迭代次數 = " << G << endl;
            if(G > 0) {
                cout << "最大變化量max_change = " << scientific << maxerror << endl;
            }
            output(G);
        }
        
        if(G > 100 && maxerror < tolerance) {
            steadystate = true;
            cout << "已經達到迭代終點，溫度場收斂!!" << endl;
            break;
        }
    }
    
    if(!steadystate) {
        cout << "達到最大迭代次數，但未達到穩態!" << endl;
    }
    
    output(G);
    cout << "格子大小 " << NX << "x" << NY << " 計算完成\n" << endl;
    return 0;
}