#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>
#include <iomanip>
#include <algorithm>
using namespace std;

//會用到的參數
int Nx[] = {80,40};
//有限體積法採用Wet node處理，故計算點數目 = 網格數目, 節點數目 = 網格數目+1
int NX , NY; //NX = NY => dx = dy
int n;

vector<double> x;
vector<double> x_old;
vector<vector<double> > T;
double dx, dy; //一個網格間距
//邊界條件定義
const double T_left = 1.0;    //左邊界Dirichlet條件
const double T_right = 0.0;   //右邊界Dirichlet條件  
const double T_bottom = 0.0;  //下邊界Dirichlet條件
//迭代參數
const double lamda = 0.05; //逐次超鬆弛迭代法

int G, max_G = 5000000;
const double w = 1e-6 ;
double maxerror ;
const float tolerance = 1e-10;
bool steadystate;
double MUSCL_w(int i,int j){
	double T_D = T[i][j] ;
	double T_U = T[i-1][j] ;
	double T_UU = T[i-2][j];
	double x = (T_U-T_UU)/(T_D-T_U) ;
    if(x>= 0 && x<0.8) return (3+x)/4.0 ;  
    else if(x<0) return 0.0 ;
    else return 0.8 ;   
}
double MUSCL_e(int i,int j){
	double T_D = T[i+1][j] ;
	double T_U = T[i][j] ;
	double T_UU = T[i-1][j];
	double x = (T_U-T_UU)/(T_D-T_U) ;
	if(x>= 0 && x<0.8) return (3+x)/4.0 ;  
    else if(x<0) return 0.0 ;
    else return 0.8 ;   
}
double MUSCL_s(int i,int j){
	double T_D = T[i][j] ;
	double T_U = T[i][j-1] ;
	double T_UU = T[i][j-2];
	double x = (T_U-T_UU)/(T_D-T_U) ;
	if(x>= 0 && x<0.8) return (3+x)/4.0 ;  
    else if(x<0) return 0.0 ;
    else return 0.8 ;   
}
double MUSCL_n(int i,int j){
	double T_D = T[i][j+1] ;
	double T_U = T[i][j] ;
	double T_UU = T[i][j-1];
	double x = (T_U-T_UU)/(T_D-T_U) ;
	if(x>= 0 && x<0.8) return (3+x)/4.0 ;  
    else if(x<0) return 0.0 ;
    else return 0.8 ;   
}
//TVD-MUSCL格式為在東西方向為三點格式,南北方向為三點格式

//顯示迭代法（考慮Neumann邊界條件）
void Gamma0(vector<double>& x) { //因為係數式時變的,每一次迭代都必須要做.....
    if(G == 0){
    	for(int j = 1 ; j <= NX ; j++){
    		T[0][j]=  2*T_left-T[1][j] ;
		}
		for(int i = 1 ; i <= NX ; i++){
			T[i][0] = 2*T_bottom - T[i][1] ;
		}
	}
    int index1 = (1-1)*NX+1 ;
    x_old[index1] = x[index1] ;
    double x_new1 =  ((dx*T_left + dx*T_bottom) - x[index1+1]*dx*(MUSCL_e(1,1))/(2.0) - x[index1+NX]*dx*(MUSCL_n(1,1))/2.0)/ (dx*(1-(MUSCL_e(1,1))) + dx*(1-(MUSCL_n(1,1))/2.0)) ;
    x[index1] = x_old[index1] + lamda * (x_new1 - x_old[index1]);
    T[1][1] = x[index1] ;
    T[0][1] = 2*T_left-T[1][1] ;
    T[1][0] = 2*T_bottom-T[1][1] ;
    int index2 = (1-1)*NX+NX ;      
    x_old[index2] = x[index2] ;
    double x_new2 =  (dx*x[index2-1]*(1-(MUSCL_w(NX,1))/2.0)   -dx*x[index2+NX]*MUSCL_n(NX,1)/2.0      +dx*x[index2-NX]*(1-(MUSCL_s(NX,1)/2.0)))/(dx*(1-(MUSCL_w(NX,1))/2.0)+dx*(1-(MUSCL_n(NX,1))/2.0-(MUSCL_s(NX,1))/2.0)) ;
    x[index2] = x_old[index2] + lamda * (x_new2 - x_old[index2]);
    T[NX][1] = x[index2] ;
    T[NX+1][1] = 0.0 ;
    T[NX][0] = 2*T_bottom-T[NX][1] ;
    int index3 = (NX-1)*NX+1 ;
    x_old[index3] = x[index3] ;
    double x_new3 =  (dx*T_left - x[index3+1]*dx*(MUSCL_e(1,NX))/2.0 + x[index3-NX]*dx*(1-MUSCL_s(1,NX)/2.0) )/(dx*(1-(MUSCL_e(1,NX)/2.0)) + dx*(1-MUSCL_s(1,NX)/2.0)) ;
    x[index3] = x_old[index3] + lamda * (x_new3 - x_old[index3]);
    T[1][NX] = x[index3] ;
	T[0][NX] = 2*T_left-T[1][NX] ; 
	T[1][NX+1] = 0.0 ;   
    int index4 = (NX-1)*NX+NX ;        
    x_old[index4] = x[index4] ;
    double x_new4 =  (dx*x[index4-1]*(1-MUSCL_w(NX,NX)/2.0))+dx*x[index4-NX]*(1-MUSCL_s(NX,NX)/2.0)/(dx*(1-MUSCL_w(NX,NX)/2.0)+dx*(1-MUSCL_s(NX,NX)/2.0)) ;
    x[index4] = x_old[index4] + lamda * (x_new4 - x_old[index4]);
    T[NX][NX] = x[index4] ;
	T[NX][NX+1] = 0.0 ; 
	T[NX+1][NX] = 0.0 ;          
    // 邊界處理（不含角點）
    // 下邊界 (除角點外)(i = 2~NX-1 ; j = 1)
    for(int i = 2 ; i <= NX-1; i++) {
        int index = (1-1)*NX + i;
        x_old[index] = x[index] ;
        double x_new = (dx*T_bottom - x[index+1]*dx*MUSCL_e(i,1)/2.0 + x[index-1]*(dx*(1-MUSCL_w(i,1)/2.0)) -x[index+NX]*dx*MUSCL_n(i,1)/2.0 )/(dx*(1-MUSCL_e(i,1)/2.0-MUSCL_w(i,1)/2.0)+dx*(1-MUSCL_n(i,1)/2.0));
        x[index] = x_old[index] + lamda * (x_new - x_old[index]);
        T[i][1] = x[index] ;
        T[i][0] = 2*T_bottom-T[i][1] ;
    }
    //上邊界(i = 2~NX-1 ; j = NX)
    for(int i = 2 ; i <= NX-1; i++) {
    int index = (NX-1)*NX + i;
        x_old[index] = x[index] ;
        double x_new = (0.0 - x[index+1]*dx*MUSCL_e(i,NX)/2.0 + x[index-1]*(dx*(1-MUSCL_w(i,NX)/2.0)) + x[index-NX]*(dx*(1-MUSCL_s(i,NX)/2.0)))/(dx*(1-MUSCL_e(i,NX)/2.0-MUSCL_w(i,NX)/2.0)+dx*(1-MUSCL_s(i,NX)/2.0)) ;
        x[index] = x_old[index] + lamda * (x_new - x_old[index]);
        T[i][NX] = x[index] ;
        T[i][NX+1] = 0.0 ;
    }

    // 左邊界 (除角點外)
    for(int j = 2; j <= NX-1; j++) {
        int index = (j-1)*NX + 1;
        x_old[index] = x[index] ;
        double x_new = (dx*T_left - x[index+1]*dx*MUSCL_e(1,j)/2.0 - x[index+NX]*dx*MUSCL_n(1,j)/2.0 +x[index-NX]*(dx*(1-MUSCL_s(1,j)/2.0)) )/(dx*(1-MUSCL_e(1,j)/2.0)+dx*(1-MUSCL_n(1,j)/2.0-MUSCL_s(1,j)/2.0)) ;
        x[index] = x_old[index] + lamda * (x_new - x_old[index]);
        T[1][j] = x[index] ;
		T[0][j] = 2*T_left-T[1][j] ; 
    }
   
    // 右邊界 (除角點外)
    for(int j = 2; j <= NX-1; j++) {
        int index = (j-1)*NX + NX;
        x_old[index] = x[index] ;
        double x_new = (dx*x[index-1]*(1-(MUSCL_w(NX,j))/2.0)-dx*x[index+NX]*(MUSCL_n(NX,j)/2.0) +dx*x[index-NX]*(1-(MUSCL_s(NX,j))/2.0))/(dx*(1-(MUSCL_w(NX,j))/2.0)+dx*(1-(MUSCL_n(NX,j))/2.0-(MUSCL_s(NX,j))/2.0)) ;
        x[index] = x_old[index] + lamda * (x_new - x_old[index]);
        T[NX][j] = x[index] ;
        T[NX+1][j] = 0.0 ;
    }
    // 內點
    for(int i = 2; i <= NX-1; i++) {
        for(int j = 2; j <= NX-1; j++) {
            int index = (j-1)*NX+i;
            x_old[index] = x[index] ;
            double x_new = (0.0 - x[index+1]*dx*MUSCL_e(i,j)/2.0 + x[index-1]*(dx*(1-MUSCL_w(i,j)/2.0)) - x[index+NX]*dx*MUSCL_n(i,j)/2.0 + x[index-NX]*(dx*(1-MUSCL_s(i,j)/2.0)))/(dx*(1-MUSCL_e(i,j)/2.0-MUSCL_w(i,j)/2.0)+dx*(1-MUSCL_n(i,j)/2.0-MUSCL_s(i,j)/2.0)) ;
            x[index] = x_old[index] + lamda * (x_new - x_old[index]);
            T[i][j] = x[index] ;
        }
    }
}

void MaxError(vector<double>& x, int n) {
    //每一次迭代完成都更新溫度場
    //計算當前步最大誤差
    maxerror = 0;
    for(int k = 1; k <= n; k++) {
        double error = fabs(x[k] - x_old[k]);
        if(maxerror < error) {
            maxerror = error;
        }
    }//計算i = 1 : NX ; j = 1 : NX-1 的最大誤差
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
    name << "19.TVD_MUSCL_Pe = infintie" << "," << "Grid "<< nodes_i << "x" << nodes_j << ","<< setfill('0') << setw(6) << m << ".vtk";
    ofstream out(name.str().c_str());
    
    // VTK 文件頭：使用 STRUCTURED_GRID
    out << "# vtk DataFile Version 3.0\n";
    out << "19.TVD_MUSCL_Pe = infintie" << "," << "Grid "<< nodes_i << "x" << nodes_j << ","<< setfill('0') << setw(6) << m << ".vtk\n";
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
        dy = dx ;
        NY = NX ; 
        n = (NX)*(NX);
        // 重新調整向量大小
        x.assign(n+2, 0.0);
        x_old.assign(n+2, 0.0); //向量大小必須唯一
        T.assign(NX+2, vector<double>(NX+2, 0.0));
        cout << "========================================" << endl;
        cout << "TVD_MUSCL Pe = infinite is coming ...." << endl;
        cout << "mesh size : " << NX+1 << " x " << NX+1 << endl;
        cout << "Grid space: dx = dy = " << dx << endl;
        cout << "Boundary conditions:" << endl;
        cout << "  Left boundary (Dirichlet): T = " << T_left << endl;
        cout << "  Right boundary (Dirichlet): T = " << T_right << endl;
        cout << "  Bottom boundary (Dirichlet): T = " << T_bottom << endl;
        cout << "  Top boundary (Neumann):  T_[i][NY] = T_[i][NY+1] (virtual point)" << endl;
        cout << "Relaxation factor: λ = " << lamda << endl;
        cout << "Convergence criterion: " << tolerance << endl;
        cout << "========================================" << endl;
        steadystate = false;
        for(int i = 1 ; i <= n ; i++){
        x[i] = 0.0 ;
        }
       x[0] = 0.0 ;
       x[n+1] = 0.0 ;
        for(G = 0; G < max_G; G++) {
            Gamma0(x);  //設定係數矩陣（含邊界條件）
            MaxError(x, n);  //計算最大誤差
            if(G % 1000000 == 0) {
                cout << "Itteration number = " << G;
                if(G > 0) {
                    cout << ", Max error = " << scientific << maxerror;
                }
                cout << endl;
                
            }
           
            if(G > 1000 && maxerror < tolerance) {
                steadystate = true;
                cout << "\nHave achieved steady state!" << endl;
                cout << "Final iteration: " << G << endl;
                cout << "Final error: " << scientific << maxerror << endl;
                break;
            }
        }
       
        if(!steadystate) {
            cout << "\nWarning: Reached maximum iterations, but did not achieve steady state!" << endl;
        }
       
        output(G);  //輸出最終結果
        cout << "Grid " << NX+1 << "x" << NX+1 << " computation completed\n" << endl;
    }
    cout << "All computations completed!" << endl;
    return 0;
}