//利用三階精度迎風QUICK離散格式求解二維穩態熱擴散方程式 
//增加對角線兩端預設值的係數矩陣處裡
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
const double Pelect = 100.0 ;
const double  Gamma = (sqrt(2)/Pelect) ;//熱擴散係數 alpha
double dx ,dy  ;  // 修正：網格間距應為1/(NX-1)
const double u = 1.0;  // x方向速度
const double v = 1.0;  // y方向速度
const double T_l = 1.0 ;
const double T_b = 0.0 ;
const double lamda = 0.008 ;
int G, max_G = 10000000000;
double maxerror; 
const float tolerance = 1e-10;
bool steadystate;

//初始化矩陣
void initial(vector<vector<double> >& a, vector<double>& b , vector<double>& x) {
    // 初始化所有元素為0
    for(int i = 1; i <= n; i++) {
        for(int j = 1; j <= n; j++) {
        	if(i == 0 && j == 0){
        		a[i][j] = 1.0 ;
			}else if (i == n+1 && j == n+1){
			
                a[i][j] = 1.0;
            }else{
        	    a[i][j] = 0.0 ;
		    }
		}
        b[i] = 0.0;
        x[i] = 0.0;
    }
    ///////////////////////////////////////////////////////////////
    // 角點處理
    //第一個點T_{0,0} -> a[1][1] (左下角)
    a[1][1] = 6.0 * Gamma + dx * (7.0/8.0) * u + dy * (7.0/8.0) * v ;
    a[1][2] = (-1.0) * Gamma + dy * (3.0/8.0) * u;  // T_{1,0}
    a[1][NX+1] = (-1.0) * Gamma + dx * (3.0/8.0) * v;  // T_{0,1}
    b[1] = T_l * (2.0 * Gamma + dy * (5.0/4.0) * u ) + T_b * (2.0 * Gamma + dx * (4.0/5.0) * v);  // 左邊界=1.0

    //第二個點T_{1,0} -> a[2][2] (下邊界第二排)
    a[2][2] = 5.0 * Gamma + dy * (3.0/8.0) * u + dx * (7.0/8.0) * v ;
    a[2][1] = -1.0 * Gamma - dy * u ;  // T_{0,0}
    a[2][3] = -1.0 * Gamma + dy * (3.0/8.0) * u ;  // T_{2,0}
    a[2][NX+2] = (-1.0) * Gamma + dx * (3.0/8.0) * v ;  // T_{1,1}
    b[2] = -dy * (1.0/4.0) * u * T_l + (2.0 * Gamma + dx * (10.0/8.0) * u) * T_b ;  // 下邊界=0.0

    //第三個點T_{NX-1,0} -> a[NX][NX] (右下角)
    a[NX][NX] = 4.0 * Gamma + dy * (5.0/8.0) * u + dx * (7.0/8.0) * v ;
    a[NX][NX-1] = -Gamma + dy * (-6.0/8.0) * u ;  // T_{NX-2,0}
    a[NX][NX-2] = dy * (1.0/8.0) * u ;  // T_{NX-3,0}
    a[NX][2*NX] = (-1.0) * Gamma + dx * (3.0/8.0) * v ;  // T_{NX-1,1}
    b[NX] =((10.0/8.0)*dx+Gamma*(2.0))* T_b ;  // 右邊界和下邊界都=0.0

    //第四個點T_{0,1} -> a[NX+1][NX+1] (左邊界第二排)
    a[NX+1][NX+1] = 5.0 * Gamma + dy * (7.0/8.0) * u  + dx * (3.0/8.0) * v ;
    a[NX+1][NX+2] = (-1.0) * Gamma + dy * (3.0/8.0) * u ;  // T_{1,1}
    a[NX+1][1] = -1.0 * Gamma - dx * v ;  // T_{0,0}
    a[NX+1][2*NX+1] = -1.0 * Gamma + dx * (3.0/8.0) * v ;  // T_{0,2}
    b[NX+1] = ((2.0) * Gamma + dy * (10.0/8.0) * u) * T_l + (-dx * (2.0/8.0) * v)*T_b ;

    //第五個點T_{1,1} -> a[NX+2][NX+2] (內部第二排)
    a[NX+2][NX+2] = 4.0 * Gamma + dy * (3.0/8.0) * u + dx * (3.0/8.0) * v;
    a[NX+2][NX+1] = -1.0 * Gamma - dy * u ;  // T_{0,1}
    a[NX+2][NX+3] = -1.0 * Gamma + dy * (3.0/8.0) * u ;  // T_{2,1}
    a[NX+2][2] = -1.0 * Gamma - dx * v ;  // T_{1,0}
    a[NX+2][2*NX+2] = -1.0 * Gamma + dx * (3.0/8.0) * v;  // T_{1,2}
    b[NX+2] = -dy * (2.0/8.0) * u * T_l + (-dx * (2.0/8.0) * v) * T_b;  // 左邊界影響

    // 繼續其他角點...
    //第六個點T_{NX-1,1} -> a[2*NX][2*NX] (右邊界第二排)
    a[2*NX][2*NX] = 3*Gamma + dy * (5.0/8.0) * u  + dx * (3.0/8.0) * v ;
    a[2*NX][2*NX-1] = -Gamma + dy * (-6.0/8.0) * u;  // T_{NX-2,1}
    a[2*NX][2*NX-2] = dy * (1.0/8.0) * u ;  // T_{NX-3,1}
    a[2*NX][NX] = -1.0 * Gamma - dx * v ;  // T_{NX-1,0}
    a[2*NX][3*NX] = -1.0 * Gamma + dx * (3.0/8.0) * v ;  // T_{NX-1,2}
    b[2*NX] = dx*(-(2.0/8.0))*T_b ;

    //第七個點T_{0,NY-1} -> a[n-NX+1][n-NX+1] (左上角)
    a[n-NX+1][n-NX+1] = 4*Gamma + dy*(7.0/8.0) + dx*(5.0/8.0);
    a[n-NX+1][n-NX+2] = (-1.0)*Gamma+(3.0/8.0)*dy ;  // T_{1,NY-1}
    a[n-NX+1][n-2*NX+1] = -Gamma+(-6.0/8.0)*dx;  // T_{0,NY-2}
    a[n-NX+1][n-3*NX+1] = (1.0/8.0)*dx;  // T_{0,NY-3}
    b[n-NX+1] = T_l*((5.0/4.0)*dy+(2.0)*Gamma);  // 左邊界和上邊界都=1.0

    //第八個點T_{1,NY-1} -> a[n-NX+2][n-NX+2] (上邊界第二排)
    a[n-NX+2][n-NX+2] = 3*Gamma + (3.0/8.0)*dy + (5.0/8.0)*dx ;
    a[n-NX+2][n-NX+1] = -Gamma + (-8.0/8.0)*dy ;
    a[n-NX+2][n-NX+3] = -Gamma + (3.0/8.0)*dy ;
    a[n-NX+2][n-2*NX+2] = -Gamma + (-6.0/8.0)*dx;
    a[n-NX+2][n-3*NX+2] = dx*(1.0/8.0);
    b[n-NX+2] = T_l * (-1.0/4.0)*dy ;

    //第九個點T_{NX-1,NY-1} -> a[n][n] (右上角)
    a[n][n] = 2*Gamma+(5.0/8.0)*dy+(5.0/8.0)*dx ;
    a[n][n-1] = -Gamma+(-6.0/8.0)*dy;  // T_{NX-2,NY-1}
    a[n][n-2] = (1.0/8.0)*dy;  // T_{NX-3,NY-1}
    a[n][n-NX] =-Gamma+(-6.0/8.0)*dx;  // T_{NX-1,NY-2}
    a[n][n-2*NX] = (1.0/8.0)*dx;  // T_{NX-1,NY-3}
    b[n] = 0.0;

    /////////////////////////////////////////////////////////////////////
    // 邊界處理
    // 下邊界 (除角點外)
    for(int i = 3; i <= NX-1; i++) {
        a[i][i] = 5* Gamma + dy * (3.0/8.0) * u + dx * (7.0/8.0) * v;
        a[i][i-1] = -1.0 * Gamma - dy * (7.0/8.0) * u;
        a[i][i+1] = -1.0 * Gamma + dy * (3.0/8.0) * u;
        a[i][i-2] = dy * (1.0/8.0) * u;
        a[i][i+NX] = (-1) * Gamma + dx * (3.0/8.0) * v;
        b[i] = 2.0 * Gamma * T_b + dx * (10.0/8.0) * v * T_b ;  // 下邊界=0.0
    }
    
    // 上邊界 (除角點外)
    for(int i = n-NX+3; i <= n-1 ; i++) {
        a[i][i] = 3*Gamma + (3.0/8.0)*dy + (5.0/8.0)*dx;
        a[i][i-1] = -Gamma + dy*(-7.8/8.0);//西計算點 
        a[i][i+1] = -Gamma + dy*(3.0/8.0) ;//東計算點 
        a[i][i-2] = (1.0/8.0)*dy ;//西西計算點 
        a[i][i-NX] = -Gamma + dx*(-6.0/8.0);//南計算點 
        a[i][i-2*NX] = dx*(1.0/8.0 );//南南計算點 
        b[i] = 0.0 ;  // 上邊界=1.0
    }
    
    // 左邊界 (除角點外)
    for(int i = 3 ; i <= NX-1 ; i++) {
        int idx = NX * (i-1) + 1;
        a[idx][idx] = 5.0 * Gamma + dy * (7.0/8.0) * u  + dx * (3.0/8.0) * v ;
        a[idx][idx+1] = (-1.0) * Gamma + dy * (3.0/8.0) * u ;
        a[idx][idx-NX] = -1.0 * Gamma - dx * (7.0/8.0) * v ;
        a[idx][idx+NX] = -1.0 * Gamma + dx * (3.0/8.0) * v ;
        a[idx][idx-2*NX] = dx * (1.0/8.0) * v ;
        b[idx] = (2.0) * Gamma * T_l + dy * (10.0/8.0) * u * T_l ;  // 左邊界=1.0
    }
    
    // 右邊界 (除角點外)
    for(int i = 3; i <= NX-1; i++) {
        int idx = NX * i;
        a[idx][idx] = 3.0 * Gamma + dy * (5.0/8.0) * u + dx * (3.0/8.0) * v ;
        a[idx][idx-1] = -Gamma + dy * (-6.0/8.0) * u;
        a[idx][idx-2] = dy * (1.0/8.0) * u ;
        a[idx][idx-NX] = -1.0 * Gamma + dx * (-7.0/8.0) * v ;
        a[idx][idx+NX] = -1.0 * Gamma + dx * (3.0/8.0) * v ;
        a[idx][idx-2*NX] = dx * (1.0/8.0) * v ;
        b[idx] = 0.0 ;  // 右邊界=0.0
    }
    
    ////////////////////////////////////////////////////////////////////
    // 第二排邊界
    // 左邊界第二排 (除角點外)
    for(int i = 2; i <= NX-2; i++) {
        int idx = NX * i + 2;
        a[idx][idx] = 4.0 * Gamma + dy * (3.0/8.0) * u + dx * (3.0/8.0) * v;
        a[idx][idx-1] = -1.0 * Gamma - dy * u ;
        a[idx][idx+1] = -1.0 * Gamma + dy * (3.0/8.0) * u ;
        a[idx][idx-NX] = -1.0 * Gamma - dx * (7.0/8.0) * v ;
        a[idx][idx+NX] = -1.0 * Gamma + dx * (3.0/8.0) * v ;
        a[idx][idx-2*NX] = dx * (1.0/8.0) * v;
        b[idx] = -dy * (2.0/8.0) * u * T_l ;  // 左邊界影響
    }

    // 下邊界第二排 (除角點外)
    for(int i = 3; i <= NX-1; i++) {
        int k = i + NX;
        a[k][k] = 4.0 * Gamma + dy * (3.0/8.0) * u + dx * (3.0/8.0) * v ;
        a[k][k-1] = -1.0 * Gamma - dy * (7.0/8.0) * u ;
        a[k][k+1] = -1.0 * Gamma + dy * (3.0/8.0) * u ;
        a[k][k-2] = dy * (1.0/8.0) * u ;
        a[k][k-NX] = -1.0 * Gamma - dx * v ;
        a[k][k+NX] = -1.0 * Gamma + dx * (3.0/8.0) * v ;
        b[k] = -dx * (2.0/8.0) * v * T_b ;  // 下邊界影響
    }
    
    /////////////////////////////////////////////////////////////////////
    // 內點
    for(int j = 3; j <= NX-1; j++) {
        for(int i = 3; i <= NX-1; i++) {
            int idx = (j-1) * NX + i;
            a[idx][idx] = 4.0 * Gamma + dy * (3.0/8.0) * u + dx * (3.0/8.0) * v;
            a[idx][idx-1] = -1.0 * Gamma - dy * (7.0/8.0) * u ;  // 西點
            a[idx][idx+1] = -1.0 * Gamma + dy * (3.0/8.0) * u ;  // 東點
            a[idx][idx-2] = dy * (1.0/8.0) * u ;  // 西西點
            a[idx][idx-NX] = -1.0 * Gamma - dx * (7.0/8.0) * v ;  // 南點
            a[idx][idx+NX] = -1.0 * Gamma + dx * (3.0/8.0) * v ;  // 北點
            a[idx][idx-2*NX] = dx * (1.0/8.0) * v ;  // 南南點
            b[idx] = 0.0;  // 內點無熱源
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
    name << "QUICK_2JC(BW)0925_Pe = " << Pelect <<","<<"Grid "<< nodes_i << "x" << nodes_j << ","<< setfill('0') << setw(6) << m << ".vtk";
    ofstream out(name.str().c_str());
    
    // VTK 文件頭：使用 STRUCTURED_GRID
    out << "# vtk DataFile Version 3.0\n";
    out << "QUICK_2JC(BW)0925\n";
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
        cout << "Boundary condition : up & right : Neumann n dot Gradient of T = 0 " << endl;
        cout << "Boundary condition : left = 1.0,  and bottom = 0.0" << endl;


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