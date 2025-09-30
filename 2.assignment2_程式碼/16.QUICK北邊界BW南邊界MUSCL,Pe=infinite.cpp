//利用QUICK格式求解二維穩態純對流方程式
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
int Nx[] = {40,80};
//有限體積法採用Wet node處理，故計算點數目 = 網格數目, 節點數目 = 網格數目+1
int NX; //NX = NY => dx = dy
int n;
vector<vector<double> > a; //a二維矩陣為等式左側係數矩陣
vector<double> b; //b一維矩陣為等式右側之外源項，方程組之非齊性項
vector<double> x;
vector<double> x_old;
vector<vector<double> > T;
vector<double> r_s;
double dx; //一個網格間距
//邊界條件定義
const double T_left = 1.0;    //左邊界Dirichlet條件
const double T_right = 0.0;   //右邊界Dirichlet條件  
const double T_bottom = 0.0;  //下邊界Dirichlet條件
//上邊界為Neumann條件: ?T/?n = 0
//QUICK格式係數
double a1,a2,a3,a4,a5;
//迭代參數
const double lamda = 0.3; //逐次超鬆弛迭代法
int G, max_G = 10000000;
double maxerror ; 
const float tolerance = 1e-10;
bool steadystate;


double D(double x){
	if(x<0)return  0.0 ;
	else if(x >= 0 && x < (1.0/3.0)) return 0.0 ;
	else if(x >= (1.0/3.0) && x < 3.0) return 0.25 ;
	else return 1.0 ;
}
double U(double x){
	if(x<0)return 1.0 ;
	else if(x >= 0 && x < (1.0/3.0)) return  2.0 ;
	else if(x >= (1.0/3.0) && x < 3.0) return 1.0 ;
	else return 0.0 ;
}
double UU(double x){
	if(x<0)return  0.0 ;
	else if(x >= 0 && x < (1.0/3.0)) return  -1.0 ;
	else if(x >= (1.0/3.0) && x < 3.0) return  -0.25 ;
	else return 0.0 ;
}

//係數矩陣賦值（考慮Neumann邊界條件）
void Gamma0(vector<vector<double> >& a, vector<double>& b , vector<double>& x) {
	// 初始化所有元素為0
    for(int i = 1; i <= NX*NX ; i++) {
        for(int j = 1; j <= NX*NX; j++) {
            a[i][j] = 0.0 ;
        }  
        b[i] = 0.0 ;
        x[i] = 0.0 ;
    }//記得不要每一次迭代都把 解賦值為0 
    a[0][0] = 0.0 ;
    a[NX*NX+1][NX*NX+1] = 1.0 ;
    ///////////////////////////////////////////////////////////////  
    //第一角點 (1,1) - 左下角
    int p1 = (1-1)*NX + 1;
    a[p1][p1] = a4+a4 ;     //本點
    a[p1][p1+1] = a2;      //東側計算點
    a[p1][p1+NX] = a2;     //北側計算點
    b[p1] = dx*(10.0/8.0)*(T_bottom+T_left);
    
    //第二角點 (1,2)
    int p2 = (1-1)*NX + 2;
    a[p2][p2] = a2+a4;     //本點
    a[p2][p2+1] = a2;      //東側計算點
    a[p2][p2-1] = -dx;     //西側計算點
    a[p2][p2+NX] = a2;     //北側計算點
    b[p2] = dx*(-2.0/8.0)*T_left + dx*(10.0/8.0)*T_bottom;
    
    //第三角點 (NX,1) - 右下角
    int p3 = (1-1)*NX + NX;
    a[p3][p3] = dx*(6.0/8.0) + dx*(7.0/8.0) ;    //本點
    a[p3][p3-1] = dx*(-7.0/8.0) ;     //西側計算點
    a[p3][p3-2] = dx*(1.0/8.0) ;      //西西側計算點
    a[p3][p3+NX] = (3.0/8.0)*dx ;      //北側計算點
    b[p3] = dx*(10.0/8.0)*T_bottom;
    ///////////////////////////////////////////////////////////////
    //第二排角點處理
    ///////////////////////////////////////////////////////////////
    
    //第四角點 (2,1)
    int p4 = (2-1)*NX + 1;
    a[p4][p4] = a4+a2;     //本點
    a[p4][p4+1] = a2;      //東側計算點
    a[p4][p4+NX] = a2;     //北側計算點
    a[p4][p4-NX] = -dx;    //南側計算點
    b[p4] = dx*(10.0/8.0)*T_left + dx*(-2.0/8.0)*T_bottom;
    
    //第五角點 (2,2)
    int p5 = (2-1)*NX + 2;
    a[p5][p5] = a2+a2;     //本點
    a[p5][p5+1] = a2;      //東側計算點
    a[p5][p5-1] = -dx;     //西側計算點
    a[p5][p5+NX] = a2;     //北側計算點
    a[p5][p5-NX] = -dx;    //南側計算點
    b[p5] = dx*(-2.0/8.0)*T_left + dx*(-2.0/8.0)*T_bottom;
    
    //第六角點 (2,NX)
    int p6 = (2-1)*NX + NX;
    a[p6][p6] = (6.0/8.0)*dx + (3.0/8.0)*dx ; //本點 (對角係數 = 0)
    a[p6][p6-1] = dx*(-7.0/8.0) ;     //西側計算點
    a[p6][p6-2] = dx*(1.0/8.0) ;       //西西側計算點
    a[p6][p6+NX] = (3.0/8.0)*dx ;     //北側計算點
    a[p6][p6-NX] = -dx ;    //南側計算點
    b[p6] = dx*(-2.0/8.0)*T_bottom;
    
    ///////////////////////////////////////////////////////////////
    //上邊界角點處理（最後一排）- Neumann條件
    ///////////////////////////////////////////////////////////////
    
    //第七個角點 (1,NX)
    int p7 = (NX-1)*NX + 1;
    a[p7][p7] = a4 + dx*(1-D(r_s[1-1])) ;         //本點（原始係數組合）
    a[p7][p7+1] = a2 ;            //東側計算點           
    a[p7][p7-NX] = -dx*U(r_s[1-1]);            //南側計算點
    a[p7][p7-2*NX] = -dx*UU(r_s[1-1]) ;      //南南側計算點
    b[p7] = dx*(10.0/8.0)*T_left ;
    
    ///第八個角點 (2,NX)
    int p8 = (NX-1)*NX + 2;
    // 原本北側的a2項，因為T_N = T_{N-1}，所以移到對角線
    a[p8][p8] = a2 + dx*(1-D(r_s[2-1])) ;         //本點（a2原有 + a2從北側移來）
    a[p8][p8+1] = a2;            //東側計算點
    a[p8][p8-1] = -dx;           //西側計算點
    a[p8][p8-NX] = -dx*U(r_s[2-1]);          //南側計算點
    a[p8][p8-2*NX] = -dx*UU(r_s[2-1]) ;       //南南側計算點
    b[p8] = dx*(-2.0/8.0)*T_left  ;
    
    //第九個角點 (NX, NX)
    int p9 = (NX-1)*NX + NX;
    a[p9][p9] = dx*(3.0/8.0) + dx*(1-D(r_s[NX-1])) ;         //本點（a2原有 + a2從北側移來）
    a[p9][p9-1] = dx*(-7.0/8.0) ;           //西側計算點
    a[p9][p9-2] = dx*(1.0/8.0) ; 
    a[p9][p9-NX] = -dx*U(r_s[NX-1]);      //南側計算點
    a[p9][p9-2*NX] = -dx*UU(r_s[NX-1]) ;        //南南側計算點
    b[p9] = 0.0   ;
    
    /////////////////////////////////////////////////////////////////////
    // 邊界處理（不含角點）
    /////////////////////////////////////////////////////////////////////
    
    // 下邊界 (除角點外)(i = 3~NX-1 ; j = 1) 
    for(int i = 3; i <= NX-1; i++) {
        int index = (1-1)*NX + i;
        a[index][index] = a2+a4;   //本點
        a[index][index+1] = a2;     //東側計算點
        a[index][index-1] = -a4;    //西側計算點
        a[index][index-2] = a1;     //西西側計算點
        a[index][index+NX] = a2;    //北側計算點
        b[index] = dx*(10.0/8.0)*T_bottom;
    }
    
    //上邊界(i = 3~NX-1 ; j = NX) 
    for(int i = 3; i <= NX-1; i++) {
        int index = (NX-1)*NX + i;
        a[index][index] = a2 + dx*(1-D(r_s[i-1])) ;      //本點（注意：無上側項）
        a[index][index+1] = a2;     //東側計算點
        a[index][index-1] = -a4;    //西側計算點
        a[index][index-2] = a1;     //西西側計算點
        a[index][index-NX] = -dx*U(r_s[i-1]);
        a[index][index-2*NX] = -dx*UU(r_s[i-1]) ;
        b[index] = 0 ;         
    }
    //Break
    
    // 左邊界 (除角點外)
    for(int j = 3; j <= NX-1; j++) {
        int index = (j-1)*NX + 1;
        a[index][index] = a4+a2;    //本點
        a[index][index+1] = a2;      //東側計算點
        a[index][index-NX] = -a4;    //南側計算點
        a[index][index+NX] = a2;     //北側計算點
        a[index][index-2*NX] = a1;   //南南側計算點
        b[index] = dx*(10.0/8.0)*T_left;
    }
    
    // 右邊界 (除角點外)
    for(int j = 3; j <= NX-1; j++) {
        int index = (j-1)*NX + NX;
        a[index][index] = dx*(6.0/8.0) + dx*(3.0/8.0) ;   //本點
        a[index][index-1] = dx*(-7.8/8.0) ;     //西側計算點
        a[index][index-2] = dx*(1.0/8.0) ;      //西西側計算點
        a[index][index+NX] = dx*(3.0/8.0) ;     //北側計算點
        a[index][index-NX] = dx*(-7.0/8.0) ;    //南側計算點
        a[index][index-2*NX] = dx*(1.0/8.0) ;   //南南側計算點
        b[index] = 0.0; 
    }
    
    // 左邊第二層邊界
    for(int j = 3; j <= NX-1; j++) {
        int index = (j-1)*NX + 2;
        a[index][index] = a2+a2;     //本點
        a[index][index+1] = a2;       //東側計算點
        a[index][index-1] = -dx;      //西側計算點
        a[index][index+NX] = a2;      //北側計算點
        a[index][index-NX] = -a4;     //南側計算點
        a[index][index-2*NX] = a1;    //南南側計算點
        b[index] = dx*(-2.0/8.0)*T_left;
    }
    
    // 下面第二層邊界
    for(int i = 3; i <= NX-1; i++){
        int index = (2-1)*NX + i;
        a[index][index] = a2+a2;     //本點
        a[index][index+1] = a2;       //東側計算點
        a[index][index-1] = -a4 ;      //西側計算點
        a[index][index-2] = a1;       //西西側計算點
        a[index][index+NX] = a2;      //北側計算點
        a[index][index-NX] = -dx ;     //南側計算點
        b[index] = dx*(-2.0/8.0)*T_bottom;
    }
    
    // 內點
    for(int i = 3; i <= NX-1; i++) {
        for(int j = 3; j <= NX-1; j++) {
            int index = NX*(j-1)+i;
            a[index][index] = a2+a2;     //本點
            a[index][index-1] = -a4;      //西側計算點
            a[index][index+1] = a2;       //東側計算點
            a[index][index-2] = a1;       //西西側計算點
            a[index][index-NX] = -a4;     //南側計算點
            a[index][index+NX] = a2;      //北側計算點
            a[index][index-2*NX] = a1;    //南南側計算點
            b[index] = 0.0;
        }
    }
}
//Jacobi迭代求解
void Jacobi(vector<vector<double> >& a, vector<double>& b, vector<double>& x, int n) {
    for(int i = 1 ; i <= n ; i++){
    	x_old[i] = x[i] ;
	}
    for(int k = 1; k <= n; k++) {
        if(fabs(a[k][k]) < 1e-15) continue; // 跳過奇異矩陣
        
        double sum = 0;
        for(int p = 1; p <= n; p++) {
            if(p != k) {
                sum += a[k][p] * x[p];
            }
        }
        double x_new = (b[k] - sum) / a[k][k];
        x[k] = x_old[k] + lamda * (x_new - x_old[k]);
    }
    for(int i = 1 ; i <= NX ; i++){
    	int index = (NX-1)*NX + i ;
    	r_s[i-1] = (x[index-NX] - x[index-2*NX])/(x[index] - x[index-NX]) ;
	
	} 
    //計算當前步最大誤差
    maxerror = 0;
    for(int k = 1; k <= n; k++) {
        double error = fabs(x[k] - x_old[k]);
        if(maxerror < error) {
            maxerror = error;
        }
    }//計算i = 1 : NX ; j = 1 : NX-1 的最大誤差 
}

//輸出VTK檔案
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
    name << "QUICK北邊界BW南邊界MUSCL,Pe=infinite,Grid "<< nodes_i << "x" << nodes_j << ","<< setfill('0') << setw(6) << m << ".vtk";
    ofstream out(name.str().c_str());
    
    // VTK 文件頭：使用 STRUCTURED_GRID
    out << "# vtk DataFile Version 3.0\n";
    out << "QUICK北邊界BW南邊界MUSCL,Pe=infinite,Grid "<< nodes_i << "x" << nodes_j << ","<< setfill('0') << setw(6) << m << ".vtk\n";
    out << "ASCII\n";
    out << "DATASET STRUCTURED_GRID\n";  // 改為 STRUCTURED_GRID
    out << "DIMENSIONS " << nodes_i << " " << nodes_j << " 1\n";
    
    // 輸出節點座標
    out << "POINTS " << nodes_i * nodes_j << " double\n";
    for(int j = 0; j < nodes_j; j++) {
        for(int i = 0; i < nodes_i; i++) {
            double x_coord = i * dx;
            double y_coord = j * dx;
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
        // QUICK格式係數定義
        a1 = dx*(1.0/8.0);  //係數1/8
        a2 = dx*(3.0/8.0);  //係數3/8
        a3 = dx*(6.0/8.0);  //係數6/8
        a4 = dx*(7.0/8.0);  //係數7/8
        a5 = dx*(5.0/8.0); //係數12/8 
        n = (NX)*(NX);
        // 重新調整向量大小
        a.assign(n+2, vector<double>(n+2, 0.0));
        b.assign(n+2, 0.0);
        x.assign(n+2, 0.0);
        x_old.assign(n+2, 0.0); //向量大小必須唯一 
        T.assign(NX, vector<double>(NX, 0.0));
        r_s.assign(NX,0.0) ;
        cout << "========================================" << endl;
        cout << "程式碼開始執行...." << endl;
        cout << "網格大小: " << NX+1 << " x " << NX+1 << endl;
        cout << "網格間距: dx = dy = " << dx << endl;
        cout << "邊界條件設定：" << endl;
        cout << "  左邊界(Dirichlet): T = " << T_left << endl;
        cout << "  右邊界(Dirichlet): T = " << T_right << endl;
        cout << "  下邊界(Dirichlet): T = " << T_bottom << endl;
        cout << "  上邊界(Neumann):  T_[i][NY] = T_[i][NY+1](虛擬點)" << endl;
        cout << "鬆弛因子: λ = " << lamda << endl;
        cout << "收斂準則: " << tolerance << endl;
        cout << "========================================" << endl;
        steadystate = false;
        Gamma0(a, b ,x);         //設定係數矩陣（含邊界條件）
        for(G = 0; G < max_G; G++) {
            Jacobi(a, b, x, n);  //執行一次Jacobi迭代
            
            if(G % 100 == 0) {
                cout << "迭代次數 = " << G;
                if(G > 0) {
                    cout << ", 最大變化量 = " << scientific << maxerror;
                }
                cout << endl;
                output(G);
            }
            
            if(G > 100 && maxerror < tolerance) {
                steadystate = true;
                cout << "\n已達收斂條件！" << endl;
                cout << "最終迭代次數: " << G << endl;
                cout << "最終誤差: " << scientific << maxerror << endl;
                break;
            }
        }
        
        if(!steadystate) {
            cout << "\n警告：達到最大迭代次數，但未達到穩態！" << endl;
        }
        
        output(G);  //輸出最終結果
        cout << "網格 " << NX+1 << "x" << NX+1 << " 計算完成\n" << endl;
    }
    cout << "所有計算完成！" << endl;
    return 0;
}