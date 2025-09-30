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
int Nx[] = {40,80};
//有限體積法採用Wet node處理，故計算點數目 = 網格數目, 節點數目 = 網格數目+1
int NX; //NX = NY => dx = dy
int n;

vector<double> x;
vector<double> x_old;
vector<vector<double> > T;
double dx; //一個網格間距
//邊界條件定義
const double T_left = 1.0;    //左邊界Dirichlet條件
const double T_right = 0.0;   //右邊界Dirichlet條件  
const double T_bottom = 0.0;  //下邊界Dirichlet條件
//迭代參數
const double lamda = 0.1; //逐次超鬆弛迭代法

int G, max_G = 10000000;
const double w = 1e-6 ;
double maxerror ;
const float tolerance = 1e-10;
bool steadystate;
double MUSCL(double r){
if (r <= 0.0) return 0.0001;
    return std::min(std::min(1.0, r), 0.55);  // 額外限制到 0.55
}
//TVD-MUSCL格式為在東西方向為三點格式,南北方向為三點格式

//顯示迭代法（考慮Neumann邊界條件）
void Gamma0(vector<double>& x) { //因為係數式時變的,每一次迭代都必須要做.....

    int index1 = (1-1)*NX+1 ;
    x_old[index1] = x[index1] ;
    double x_new1 =  ((dx*T_left + dx*T_bottom) - x[index1+1]*dx*(MUSCL(2*(T[1][1]-T_left)/(T[2][1]-T[1][1])))/(2.0) - x[index1+NX]*dx*(MUSCL(2.0*(T[1][1]-T_bottom)/(T[1][2]-T[1][1]))/2.0))/ (dx*(1-(MUSCL(2*(T[1][1]-T_left)/(T[2][1]-T[1][1])))/(2.0)) + dx*(1-(MUSCL(2.0*(T[1][1]-T_bottom)/(T[1][2]-T[1][1]))/2.0))) ;
    x[index1] = x_old[index1] + lamda * (x_new1 - x_old[index1]);
   
    int index2 = (1-1)*NX+NX ;      
    x_old[index2] = x[index2] ;
    double x_new2 =  ((-dx*T_right + dx*T_bottom) + x[index2-1]*dx*(1-(MUSCL((T[NX-1][1]-T[NX-2][1])/(T[NX][1]-T[NX-1][1]))/2.0)) - x[index2+NX]*dx*(MUSCL(2.0*(T[1][1]-T_bottom)/(T[1][2]-T[1][1]))/2.0) )/((-dx*(MUSCL(T[NX-1][1]-T[NX-2][1])/(T[NX][1]-T[NX-1][1]))/2.0) + dx*(1-(MUSCL(2.0*(T[1][1]-T_bottom)/(T[1][2]-T[1][1]))/2.0)));
    x[index2] = x_old[index2] + lamda * (x_new2 - x_old[index2]);

    int index3 = (NX-1)*NX+1 ;
    x_old[index3] = x[index3] ;
    double x_new3 =  (dx*T_left - x[index3+1]*dx*(MUSCL(2*(T[1][NX]-T_left)/(T[2][NX]-T[1][NX]))/2.0) + x[index3-NX]*(dx*MUSCL((T[1][NX-1]-T[1][NX-2])/(T[1][NX]-T[1][NX-1]))/2.0) )/(dx*(1-(MUSCL(2*(T[1][NX]-T_left)/(T[2][NX]-T[1][NX]))/2.0)) + dx*(1-MUSCL((T[1][NX-1]-T[1][NX-2])/(T[1][NX]-T[1][NX-1]))/2.0)) ;
    x[index3] = x_old[index3] + lamda * (x_new3 - x_old[index3]);
           
    int index4 = (NX-1)*NX+NX ;        
    x_old[index4] = x[index4] ;
    double x_new4 =  ((-dx*T_right) + x[index4-1]*(dx*(1-MUSCL((T[NX-1][NX]-T[NX-2][NX])/(T[NX][NX]-T[NX-1][NX]+w))/2.0)) + x[index4-NX]*(dx*(1-MUSCL((T[NX][NX-1]-T[NX][NX-2])/(T[NX][NX]-T[NX][NX-1]+w))/2.0)))/(dx*(-MUSCL((T[NX-1][NX]-T[NX-2][NX])/(T[NX][NX]-T[NX-1][NX]+w))/2.0) + dx*(1-MUSCL((T[NX][NX-1]-T[NX][NX-2])/(T[NX][NX]-T[NX][NX-1]+w))/2.0)) ;
    x[index4] = x_old[index4] + lamda * (x_new4 - x_old[index4]);
           
    // 邊界處理（不含角點）
    // 下邊界 (除角點外)(i = 2~NX-1 ; j = 1)
    for(int i = 2 ; i <= NX-1; i++) {
        int index = (1-1)*NX + i;
        x_old[index] = x[index] ;
        double x_new = (dx*T_bottom - x[index+1]*dx*MUSCL((T[i][1]-T[i-1][1])/(T[i+1][1]-T[i][1]+w))/2.0 + x[index-1]*(dx*(1-MUSCL((T[i-1][1]-T[i-2][1])/(T[i][1]-T[i-1][1]+w))/2.0)) -x[index+NX]*dx*MUSCL((T[i][1]-T[i][0])/(T[i][2]-T[i][1]+w))/2.0 )/(dx*(1-MUSCL((T[i][1]-T[i-1][1])/(T[i+1][1]-T[i][1]+w))/2.0-MUSCL((T[i-1][1]-T[i-2][1])/(T[i][1]-T[i-1][1]+w))/2.0)+dx*(1-MUSCL((T[i][1]-T[i][0])/(T[i][2]-T[i][1]+w))/2.0));
        x[index] = x_old[index] + lamda * (x_new - x_old[index]);
    }
    //上邊界(i = 2~NX-1 ; j = NX)
    for(int i = 2 ; i <= NX-1; i++) {
    int index = (NX-1)*NX + i;
        x_old[index] = x[index] ;
        double x_new = (0.0 - x[index+1]*dx*MUSCL((T[i][NX]-T[i-1][NX])/(T[i+1][NX]-T[i][NX]+w))/2.0 + x[index-1]*(dx*(1-MUSCL((T[i-1][NX]-T[i-2][NX])/(T[i][NX]-T[i-1][NX]+w))/2.0)) + x[index-NX]*(dx*(1-MUSCL((T[i][NX-1]-T[i][NX-2])/(T[i][NX]-T[i][NX-1]+w))/2.0)))/(dx*(1-MUSCL((T[i][NX]-T[i-1][NX])/(T[i+1][NX]-T[i][NX]+w))/2.0-MUSCL((T[i-1][NX]-T[i-2][NX])/(T[i][NX]-T[i-1][NX]+w))/2.0)+dx*(1-MUSCL((T[i][NX-1]-T[i][NX-2])/(T[i][NX]-T[i][NX-1]+w))/2.0)) ;
        x[index] = x_old[index] + lamda * (x_new - x_old[index]);
    }

    // 左邊界 (除角點外)
    for(int j = 2; j <= NX-1; j++) {
        int index = (j-1)*NX + 1;
        x_old[index] = x[index] ;
        double x_new = (dx*T_left - x[index+1]*dx*MUSCL((T[1][j]-T[0][j])/(T[2][j]-T[1][j]+w))/2.0 - x[index+NX]*dx*MUSCL((T[1][j]-T[1][j-1])/(T[1][j+1]-T[1][j]+w))/2.0 +x[index-NX]*(dx*(1-MUSCL((T[1][j-1]-T[1][j-2])/(T[1][j]-T[1][j-1]+w))/2.0)) )/(dx*(1-MUSCL((T[1][j]-T[0][j])/(T[2][j]-T[1][j]+w))/2.0)+dx*(1-MUSCL((T[1][j]-T[1][j-1])/(T[1][j+1]-T[1][j]+w))/2.0-MUSCL((T[1][j-1]-T[1][j-2])/(T[1][j]-T[1][j-1]+w))/2.0)) ;
        x[index] = x_old[index] + lamda * (x_new - x_old[index]);
    }
   
    // 右邊界 (除角點外)
    for(int j = 2; j <= NX-1; j++) {
        int index = (j-1)*NX + NX;
        x_old[index] = x[index] ;
        double x_new = (-dx*T_right +x[index-1]*(dx*(1-MUSCL((T[NX-1][j]-T[NX-2][j])/(T[NX][j]-T[NX-1][j]+w))/2.0)) -x[index+NX]*dx*MUSCL((T[NX][j]-T[NX][j-1])/(T[NX][j+1]-T[NX][j]+w))/2.0 +x[index-NX]*(dx*(1-MUSCL((T[NX][j-1]-T[NX][j-2])/(T[NX][j]-T[NX][j-1]+w))/2.0)) )/(dx*(-MUSCL((T[NX-1][j]-T[NX-2][j])/(T[NX][j]-T[NX-1][j]+w))/2.0)+dx*(1-MUSCL((T[NX][j]-T[NX][j-1])/(T[NX][j+1]-T[NX][j]+w))/2.0-MUSCL((T[NX][j-1]-T[NX][j-2])/(T[NX][j]-T[NX][j-1]+w))/2.0)) ;
        x[index] = x_old[index] + lamda * (x_new - x_old[index]);
    }
    // 內點
    for(int i = 2; i <= NX-1; i++) {
        for(int j = 2; j <= NX-1; j++) {
            int index = (j-1)*NX+i;
            x_old[index] = x[index] ;
            double x_new = (0.0 - x[index+1]*dx*MUSCL((T[i][j]-T[i-1][j])/(T[i+1][j]-T[i][j]+w))/2.0 + x[index-1]*(dx*(1-MUSCL((T[i-1][j]-T[i-2][j-1])/(T[i][j]-T[i-1][j]+w))/2.0)) - x[index+NX]*dx*MUSCL((T[i][j]-T[i][j-1])/(T[i][j+1]-T[i][j]+w))/2.0 + x[index-NX]*(dx*(1-MUSCL((T[i][j-1]-T[i][j-2])/(T[i][j]-T[i-1][j]+w))/2.0)))/(dx*(1-MUSCL((T[i][j]-T[i-1][j])/(T[i+1][j]-T[i][j]+w))/2.0-MUSCL((T[i-1][j]-T[i-2][j])/(T[i][j]-T[i-1][j]+w))/2.0)+dx*(1-MUSCL((T[i][j]-T[i][j-1])/(T[i][j+1]-T[i][j]+w))/2.0-MUSCL((T[i][j-1]-T[i][j-2])/(T[i][j]-T[i-1][j]+w))/2.0)) ;
            x[index] = x_old[index] + lamda * (x_new - x_old[index]);
        }
    }
}

void SOU(vector<double>& x, int n) {
    // 將一維解轉換為二維溫度場
    for(int j = 1; j <= NX; j++) {
        for(int i = 1; i <= NX; i++) {
            T[i][j] = x[(j-1)*NX + i];
        }
    }
    for(int j  = 1 ; j <= NX ; j++){
        T[0][j] = 2*T_left - T[1][j] ; //左邊界虛擬點
        T[NX+1][j] = 0.0 ; 
    }
    for(int i = 1;  i<= NX ; i++){
        T[i][0] = 2*T_bottom - T[i][1] ; //下邊界虛擬點
        T[i][NX+1] = 0.0 ; 
    }
    T[0][0] = 0.0 ;
    T[NX+1][0] = 0.0 ;
    T[0][NX+1] = 0.0 ;
    T[NX+1][NX+1] = 0.0 ; //溫度更新完畢
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

//輸出VTK檔案
void output(int m) {
    ostringstream name;
    name << "MUSCL"<< NX+1 << "x" << NX+1 << "步數 = "<< setfill('0') << setw(6)  <<  m << ".vtk";
    ofstream out(name.str().c_str());
   
    // VTK 文件頭
    out << "# vtk DataFile Version 3.0\n";
    out << "MUSCL\n";
    out << "ASCII\n";
    out << "DATASET STRUCTURED_POINTS\n";
    out << "DIMENSIONS " << NX << " " << NX << " 1\n";
    out << "ORIGIN 0 0 0\n";
    out << "SPACING " << dx << " " << dx << " 1\n";
    out << "POINT_DATA " << NX*NX << "\n";
   
    // 輸出溫度場
    out << "SCALARS Temperature double 1\n";
    out << "LOOKUP_TABLE default\n";
    for(int j = 0; j < NX; j++) {
        for(int i = 0; i < NX; i++) {
            out << scientific << setprecision(6) << T[i][j] << "\n";
        }
    }
    out.close();
    cout << "VTK 文件已輸出: " << name.str() << endl;
}

int main() {
    for(int grid_idx = 0; grid_idx < 2; grid_idx++){
        NX = Nx[grid_idx];
        dx = 1.0/double(NX);    
        n = (NX)*(NX);
        // 重新調整向量大小
        x.assign(n+2, 0.0);
        x_old.assign(n+2, 0.0); //向量大小必須唯一
        T.assign(NX+2, vector<double>(NX+2, 0.0));
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
        for(int i = 1 ; i <= n ; i++){
        x[i] = 0.0 ;
        }
       x[0] = 0.0 ;
       x[n+1] = 0.0 ;
        for(G = 0; G < max_G; G++) {
            Gamma0(x);  //設定係數矩陣（含邊界條件）
            SOU(x, n);  //計算最大誤差
            if(G % 10000 == 0) {
                cout << "迭代次數 = " << G;
                if(G > 0) {
                    cout << ", 最大變化量 = " << scientific << maxerror;
                }
                cout << endl;
                output(G);
            }
           
            if(G > 1000 && maxerror < tolerance) {
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