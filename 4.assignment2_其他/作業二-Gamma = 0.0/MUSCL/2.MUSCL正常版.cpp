//利用TVD-MUSCL格式求解二維穩態純對流方程式
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
///////////////////////////////////////////////////////////////  
vector<vector<double> >rf_w  ;
vector<vector<double> >coe_w  ;
vector<vector<double> >rf_e  ;
vector<vector<double> >coe_e  ;
vector<vector<double> >rf_s  ;
vector<vector<double> >coe_s  ;
vector<vector<double> >rf_n  ;
vector<vector<double> >coe_n ;
double dx; //一個網格間距
//邊界條件定義
const double T_left = 1.0;    //左邊界Dirichlet條件
const double T_right = 0.0;   //右邊界Dirichlet條件  
const double T_bottom = 0.0;  //下邊界Dirichlet條件
//上邊界為Neumann條件: ?T/?n = 0
//對上邊界的Nuemann採用一階精度後項差分,將邊界上的溫度賦值為中心計算點溫度 
//QUICK格式係數
double psie , psiw , psin , psis ;
//迭代參數
const double lamda = 0.3; //逐次超鬆弛迭代法
int G, max_G = 10000000;
double maxerror ; 
const float tolerance = 1e-10;
bool steadystate;

double MUSCL(double x){
   
    return std::min(0.0, std::max(-2.0*x, std::max(-(x+1.0)/2.0, -2.0)));
}

// 更新T矩陣從x向量
void updateT() {
    for(int j = 1; j <= NX; j++) {
        for(int i = 1; i <= NX; i++) {
            T[i-1][j-1] = x[(j-1)*NX + i];
        }
    }
}

//係數矩陣賦值（考慮Neumann邊界條件）
void Gamma0(vector<vector<double> >& a, vector<double>& b , vector<double>& x) {
    // 更新T矩陣 - 這是關鍵修正！
    updateT();
    
    // 初始化矩陣
    for(int i = 1; i <= NX*NX ; i++) {
        for(int j = 1; j <= NX*NX; j++) {
            a[i][j] = 0.0 ;
        }  
        b[i] = 0.0 ;
    }
    
    // 初始化所有TVD係數為0
    for(int i = 0; i < NX; i++) {
        for(int j = 0; j < NX; j++) {
            rf_w[i][j] = 0.0;
            coe_w[i][j] = 0.0;
            rf_e[i][j] = 0.0;
            coe_e[i][j] = 0.0;
            rf_s[i][j] = 0.0;
            coe_s[i][j] = 0.0;
            rf_n[i][j] = 0.0;
            coe_n[i][j] = 0.0;
        }
    }
    
    //該次迭代,該點位置，係數更新
    //西邊界係數
    for(int j = 1 ; j <= NX ; j++){
        for(int i = 3 ; i <= NX ; i++){
            double T_D1 = T[i-1][j-1] ;//本點 
            double T_C1  = T[i-2][j-1] ;//西點 
            double T_U1 = T[i-3][j-1] ;//西西點
            if(fabs(T_D1-T_C1) > 1e-10) { 
                rf_w[i-1][j-1]= (T_C1-T_U1)/(T_D1-T_C1) ;
                coe_w[i-1][j-1] = MUSCL(rf_w[i-1][j-1])/2.0 ;
            }
        }
        //單獨處理i = 2 的情況
        if(NX >= 2) {
            double T_D1 = T[1][j-1] ;//本點 
            double T_C1  = T[0][j-1] ;//西點 
            double T_U1prime = 2*T_left - T[0][j-1] ;
            if(fabs(T_D1-T_C1) > 1e-10) {
                rf_w[1][j-1]= (T_C1-T_U1prime)/(T_D1-T_C1) ;
                coe_w[1][j-1] = MUSCL(rf_w[1][j-1])/2.0 ;
            }
        }
    }
    
    //南邊界係數
    for(int i = 1 ; i <= NX ; i++){
        for(int j = 3 ; j <= NX ; j++){
            double T_D2 = T[i-1][j-1] ;//本點 
            double T_C2  = T[i-1][j-2] ;//南點 
            double T_U2 = T[i-1][j-3] ;//南南點
            if(fabs(T_D2-T_C2) > 1e-10) { 
                rf_s[i-1][j-1]= (T_C2-T_U2)/(T_D2-T_C2) ;
                coe_s[i-1][j-1] = MUSCL(rf_s[i-1][j-1])/2.0 ;
            }
        }
        //單獨處理j = 2的情況		
        if(NX >= 2) {
            double T_D2 = T[i-1][1] ;//本點 
            double T_C2  = T[i-1][0] ;//南點 
            double T_U2prime = 2*T_bottom - T[i-1][0] ;
            if(fabs(T_D2-T_C2) > 1e-10) {
                rf_s[i-1][1]= (T_C2-T_U2prime)/(T_D2-T_C2) ;
                coe_s[i-1][1] = MUSCL(rf_s[i-1][1])/2.0 ;
            }
        }
    }
    
    //東邊界係數
    for(int j = 1 ; j <= NX ; j++){
        for(int i = 2 ; i <= NX-1 ; i++){
            double T_D3 = T[i][j-1] ;//東點 
            double T_C3  = T[i-1][j-1] ;//本點 
            double T_U3 = T[i-2][j-1] ;//西點
            if(fabs(T_D3-T_C3) > 1e-10) { 
                rf_e[i-1][j-1]= (T_C3-T_U3)/(T_D3-T_C3) ;
                coe_e[i-1][j-1] = MUSCL(rf_e[i-1][j-1])/2.0 ;
            }
        }
        //單獨處理i = 1的情況 
        double T_D3 = T[1][j-1] ;//東點 
        double T_C3  = T[0][j-1] ;//本點 
        double T_U3prime = 2*T_left - T[0][j-1] ;
        if(fabs(T_D3-T_C3) > 1e-10) {
            rf_e[0][j-1]= (T_C3-T_U3prime)/(T_D3-T_C3) ;
            coe_e[0][j-1] = MUSCL(rf_e[0][j-1])/2.0 ;
        }
    }
    
    //北邊界係數
    for(int i = 1 ; i <= NX ; i++){
        for(int j = 2 ; j <= NX-1 ; j++){
            double T_D4 = T[i-1][j] ;//北點 
            double T_C4  = T[i-1][j-1] ;//本點 
            double T_U4 = T[i-1][j-2] ;//南點
            if(fabs(T_D4-T_C4) > 1e-10) { 
                rf_n[i-1][j-1] = (T_C4-T_U4)/(T_D4-T_C4) ;
                coe_n[i-1][j-1] = MUSCL(rf_n[i-1][j-1])/2.0 ;
            }
        }
        //單獨處理j = 1的情況
        double T_D4 = T[i-1][1] ;//北點 
        double T_C4  = T[i-1][0] ;//本點 
        double T_U4prime = 2*T_bottom - T[i-1][0] ;
        if(fabs(T_D4-T_C4) > 1e-10) {
            rf_n[i-1][0]= (T_C4-T_U4prime)/(T_D4-T_C4) ;
            coe_n[i-1][0] = MUSCL(rf_n[i-1][0])/2.0 ;
        }
    }
    
    // 組裝矩陣
    for(int i = 1 ; i <= NX ; i++){
        for(int j = 1 ; j <= NX ; j++){
            int index = (j-1)*NX + i ; 
            
            // 角點處理
            if(i == 1 && j == 1){  // 左下角
                a[index][index] = dx*(1-coe_e[0][0]) + dx*(1-coe_n[0][0]);
                a[index][index+1] = dx*coe_e[0][0];
                a[index][index+NX] = dx*coe_n[0][0];
                b[index] = dx*T_left + dx*T_bottom;
            }
            else if(i == NX && j == 1){  // 右下角
                a[index][index] = dx*(-coe_w[NX-1][0]) + dx*(1-coe_n[NX-1][0]);
                a[index][index-1] = -dx*(1-coe_w[NX-1][0]);
                a[index][index+NX] = dx*coe_n[NX-1][0];
                b[index] = -dx*T_right + dx*T_bottom;	
            }
            else if(i == 1 && j == NX){  // 左上角（Neumann邊界）
                a[index][index] = dx*(1-coe_e[0][NX-1]) + dx*(1-coe_s[0][NX-1]);
                a[index][index+1] = dx*coe_e[0][NX-1];
                a[index][index-NX] = -dx*(1-coe_s[0][NX-1]);
                b[index] = dx*T_left; 
            }
            else if(i == NX && j == NX){  // 右上角（Neumann邊界）
                a[index][index] = dx*(-coe_w[NX-1][NX-1]) + dx*(1-coe_s[NX-1][NX-1]);
                a[index][index-1] = -dx*(1-coe_w[NX-1][NX-1]);
                a[index][index-NX] = -dx*(1-coe_s[NX-1][NX-1]);
                b[index] = -dx*T_right;
            }
            // 下邊界 (除角點外)
            else if(j == 1 && i > 1 && i < NX) {
                a[index][index] = dx*(1-coe_e[i-1][0]-coe_w[i-1][0])+dx*(1-coe_n[i-1][0]);
                a[index][index+1] = dx*coe_e[i-1][0];
                a[index][index-1] = -dx*(1-coe_w[i-1][0]);
                a[index][index+NX] = dx*coe_n[i-1][0];
                b[index] = dx*T_bottom;
            }
            // 上邊界 (除角點外)
            else if(j == NX && i > 1 && i < NX) {
                a[index][index] = dx*(1-coe_e[i-1][NX-1]-coe_w[i-1][NX-1])+dx*(1-coe_s[i-1][NX-1]);
                a[index][index+1] = dx*coe_e[i-1][NX-1];
                a[index][index-1] = -dx*(1-coe_w[i-1][NX-1]);
                a[index][index-NX] = -dx*(1-coe_s[i-1][NX-1]);
                b[index] = 0.0;
            }
            // 左邊界 (除角點外)
            else if(i == 1 && j > 1 && j < NX) {
                a[index][index] = dx*(1-coe_e[0][j-1])+dx*(1-coe_n[0][j-1]-coe_s[0][j-1]);
                a[index][index+1] = dx*coe_e[0][j-1];
                a[index][index+NX] = dx*coe_n[0][j-1];
                a[index][index-NX] = -dx*(1-coe_s[0][j-1]);
                b[index] = dx*T_left;
            }
            // 右邊界 (除角點外)
            else if(i == NX && j > 1 && j < NX) {
                a[index][index] = dx*(-coe_w[NX-1][j-1])+dx*(1-coe_n[NX-1][j-1]-coe_s[NX-1][j-1]);
                a[index][index-1] = -dx*(1-coe_w[NX-1][j-1]);
                a[index][index+NX] = dx*coe_n[NX-1][j-1];
                a[index][index-NX] = -dx*(1-coe_s[NX-1][j-1]);
                b[index] = -dx*T_right;
            }
            // 內點
            else {
                a[index][index] = dx*(1-coe_e[i-1][j-1]-coe_w[i-1][j-1])+dx*(1-coe_n[i-1][j-1]-coe_s[i-1][j-1]);
                a[index][index+1] = dx*coe_e[i-1][j-1];
                a[index][index-1] = -dx*(1-coe_w[i-1][j-1]);
                a[index][index+NX] = dx*coe_n[i-1][j-1];
                a[index][index-NX] = -dx*(1-coe_s[i-1][j-1]);
                b[index] = 0.0;
            }
        }
    }
}

void SOU(vector<vector<double> >& a, vector<double>& b, vector<double>& x, int n) {
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
    
    //計算當前步最大誤差
    maxerror = 0;
    for(int k = 1; k <= n; k++) {
        double error = fabs(x[k] - x_old[k]);
        if(maxerror < error) {
            maxerror = error;
        }
    }
}

// 初始化溫度場
void initializeField() {
    // 設定初始猜測值
    for(int j = 1; j <= NX; j++) {
        for(int i = 1; i <= NX; i++) {
            int index = (j-1)*NX + i;
            // 線性插值作為初始值
            double x_pos = (i-1)*dx;
            double y_pos = (j-1)*dx;
            x[index] = T_left * (1.0 - x_pos);
        }
    }
    updateT();
}

//輸出VTK檔案
void output(int m) {
    // 更新T矩陣
    updateT();
    
    ostringstream name;
    name << "MUSCL"<< NX+1 << "x" << NX+1 << "_step" << setfill('0') << setw(6)  <<  m << ".vtk";
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
        n = NX*NX;
        
        // 重新調整向量大小
        a.assign(n+2, vector<double>(n+2, 0.0));
        b.assign(n+2, 0.0);
        x.assign(n+2, 0.0);
        x_old.assign(n+2, 0.0);
        T.assign(NX, vector<double>(NX, 0.0));
        
        rf_w.assign(NX, vector<double>(NX, 0.0));
        coe_w.assign(NX, vector<double>(NX, 0.0));
        rf_e.assign(NX, vector<double>(NX, 0.0));
        coe_e.assign(NX, vector<double>(NX, 0.0));
        rf_s.assign(NX, vector<double>(NX, 0.0));
        coe_s.assign(NX, vector<double>(NX, 0.0));
        rf_n.assign(NX, vector<double>(NX, 0.0));
        coe_n.assign(NX, vector<double>(NX, 0.0));
        
        cout << "========================================" << endl;
        cout << "程式碼開始執行...." << endl;
        cout << "網格大小: " << NX+1 << " x " << NX+1 << endl;
        cout << "網格間距: dx = dy = " << dx << endl;
        cout << "邊界條件設定：" << endl;
        cout << "  左邊界(Dirichlet): T = " << T_left << endl;
        cout << "  右邊界(Dirichlet): T = " << T_right << endl;
        cout << "  下邊界(Dirichlet): T = " << T_bottom << endl;
        cout << "  上邊界(Neumann):  ?T/?n = 0" << endl;
        cout << "鬆弛因子: λ = " << lamda << endl;
        cout << "收斂準則: " << tolerance << endl;
        cout << "========================================" << endl;
        
        // 初始化溫度場
        initializeField();
        
        steadystate = false;
        for(G = 0; G < max_G; G++) {
            Gamma0(a, b, x);         //設定係數矩陣（含邊界條件）
            SOU(a, b, x, n);        //執行一次SOU迭代
            
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
