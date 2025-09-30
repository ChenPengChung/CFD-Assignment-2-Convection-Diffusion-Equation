//利用空間二階精度MINMOD格式求解二維穩態擴散對流方程
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
int Nx[] = {80,80};
//有限體積法採用Wet node處理，故計算點數目 = 網格數目, 節點數目 = 網格數目+1
int NX; //NX = NY => dx = dy
int n;
vector<vector<double> > a; //a二維矩陣為等式左側係數矩陣
vector<double> b; //b一維矩陣為等式右側之外源項，方程組之非齊性項
vector<double> x;
vector<double> x_old;
vector<vector<double> > T;
vector<vector<double> > r_w;
vector<vector<double> > r_e;
vector<vector<double> > r_s;
vector<vector<double> > r_n;
double dx; //一個網格間距
//邊界條件定義
const double T_left = 1.0;    //左邊界Dirichlet條件
const double T_right = 0.0;   //右邊界Dirichlet條件  
const double T_bottom = 0.0;  //下邊界Dirichlet條件
const double Pelect = 1000.0 ;
const double Gamma = sqrt(2.0)/Pelect ;//熱擴散係數
//迭代參數
const double lamda = 0.0001; //逐次超鬆弛迭代法
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
void initial(vector<vector<double> >& a, vector<double>& b , vector<double>& x){
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
    //初始化溫度場
    //左邊界計算點之左延伸
    for(int j = 1 ; j <= NX ; j++){
        T[0][j] = 2*T_left - T[1][j] ;
    }
    for(int j = 1 ; j <= NX ; j++){
        T[NX+1][j] = 2*T_right - T[NX][j] ;
    }
    for(int i = 1 ; i <= NX ; i++){
        T[i][0] = 2*T_bottom - T[i][1] ;
    }
    for(int i = 1 ; i <= NX ; i++){
        T[i][NX+1] = T[i][NX] ;
    }
}

//係數矩陣賦值（考慮Neumann邊界條件）
void Gamma0(vector<vector<double> >& a, vector<double>& b , vector<double>& x) { //每一次迭代要做的事情就是將a,b矩陣重新賦值
    for(int i = 1; i <= NX*NX ; i++) {
        for(int j = 1; j <= NX*NX; j++) {
            a[i][j] = 0.0 ;
        }  
        b[i] = 0.0 ;
    }//記得不要每一次迭代都把 解賦值為0
    a[0][0] = 0.0 ;
    a[NX*NX+1][NX*NX+1] = 1.0 ;
    ///////////////////////////////////////////////////////////////  
    //第一角點 (1,1) - 左下角
    int p1 = (1-1)*NX + 1;
    a[p1][p1] =(2*Gamma+dx*(U(r_n[1-1][1-1])))+(2*Gamma + dx*U(r_e[1-1][1-1]));//本點
    //a[p1][p1-1] =- ;//西側計算點
    a[p1][p1+1] =-(Gamma - dx*D(r_e[1-1][1-1])) ;//東側計算點
    //a[p1][p1-NX] =-;//南側計算點
    a[p1][p1+NX] =-(Gamma-dx*D(r_n[1-1][1-1]));//北側計算點
    b[p1] = dx*T_bottom + dx*T_left + (Gamma - dx*UU(r_e[1-1][1-1])) * T[0][1] + (Gamma-dx*UU(r_n[1-1][1-1])) * T[1][0];
    //第二角點 (NX,1) - 右下角
    int p2 = (1-1)*NX + NX;
    a[p2][p2] =(2*Gamma+dx*(U(r_n[NX-1][1-1])))+(2*Gamma+dx*(-D(r_w[NX-1][1-1]))) ;//本點
    a[p2][p2-1] =-(Gamma-dx*(-U(r_w[NX-1][1-1]))) ;//西側計算點
    //a[p2][p2+1] =- ;//東側計算點
    a[p2][p2-2] =-(-dx*(-UU(r_w[NX-1][1-1]))) ;//西西側計算點
    //a[p2][p2-NX] =-;//南側計算點
    a[p2][p2+NX] =-(Gamma-dx*D(r_n[NX-1][1-1]));//北側計算點
    b[p2] = dx*T_bottom-dx*T_right + (Gamma) * T[NX+1][1] + (Gamma-dx*UU(r_n[NX-1][1-1])) * T[NX][0];
    //第三個角點 (1,NX) - 左上角
    int p3 = (NX-1)*NX + 1;
    a[p3][p3] =(1*Gamma+dx*(U(r_n[1-1][NX-1])-D(r_s[1-1][NX-1]))) + (2*Gamma + dx*U(r_e[1-1][NX-1]));//本點
    //a[p3][p3-1] =- ;//西側計算點
    a[p3][p3+1] =-(Gamma - dx*D(r_e[1-1][NX-1])) ;//東側計算點
    a[p3][p3-NX] =-(Gamma-dx*(UU(r_n[1-1][NX-1])-U(r_s[1-1][NX-1])));//南側計算點
    //a[p3][p3+NX] =-;//北側計算點
    a[p3][p3-2*NX] =-(-dx*(-UU(r_s[1-1][NX-1])));//南南側計算點
    b[p3] = 0.0 + dx*T_left + (Gamma - dx*UU(r_e[1-1][NX-1])) * T[0][NX] + (-dx*D(r_n[1-1][NX-1]))*T[1][NX+1];
    //第四個角點 (NX, NX) - 右上角
    int p4 = (NX-1)*NX + NX;
    a[p4][p4] =(1*Gamma+dx*(U(r_n[NX-1][NX-1])-D(r_s[NX-1][NX-1])))+ (2*Gamma+dx*(-D(r_w[NX-1][NX-1]))); //本點
    a[p4][p4-1] =-(Gamma-dx*(-U(r_w[NX-1][NX-1]))) ;//西側計算點
    //a[p4][p4+1] =- ;//東側計算點
    a[p4][p4-2] =-(-dx*(-UU(r_w[NX-1][NX-1]))) ;//西西側計算點
    a[p4][p4-NX] =-(Gamma-dx*(UU(r_n[NX-1][NX-1])-U(r_s[NX-1][NX-1])));//南側計算點
    //a[p4][p4+NX] =-;//北側計算點
    a[p4][p4-2*NX] =-(-dx*(-UU(r_s[NX-1][NX-1])));//南南側計算點
    b[p4] = 0.0 -dx*T_right + (Gamma) * T[NX+1][NX] + (-dx*D(r_n[NX-1][NX-1])) * T[NX][NX+1];

    /////////////////////////////////////////////////////////////////////
    // 邊界處理（不含角點）
    /////////////////////////////////////////////////////////////////////
   
    // 下邊界 (除角點外)(i = 3~NX-1 ; j = 1)
    for(int i = 2; i <= NX-1; i++) {
        int index = (1-1)*NX + i;
        a[index][index] =(4*Gamma+dx*(U(r_e[i-1][1-1])-D(r_w[i-1][1-1]))+dx*(U(r_n[i-1][1-1])));//本點
        a[index][index-1] =-(Gamma-dx*(UU(r_e[i-1][1-1])-U(r_w[i-1][1-1])));//西側計算點
        a[index][index+1] =-(Gamma - dx*D(r_e[i-1][1-1]));//東側計算點
        //a[index][index-NX] =-;//南側計算點
        a[index][index+NX] =-(Gamma-dx*D(r_n[i-1][1-1]));//北側計算點
        //a[index][index-2] =-;//西西側計算點
        b[index] = dx*T_bottom + (Gamma-dx*UU(r_n[i-1][1-1])) * T[i][0] + (-dx*(-UU(r_w[i-1][1-1]))) * T[i-2][1];
    }
   
    //上邊界(i = 3~NX-1 ; j = NX)
    for(int i = 2; i <= NX-1; i++) {
        int index = (NX-1)*NX + i;
        a[index][index] =(3*Gamma+dx*(U(r_e[i-1][NX-1])-D(r_w[i-1][NX-1]))+dx*(U(r_n[i-1][NX-1])-D(r_s[i-1][NX-1]))) ;//本點
        a[index][index-1] =-(Gamma-dx*(UU(r_e[i-1][NX-1])-U(r_w[i-1][NX-1])));//西側計算點
        a[index][index+1] =-(Gamma-dx*(D(r_e[i-1][NX-1])));//東側計算點
        a[index][index-NX] =-(Gamma-dx*(UU(r_n[i-1][NX-1])-U(r_s[i-1][NX-1])));//南側計算點
        //a[index][index+NX] =-;//北側計算點
        //a[index][index-2] =-;//西西側計算點
        a[index][index-2*NX] =-(-dx*(-UU(r_s[i-1][NX-1])));//南南側計算點
        b[index] = 0.0 + (-dx*D(r_n[i-1][NX-1])) * T[i][NX+1] + (-dx*(-UU(r_w[i-1][NX-1])))*T[i-2][NX];    
    }
    //Break
   
    // 左邊界 (除角點外)
    for(int j = 2; j <= NX-1; j++) {
        int index = (j-1)*NX + 1;
        a[index][index] =4*Gamma + dx*U(r_e[1-1][j-1]) + dx*(U(r_n[1-1][j-1])-D(r_s[1-1][j-1])) ;//本點
        //a[index][index-1] =- ;//西側計算點
        a[index][index+1] =-(Gamma - dx*D(r_e[1-1][j-1])) ;//東側計算點
        a[index][index-NX] =-(Gamma - dx*(UU(r_n[1-1][j-1]) - U(r_s[1-1][j-1])) ) ; //南側計算點
        a[index][index+NX] =-(Gamma - dx*D(r_n[1-1][j-1]) ) ;//北側計算點
        //a[index][index-2*NX] =- ;//南南側計算點
        b[index] = dx*T_left + (Gamma - dx*UU(r_e[1-1][j-1])) * T[0][j] + (-dx*(-UU(r_s[1-1][j-1])) ) * T[1][j-2];
    }
   
    // 右邊界 (除角點外)
    for(int j = 2; j <= NX-1; j++) {
        int index = (j-1)*NX + NX;
        a[index][index] =(4*Gamma+dx*(-D(r_w[NX-1][j-1]))+dx*(U(r_n[NX-1][j-1])-D(r_s[NX-1][j-1]))) ; //本點
        a[index][index-1] =-(Gamma-dx*(-U(r_w[NX-1][j-1]))) ;//西側計算點
        //a[index][index+1] =- ;//東側計算點
        a[index][index-NX] =-(Gamma-dx*(UU(r_n[NX-1][j-1])-U(r_s[NX-1][j-1]))) ;//南側計算點
        a[index][index+NX] =-(Gamma-dx*(D(r_n[NX-1][j-1]))) ;//北側計算點
        a[index][index-2] =-(-dx*(-UU(r_w[NX-1][j-1]))) ;//西西側計算點
        //a[index][index-2*NX] =- ;//南南側計算點
        b[index] =-dx*T_right + (Gamma)*T[NX+1][j] + (-dx*(-UU(r_s[NX-1][j-1]))) * T[NX][j-2];
    }
   
   
    // 內點
    for(int i = 2; i <= NX-1; i++) {
        for(int j = 2; j <= NX-1; j++) {
            int index = NX*(j-1)+i;
            a[index][index] =4.0*Gamma + dx*(U(r_e[i-1][j-1]) - D(r_w[i-1][j-1])) + dx*(U(r_n[i-1][j-1]) - D(r_s[i-1][j-1]));//本點
            a[index][index+1] = -(Gamma - dx*(D(r_e[i-1][j-1]))) ;//東側計算點
            a[index][index-1] = -(Gamma - dx*(UU(r_e[i-1][j-1]) - U(r_w[i-1][j-1]))) ;//西側計算點
            a[index][index+NX] = -(Gamma - dx*(D(r_n[i-1][j-1]))) ;//北側計算點
            a[index][index-NX] = -(Gamma - dx*(UU(r_n[i-1][j-1]) - U(r_s[i-1][j-1]))) ;//南側計算點
            //a[index][index-2] =-;//西西側計算點
            //a[index][index-2*NX] =-;//南南側計算點
            b[index] = 0.0+dx*UU(r_w[i-1][j-1])*T[i-2][j] +dx*UU(r_s[i-1][j-1])*T[i][j-2] ;
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
     //更新溫度場
    // 將一維解轉換為二維溫度場
    for(int j = 1; j <= NX; j++) {
        for(int i = 1; i <= NX; i++) {
            T[i][j] = x[(j-1)*NX + i];
        }
    }
    //更新虛擬節點
    for(int j = 1 ; j <= NX ; j++){
        int index = (j-1)*NX + 1 ;
         T[0][j] = 2*T_left - x[index] ;
    }
    for(int j = 1 ; j <= NX ; j++){
        int index = (j-1)*NX + NX ;
         T[NX+1][j] = 2*T_right - x[index] ;
    }
    for(int i = 1 ; i <= NX ; i++){
        int index = (1-1)*NX + i ;
         T[i][0] = 2*T_bottom - x[index] ;
    }
    for(int i = 1 ; i <= NX ; i++){
        int index = (NX-1)*NX + i ;
         T[i][NX+1] = x[index] ;
    }
    //溫度場一更新完，下矩陣就馬上跟著更新
    for(int i = 1 ; i <= NX ; i++){
        for(int j = 1 ; j<= NX ; j++){
            if(i != 1){
            r_w[i-1][j-1] = (T[i-1][j]-T[i-2][j])/(T[i][j]-T[i-1][j]);
}
            if(i != NX){
            r_e[i-1][j-1] = (T[i][j]-T[i-1][j])/(T[i+1][j]-T[i][j]);
}
            if(j != 1){
            r_s[i-1][j-1] = (T[i][j-1]-T[i][j-2])/(T[i][j]-T[i][j-1]);
}
            r_n[i-1][j-1] = (T[i][j]-T[i][j-1])/(T[i][j+1]-T[i][j]);
        }
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
    ostringstream name;
    name << "12.TVD_MINMOD"<< NX+1 << "x" << NX+1 << "步數 = "<< setfill('0') << setw(6)  <<  m << ".vtk";
    ofstream out(name.str().c_str());
   
    // VTK 文件頭
    out << "# vtk DataFile Version 3.0\n";
    out << "QUICK_Gamma=0.0\n";
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
        a.assign(n+2, vector<double>(n+2, 0.0));
        b.assign(n+2, 0.0);
        x.assign(n+2, 0.0);
        x_old.assign(n+2, 0.0); //向量大小必須唯一
        T.assign(NX+2, vector<double>(NX+2, 0.0));
        r_w.assign(NX, vector<double>(NX, 1.0));
        r_e.assign(NX, vector<double>(NX, 1.0));
        r_s.assign(NX, vector<double>(NX, 1.0));
        r_n.assign(NX, vector<double>(NX, 1.0)); //四個矩陣初始值為0
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
        cout << "Pelect Number = " << Pelect << endl ;
        cout << "========================================" << endl;
        steadystate = false;
        initial(a , b , x) ; //初始化，係數矩陣，非齊性項，虛擬節點
        for(G = 0; G < max_G; G++) {
            Gamma0(a, b ,x);         //設定係數矩陣（含邊界條件）
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