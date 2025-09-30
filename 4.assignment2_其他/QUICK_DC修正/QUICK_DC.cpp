//利用空間二階精度QUICK的(DC修正)求解二維穩太擴散對流問題
//特點:引入(DC修正)
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
double a[n+2][n+2], b[n+2]; //a二維矩陣為等式左側係數矩陣;b一維矩陣為等式右側之外源項，方程組之非齊性項
double x[n+2], x_old[n+2] , x_old2[n+2] , x_nuse[n+2] , T[NX][NY] ; //x_nuse為未使用鬆弛技術之新迭代溫度場
const double Pelect = 1000 ;
const double  Gamma = (sqrt(2)/Pelect) ;//熱擴散係數 alpha
const double dx = 1.0/double(NX-1) ;  // 修正：網格間距應為1/(NX-1)
const double dy = 1.0/double(NY-1) ;  // 修正：網格間距應為1/(NY-1)
const double u = 1.0;  // x方向速度
const double v = 1.0;  // y方向速度
const double T_l = 1.0 ;
const double T_r = 0.0 ;
const double T_b = 0.0 ;
const double T_u = 1.0 ;
const double lamda = 0.3 ;
int G, max_G = 10000000000;
double maxerror; 
const float tolerance = 1e-10;
bool steadystate;

//初始化係數矩陣與原本QUICK格視同
//初始化矩陣//if(G == 0)
void initial(double a[][n+2], double b[], int n) {
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
}

//後續矩陣腹值 if(G>=0)
void DC_modified(double a[][n+2], double b[], int n) {
    ///////////////////////////////////////////////////////////////
    //定義符號:(東西)
    const double D_w = Gamma ; //for(i = 1 : NX-1)//擴散符號 
    const double D_wb = 2*Gamma ; //for(i = 0)
    const double D_e  = Gamma ; //for(i = 0 ; Nx-2)
    const double D_eb = 2*Gamma ;//for(i = Nx-1)
    const double F_w = dy*u ; //for(i = 0 : NX-1)//對流符號 
    const double F_e = dy*u ; //for(i = 0 : NX-1)
    const double a_W = D_w + F_w ; //for(i = 1 : NX-1)//線組係數 
    const double a_Wb = 0.0; //for(i = 0)
    const double a_E = D_e ; //for(i = 0 : NX-2)
    const double a_Eb = 0.0; //for(i = NX-1)
    const double S_ul = D_wb + F_w ; //for(i = 0)
    const double S_ur = D_eb - F_e ;
    
    //定義延遲修正等右外源項S_{uw}^{DC}
    vector<double> S_uwdc(n+2);
    S_uwdc[0] = 0.0 ;
    S_uwdc[n+1] = 0.0 ;  // 修正：應該是n+1而不是n+2
    for(int j = 1 ; j <= NY ; j++){
        for(int i = 1 ; i <= NX ; i++){
            if(i== 1){
                S_uwdc[NX*(j-1)+i] = 0.0 ;
            }else if(i == 2){
                S_uwdc[NX*(j-1)+i] = F_w*((3.0/8.0)*x_old[NX*(j-1)+i] - (1.0/8.0)*x_old[NX*(j-1)+i-1] - (2.0/8.0)*T_l) ;
            }else{
                S_uwdc[NX*(j-1)+i] = F_w*((3.0/8.0)*x_old[NX*(j-1)+i] - (2.0/8.0)*x_old[NX*(j-1)+i-1] - (1.0/8.0)*x_old[NX*(j-1)+i-2]) ;
            }
        }
    }
    
    //定義延遲修正等右外源項S_{ue}^{DC} 
    vector<double> S_uedc(n+2);
    S_uedc[0] = 0.0 ;
    S_uedc[n+1] = 0.0 ;  // 修正：應該是n+1而不是n+2
    for(int j = 1 ; j <= NY ; j++){
        for(int i = 1 ; i <= NX ; i++){
            if(i== 1){
                S_uedc[NX*(j-1)+i] = F_e*((3.0/8.0)*x_old[NX*(j-1)+i+1] - (1.0/8.0)*x_old[NX*(j-1)+i] - (2.0/8.0)*T_l) ;
            }else if(i == NX){
                S_uedc[NX*(j-1)+i] = 0.0 ;
            }else{
                S_uedc[NX*(j-1)+i] = F_e*((3.0/8.0)*x_old[NX*(j-1)+i+1] - (2.0/8.0)*x_old[NX*(j-1)+i] - (1.0/8.0)*x_old[NX*(j-1)+i-1]) ;  // 修正：F_w改為F_e
            }
        }
    }
    
    //定義符號:(南北) 
    const double D_s = Gamma ; //for(j = 1 : NY-1)//擴散符號 
    const double D_sb = 2*Gamma ; //for(j = 0)
    const double D_n  = Gamma ; //for(j = 0 ; NY-2)
    const double D_nb = 2*Gamma ;//for(j = NY-1)
    const double F_s = dx*v ; //for(j = 0 : NY-1)//對流符號 
    const double F_n = dx*v ; //for(j = 0 : NY-1)
    const double a_S = D_s + F_s ; //for(j = 1 : NY-1)//線組係數 
    const double a_Sb = 0.0; //for(j = 0)
    const double a_N = D_n ; //for(j = 0 : NY-2)
    const double a_Nb = 0.0; //for(j = NY-1)
    const double S_ub = D_sb + F_s ; //for(j = 0)//等右外源項 
    const double S_uu = D_nb - F_n;//for(j = NY) 
    
    //定義延遲修正等右外源項S_{us}^{DC}
    vector<double> S_usdc(n+2);
    S_usdc[0] = 0.0 ;
    S_usdc[n+1] = 0.0 ;  // 修正：應該是n+1而不是n+2
    for(int j = 1 ; j <= NY ; j++){
        for(int i = 1 ; i <= NX ; i++){
            if(j == 1){
                S_usdc[NX*(j-1)+i] = 0.0 ;
            }else if(j == 2){
                S_usdc[NX*(j-1)+i] = F_s*((3.0/8.0)*x_old[NX*(j-1)+i] - (1.0/8.0)*x_old[NX*(j-2)+i] - (2.0/8.0)*T_b) ;
            }else{
                S_usdc[NX*(j-1)+i] = F_s*((3.0/8.0)*x_old[NX*(j-1)+i] - (2.0/8.0)*x_old[NX*(j-2)+i] - (1.0/8.0)*x_old[NX*(j-3)+i]) ;  // 修正：S_uwdc改為S_usdc
            }
        }
    }
    
    //定義延遲修正等右外源項S_{un}^{DC} 
    vector<double> S_undc(n+2);
    S_undc[0] = 0.0 ;
    S_undc[n+1] = 0.0 ;  // 修正：應該是n+1而不是n+2
    for(int j = 1 ; j <= NY ; j++){
        for(int i = 1 ; i <= NX ; i++){
            if(j == 1){
                S_undc[NX*(j-1)+i] = F_n*((3.0/8.0)*x_old[NX*(j)+i] - (1.0/8.0)*x_old[NX*(j-1)+i] - (2.0/8.0)*T_b) ;
            }else if(j == NY){
                S_undc[NX*(j-1)+i] = 0.0 ;
            }else{
                S_undc[NX*(j-1)+i] = F_n*((3.0/8.0)*x_old[NX*(j)+i] - (2.0/8.0)*x_old[NX*(j-1)+i] - (1.0/8.0)*x_old[NX*(j-2)+i]) ;
            }
        }
    }
    
    // 角點處理
    //第一個點T_{0,0} -> a[1][1] (左下角)
    a[1][1] =   (a_Wb + a_E + F_e) + S_ul - F_w + (a_Sb + a_N + F_n) + S_ub - F_s ;
    a[1][2] =    -a_E ;  // T_{1,0}  
    a[1][1+NX] = -a_N ;  // T_{0,1}
    b[1] =       -S_uedc[1] + S_uwdc[1] -S_undc[1] + S_usdc[1] + S_ul * T_l + S_ub * T_b ;  // 左邊界=1.0

    //第三個點T_{NX-1,0} -> a[NX][NX] (右下角)
    a[NX][NX] =    (a_W + a_Eb - F_w) + S_ur + F_e + (a_Sb + a_N + F_n) + S_ub - F_s;
    a[NX][NX-1] =  -a_W ;  // T_{NX-2,0}
    a[NX][2*NX] =  -a_N ;  // T_{NX-1,1}
    b[NX] =        -S_uedc[NX] + S_uwdc[NX] -S_undc[NX] + S_usdc[NX] + S_ur * T_r + S_ub * T_b;  // 右邊界和下邊界都=0.0

    //第七個點T_{0,NY-1} -> a[n-NX+1][n-NX+1] (左上角)
    a[n-NX+1][n-NX+1] =   (a_Wb + a_E + F_e) + S_ul - F_w + (a_S + a_Nb - F_s) + S_uu + F_n ;
    a[n-NX+1][n-NX+2] =    -a_E ;  // T_{1,NY-1}
    a[n-NX+1][n-2*NX+1] =  -a_S ;  // T_{0,NY-2}
    b[n-NX+1] =          -S_uedc[n-NX+1] + S_uwdc[n-NX+1] -S_undc[n-NX+1] + S_usdc[n-NX+1] + S_ul * T_l + S_uu * T_u ;  // 左邊界和上邊界都=1.0

    //第九個點T_{NX-1,NY-1} -> a[n][n] (右上角)
    a[n][n] =      (a_W + a_Eb - F_w) + S_ur + F_e + (a_S + a_Nb - F_s) + S_uu + F_n  ;
    a[n][n-1] =    -a_W ;  // T_{NX-2,NY-1} 
    a[n][n-NX] =   -a_S ;  // T_{NX-1,NY-2}
    b[n] =        -S_uedc[n] + S_uwdc[n] -S_undc[n] + S_usdc[n] + S_ur * T_r + S_uu * T_u ;   

    /////////////////////////////////////////////////////////////////////
    // 邊界處理
    // 下邊界 (除角點外)
    for(int i = 2; i <= NX-1; i++) {
        a[i][i] =   (a_W + a_E - F_w + F_e)+(a_Sb + a_N + F_n) + S_ub - F_s ;  // 修正：內點應該用a_W而不是a_Wb
        a[i][i-1] = - a_W;
        a[i][i+1] = - a_E;
        a[i][i+NX] = - a_N;
        b[i] =  -S_uedc[i] + S_uwdc[i] -S_undc[i] + S_usdc[i]  + S_ub * T_b ;  // 修正：索引應該是i而不是NX * (i-1) + 1
    }
    
    // 上邊界 (除角點外)
    for(int i = n-NX+2; i <= n-1 ; i++) {
        a[i][i] =  (a_W + a_E - F_w + F_e) +(a_S + a_Nb - F_s) + S_uu + F_n ;  // 修正：內點應該用a_W而不是a_Wb
        a[i][i-1] = - a_W;
        a[i][i+1] = - a_E;
        a[i][i-NX] =  -a_S;
        b[i] = -S_uedc[i] + S_uwdc[i] -S_undc[i] + S_usdc[i]  + S_uu * T_u ;  // 修正：索引應該是i
    }
    
    // 左邊界 (除角點外)(上下內點) 
    for(int i = 2 ; i <= NY-1 ; i++) {
        int index = NX * (i-1) + 1;
        a[index][index] = (a_Wb + a_E + F_e)+(a_S + a_N - F_s + F_n) + S_ul - F_w ;
        a[index][index+1]  = -a_E ;//東點 
        a[index][index-NX] = -a_S;  // 南點
        a[index][index+NX] = -a_N;  // 北點
        b[index] = -S_uedc[index] + S_uwdc[index] -S_undc[index] + S_usdc[index]  + S_ul * T_l ;
    }
    
    // 右邊界 (除角點外)(上下內點) 
    for(int i = 2; i <= NY-1; i++) {
        int index = NX * i;
        a[index][index] =   (a_W + a_Eb - F_w) + (a_S + a_N - F_s + F_n) + S_ur + F_e;
        a[index][index-1] =  - a_W ;//西點  
        a[index][index-NX] = -a_S;  // 南點
        a[index][index+NX] = -a_N;  // 北點
        b[index] = -S_uedc[index] + S_uwdc[index] -S_undc[index] + S_usdc[index] + S_ur * T_r ;
    }
    
    // 內點
    for(int j = 2; j <= NY-1 ; j++) {
        for(int i = 2; i <= NX-1; i++) {
            int index = NX*(j-1)+i;
            a[index][index] = (a_W + a_E - F_w + F_e) + (a_S + a_N - F_s + F_n)  ;
            a[index][index-1] =-a_W ;  // 西點
            a[index][index+1] =-a_E ;  // 東點
            a[index][index-NX] =-a_S;  // 南點
            a[index][index+NX] =-a_N;  // 北點
            b[index] = -S_uedc[index] + S_uwdc[index] -S_undc[index] + S_usdc[index] ;  // 內點無熱源
        }
    }
}

void giveTvalue(int m){
    if(m == 0){//經過0次迭代,要進行第一次迭代 
    // 先複製當前解到x_old
    //以及賦值解到前二解 //112
        for(int k = 1; k <= n; k++) {
            x_old2[k] = x[k];
            x_old[k] = x[k];
        }
    }else if (m == 1){//經過一次迭代，要進行第二次迭代 
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
}

void Jacobi(double a[][n+2], double b[], double x[], int n) {
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
    for(int j = 0; j < NY; j++) {
        for(int i = 0; i < NX; i++) {
            T[i][j] = x[j*NX + i + 1]; //輸入x[] : 1~n 
        }
    }
    
    ostringstream name;
    name << "QUICK_DC" << setfill('0') << setw(6) << m << ".vtk";
    ofstream out(name.str().c_str());
    
    // VTK 文件頭
    out << "# vtk DataFile Version 3.0\n";
    out << "QUICK_DC\n";
    out << "ASCII\n";
    out << "DATASET STRUCTURED_POINTS\n";
    out << "DIMENSIONS " << NX << " " << NY << " 1\n";
    out << "ORIGIN 0 0 0\n";
    out << "SPACING " << dx << " " << dy << " 1\n";  // 使用實際網格間距
    out << "POINT_DATA " << NX * NY << "\n";
    
    // 輸出溫度場
    out << "SCALARS Temperature double 1\n";
    out << "LOOKUP_TABLE default\n";
    for(int j = 0; j < NY; j++) {
        for(int i = 0; i < NX; i++) {
            out << scientific << setprecision(6) << T[i][j] << "\n";
        }
    }
    
    out.close();
    cout << "VTK 文件已輸出: " << name.str() << endl;
}

int main() {
    cout << "程式碼開始執行...." << endl;
    cout << "網格大小: " << NX << " x " << NY << endl;
    cout << "網格間距: dx=" << dx << ", dy=" << dy << endl;
    cout << "邊界條件: 左邊界和上邊界=1.0, 右邊界和下邊界=0.0" << endl;
    steadystate = false;
    initial(a,b,n) ;
    giveTvalue(0) ;
    DC_modified(a,b,n) ;
    for(G = 0; G < max_G; G++) {
        if(G>0){
             giveTvalue(G) ;
             DC_modified(a,b,n) ;
        }
        Jacobi(a, b, x, n);
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
