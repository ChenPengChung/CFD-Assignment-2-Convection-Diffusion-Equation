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
int Nx[] = {40,80};
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
const double T_bottom = 0.0;  //下邊界Dirichlet條件
const double Pelect = 1000.0 ;
const double Gamma = sqrt(2.0)/Pelect ;//熱擴散係數
//迭代參數
const double lamda = 0.0005; //逐次超鬆弛迭代法
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
    //左邊界點
    for(int j = 1 ; j <= NX ; j++){
        int index = (j-1)*NX + 1 ;
        a[index][index] = 5*Gamma + dx*(U(r_e[1-1][j-1])-UU(r_e[1-1][j-1])) + dx*(D(r_n[1-1][j-1])-U(r_s[1-1][j-1]));
        a[index][index+1] =-Gamma + dx*(D(r_e[1-1][j-1]));
        a[index][index+NX] = -Gamma+dx*(D(r_n[1-1][j-1]));
        a[index][index-NX] = -Gamma+dx*(UU(r_n[1-1][j-1])-U((r_s[1-1][j-1])));
        a[index][index-2*NX] = dx*(-U(r_s[1-1][j-1]));
        b[index] = T_left*(2*Gamma-dx*(2*UU(r_e[1-1][j-1])-1));
    }
    //左第二排
    for(int j = 3 ; j <= NX-1 ; j++){
        int index = (j-1)*NX + 2 ;
        a[index][index] = 4*Gamma + dx*(U(r_e[2-1][j-1]) - D(r_w[2-1][j-1])) + dx*(U(r_n[2-1][j-1])-D(r_s[2-1][j-1])) ;
        a[index][index-1] = -Gamma + dx*(UU(r_e[2-1][j-1])-U(r_w[2-1][j-1])+UU(r_w[2-1][j-1])) ;
        a[index][index+1] = -Gamma + dx*(D(r_e[2-1][j-1])) ;
        a[index][index-NX] = -Gamma + dx*(UU(r_n[2-1][j-1])-U(r_s[2-1][j-1])) ;
        a[index][index-NX] = -Gamma + dx*(D(r_n[2-1][j-1])) ;
        a[index][index-2*NX] = -dx*UU(r_s[2-1][j-1]) ;
        b[index] = T_left*dx*2*UU(r_w[2-1][j-1]) ; 
    }
    //右邊界點
    for(int j = 3 ; j <= NX-1 ; j++){
        int index = (j-1)*NX + NX ;
        a[index][index] = 3*Gamma + dx*(U(r_e[NX-1][j-1])-D(r_w[NX-1][j-1])) + dx*(U(r_n[NX-1][j-1])-D(r_s[NX-1][j-1])) + dx*(D(r_e[NX-1][j-1]))  ;
        a[index][index-1] = -Gamma + dx*(UU(r_e[NX-1][j-1])-U(r_w[NX-1][j-1])) ;
        a[index][index-2] = -dx*UU(r_w[NX-1][j-1]) ;
        a[index][index+NX] = -Gamma + dx*D(r_n[NX-1][j-1]) ;
        a[index][index-NX] = -Gamma+dx*(UU(r_n[NX-1][j-1])-U(r_s[NX-1][j-1])) ;
        a[index][index-2*NX] = dx*(-UU(r_s[NX-1][j-1])) ;
        b[index] = 0.0 ;
    }
    //下邊界點
    for(int i = 3 ; i <= NX-1 ; i++){
        int index = (1-1)*NX + i ;
        a[index][index] = 5*Gamma + dx*(U(r_e[i-1][1-1])-D(r_w[i-1][1-1])) + dx*(U(r_n[i-1][1-1])-UU(r_n[i-1][1-1])) ;
        a[index][index-1] = -Gamma + dx*(UU(r_e[i-1][1-1])-U(r_w[i-1][1-1])) ;
        a[index][index+1] = -Gamma + dx*(D(r_e[i-1][1-1])) ;
        a[index][index+NX] = -Gamma + dx*(D(r_n[i-1][1-1])) ;
        a[index][index-2] = -dx*(UU(r_w[i-1][1-1])) ;
        b[index] = T_bottom * (2*Gamma-(2*(UU(r_n[i-1][1-1]))-1)) ;
     }
     //下邊界第二排
     for(int i = 3 ; i <= NX-1 ; i++){
        int index = (2-1)*NX+i ;
        a[index][index] = 4*Gamma + dx*(U(r_e[i-1][2-1])-D(r_w[i-1][2-1])) + dx*(U(r_n[i-1][2-1])-D(r_s[i-1][2-1])) ;
        a[index][index-1] = -Gamma + dx*(UU(r_e[i-1][2-1])-U(r_w[i-1][2-1])) ;
        a[index][index+1] = -Gamma +dx*(D(r_e[i-1][2-1])) ;
        a[index][index-2] = dx*(-UU(r_w[i-1][2-1])) ;
        a[index][index+NX] = -Gamma + dx*(D(r_n[i-1][2-1])) ;
        a[index][index-NX] = -Gamma + dx*(UU(r_n[i-1][2-1])-U(r_s[i-1][2-1])+UU(r_s[i-1][2-1])) ;
        b[index] = T_bottom * dx * 2 * (-UU(r_s[i-1][2-1])) ;
     }
    //上邊界第一排
    for(int i = 3  ; i <= NX-1 ; i++){
        int index = (NX-1)*NX + i ;
        a[index][index] = 3*Gamma + dx*(U(r_e[i-1][NX-1])-D(r_w[i-1][NX-1])) + dx*(D(r_n[i-1][NX-1])+U(r_n[i-1][NX-1])-D(r_s[i-1][NX-1])) ;
        a[index][index+1] = -Gamma +dx*(D(r_e[i-1][NX-1])) ;
        a[index][index-1] = -Gamma + dx*(UU(r_e[i-1][NX-1])-U(r_w[i-1][NX-1])) ;
        a[index][index-NX] = -Gamma + dx*(UU(r_n[i-1][NX-1])-U(r_s[i-1][NX-1])) ;
        a[index][index-2] = dx*(-UU(r_w[i-1][NX-1])) ;
        a[index][index-2*NX] = dx*(-UU(r_s[i-1][NX-1])) ;
        b[index] = 0.0 ;
    }
    //第一角點
    int p1 = (1-1)*NX+1 ;
    a[p1][p1] = 6*Gamma + dx*(U(r_e[1-1][1-1])-UU(r_e[1-1][1-1])) + dx*(U(r_n[1-1][1-1])-UU(r_n[1-1][1-1])) ;
    a[p1][p1+1] = -Gamma + dx*D(r_e[1-1][1-1]) ;
    a[p1][p1+1] = -Gamma + dx*D(r_n[1-1][1-1]) ;
    b[p1] = -T_left*dx*(2*UU(r_e[1-1][1-1])-1)-dx*T_bottom*(2*UU(r_n[1-1][1-1])-1) + 2*Gamma*(T_left+T_bottom) ;
    //第二角點
    int p2 = (1-1)*NX + 2 ;
    a[p2][p2] = 5*Gamma + dx*(U(r_e[2-1][1-1])-D(r_w[2-1][1-1])) + dx*(U(r_n[2-1][1-1])-UU(r_n[2-1][1-1])) ;
    a[p2][p2-1] = Gamma + dx*(UU(r_e[2-1][1-1])-U(r_w[2-1][1-1])+UU(r_w[2-1][1-1])) ;
    a[p2][p2+1] = -Gamma + dx*D(r_e[2-1][1-1]) ;
    a[p2][p2+NX] = -Gamma + dx*D(r_n[2-1][1-1]) ;
    b[p2] = 2*Gamma*T_bottom + dx * T_left * (2*UU(r_w[2-1][1-1])) -dx*T_bottom*(2*UU(r_n[2-1][1-1])-1) ;
    //第三角點
    int p3 = (1-1)*NX+NX ;
    a[p3][p3] = 4*Gamma +dx*(D(r_e[NX-1][1-1])+U(r_e[NX-1][1-1])-D(r_w[NX-1][1-1])) + dx*(U(r_n[NX-1][1-1])-UU(r_n[NX-1][1-1])) ;
    a[p3][p3-1] = -Gamma + dx*(UU(r_e[NX-1][1-1])-U(r_w[NX-1][1-1])) ;
    a[p3][p3-2] = -dx*(UU(r_w[NX-1][1-1])) ;
    a[p3][p3+NX] = -Gamma + dx*(D(r_n[NX-1][1-1])) ;
    b[p3] = 2*Gamma*T_bottom -dx * T_bottom * (2*UU(r_n[NX-1][1-1])-1) ;
    //第四角點
    int p4 = (2-1)*NX +1 ;
    a[p4][p4] = 5*Gamma + dx*(U(r_e[1-1][2-1])-UU(r_e[1-1][2-1])) + dx*(U(r_n[1-1][2-1])-D(r_s[1-1][2-1])) ;
    a[p4][p4+1] = -Gamma +dx*(D(r_e[1-1][2-1])) ;
    a[p4][p4-NX] = -Gamma +dx*(UU(r_n[1-1][2-1])-U(r_s[1-1][2-1])+UU(r_s[1-1][2-1])) ;
    a[p4][p4+NX] = -Gamma +dx*(D(r_n[1-1][2-1])) ;
    b[p4] = 2*Gamma*T_left -dx*T_left*(2*UU(r_e[1-1][2-1])-1) -dx*T_bottom*(2*UU(r_s[1-1][2-1])) ;
    //第五角點
    int p5 = (2-1)*NX + 2 ;
    a[p5][p5] = 4*Gamma + dx*(U(r_e[2-1][2-1]) - D(r_w[2-1][2-1])) + dx*(U(r_n[2-1][2-1]) - D(r_s[2-1][2-1])) ;
    a[p5][p5-1] = -Gamma + dx*(UU(r_e[2-1][2-1])-U(r_w[2-1][2-1])+UU(r_w[2-1][2-1])) ;
    a[p5][p5-NX] = -Gamma + dx*(UU(r_n[2-1][2-1])-U(r_s[2-1][2-1])+UU(r_s[2-1][2-1])) ;
    a[p5][p5+1] = -Gamma + dx*(D(r_e[2-1][2-1])) ;
    a[p5][p5+NX] = -Gamma + dx*(D(r_n[2-1][2-1])) ;
    b[p5] = dx*T_bottom*(2*UU(r_s[2-1][2-1])) + dx*T_left*(2*UU(r_w[2-1][2-1]));
    //第六角點
    int p6 = (2-1)*NX + NX ;
    a[p6][p6] = 3*Gamma + dx*(D(r_e[NX-1][2-1])+U(r_e[NX-1][2-1])-D(r_w[NX-1][2-1])) + dx*(U(r_n[NX-1][2-1])-D(r_s[NX-1][2-1])) ;
    a[p6][p6-1] = -Gamma + dx*(UU(r_e[NX-1][2-1])-U(r_w[NX-1][2-1])) ;
    a[p6][p6-2] = dx*(-UU(r_w[NX-1][2-1])) ;
    a[p6][p6-NX] = -Gamma + dx*(UU(r_n[NX-1][2-1])-U(r_s[NX-1][2-1])+UU(r_s[NX-1][2-1])) ;
    a[p6][p6+NX] = -Gamma +dx*(D(r_n[NX-1][2-1])) ;
    b[p6] = dx*T_bottom*(2*UU(r_s[NX-1][2-1])) ;
    //第七角點
    int p7 = (NX-1)*NX + 1 ;
    a[p7][p7] = 4*Gamma + dx*(U(r_e[1-1][NX-1])-UU(r_e[1-1][NX-1])) + dx*(D(r_n[1-1][NX-1])+U(r_n[1-1][NX-1])-D(r_s[1-1][NX-1])) ;
    a[p7][p7-1] = -Gamma + dx*(D(r_e[1-1][NX-1])) ;
    a[p7][p7-NX] = -Gamma +dx*(UU(r_n[1-1][NX-1])-U(r_s[1-1][NX-1])) ;
    a[p7][p7-2*NX] = dx*(-UU(r_s[1-1][NX-1])) ;
    b[p7] = T_left*2*Gamma  - dx*T_left*(2*UU(r_e[1-1][NX-1])-1) ;
    //第八角點
    int p8 = (NX-1)*NX + 2 ;
    a[p8][p8] = 3*Gamma + dx*(U(r_e[2-1][NX-1])-D(r_w[2-1][NX-1])) +dx*(D(r_n[2-1][NX-1])+U(r_n[2-1][NX-1])-D(r_s[2-1][NX-1])) ;
    a[p8][p8-1] = -Gamma + dx*(UU(r_e[2-1][NX-1])-U(r_w[2-1][NX-1])+UU(r_w[2-1][NX-1])) ;
    a[p8][p8+1]  = -Gamma + dx*(D(r_e[2-1][NX-1])) ;
    a[p8][p8-NX] = -Gamma + dx*(UU(r_n[2-1][NX-1])-U(r_s[2-1][NX-1])) ;
    a[p8][p8-2*NX] = dx*(-UU(r_s[2-1][NX-1])) ;
    b[p8] = dx*T_left*(2*UU(r_w[2-1][NX-1])) ;
    //第九角點
    int p9 = (NX-1)*NX + NX ;
    a[p9][p9] = 2*Gamma + dx*(D(r_e[NX-1][NX-1])+U(r_e[NX-1][NX-1])-D(r_w[NX-1][NX-1]))+ dx*(D(r_n[NX-1][NX-1])+U(r_n[NX-1][NX-1])-D(r_s[NX-1][NX-1]));
    a[p9][p9-1] = -Gamma + dx*(UU(r_e[NX-1][NX-1])-U(r_w[NX-1][NX-1])) ;
    a[p9][p9-NX] = -Gamma + dx*(UU(r_n[NX-1][NX-1])-U(r_s[NX-1][NX-1])) ;
    a[p9][p9-2] = dx*(-UU(r_w[NX-1][NX-1])) ;
    a[p9][p9-2*NX] = dx*(-UU(r_s[NX-1][NX-1])) ;
    b[p9] = 0.0 ;
    // 內點
    for(int i = 3; i <= NX-1; i++) {
        for(int j = 3; j <= NX-1; j++) {
            int index = NX*(j-1)+i;
            a[index][index] =4.0*Gamma + dx*(U(r_e[i-1][j-1]) - D(r_w[i-1][j-1])) + dx*(U(r_n[i-1][j-1]) - D(r_s[i-1][j-1]));//本點
            a[index][index+1] = -(Gamma - dx*(D(r_e[i-1][j-1]))) ;//東側計算點
            a[index][index-1] = -(Gamma - dx*(UU(r_e[i-1][j-1]) - U(r_w[i-1][j-1]))) ;//西側計算點
            a[index][index+NX] = -(Gamma - dx*(D(r_n[i-1][j-1]))) ;//北側計算點
            a[index][index-NX] = -(Gamma - dx*(UU(r_n[i-1][j-1]) - U(r_s[i-1][j-1]))) ;//南側計算點
            a[index][index-2] =-dx*UU(r_w[i-1][j-1])  ;//西西側計算點
            a[index][index-2*NX] =-dx*UU(r_s[i-1][j-1]);//南南側計算點
            b[index] = 0.0 ;
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
            T[i-1][j-1] = x[(j-1)*NX + i];
        }
    }
    //溫度場一更新完，下矩陣就馬上跟著更新
    for(int i = 1 ; i <= NX ; i++){
        for(int j = 1 ; j<= NX ; j++){
            if(i != 1){
            	if(i == 2){
            		r_w[i-1][j-1] = (T[i-2][j]-(2*T_left-T[0][j-1]))/(T[i-1][j-1]-T[i-2][j]);
				}else{
            r_w[i-1][j-1] = (T[i-2][j]-T[i-3][j])/(T[i-1][j-1]-T[i-2][j]);}
            }

            if(i == 1){
                r_e[i-1][j-1] = (T[i-1][j-1]-(2*T_left-T[0][j-1]))/(T[i-1+1][j-1]-T[i-1][j-1]);
			}else if(i == NX){
                r_e[i-1][j-1] = -1.0 ;
			}else{
            r_e[i-1][j-1] = (T[i-1][j-1]-T[i-2][j-1])/(T[i-1+1][j-1]-T[i-1][j-1]);}
            
            if(j != 1){
            	if(j == 2){
            		r_s[i-1][j-1] = (T[i-1][j-1-1]-(2*T_bottom-T[i-1][0]))/(T[i-1][j-1]-T[i-1][j-1-1]);
				}else{
            r_s[i-1][j-1] = (T[i-1][j-1-1]-T[i-1][j-2-1])/(T[i-1][j-1]-T[i-1][j-1-1]);}
            }


            if(j == 1){
            	r_n[i-1][j-1] = (T[i-1][j-1]-(2*T_bottom-T[i-1][0]))/(T[i-1][j+1-1]-T[i-1][j-1]);
			}else if(j == NX){
            	r_n[i-1][j-1] = -1.0 ;
			}else{
            	r_n[i-1][j-1] = (T[i-1][j-1]-T[i-1][j-1-1])/(T[i-1][j+1-1]-T[i-1][j-1]);}
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
    // 將一維解轉換為二維溫度場
    for(int j = 1; j <= NX; j++) {
        for(int i = 1; i <= NX; i++) {
            T[i-1][j-1] = x[(j-1)*NX+1] ; //輸入x[] : 1~n 
        }
    }
    
    int nodes_i = NX + 1;  // x方向節點數
    int nodes_j = NX + 1;  // y方向節點數
    int cells_i = NX;      // x方向網格數
    int cells_j = NX;      // y方向網格數
    
    ostringstream name;
    name << "13.TVD_MUSCL schemePe = " << Pelect <<","<<"Grid "<< nodes_i << "x" << nodes_j << ","<< setfill('0') << setw(6) << m << ".vtk";
    ofstream out(name.str().c_str());
    
    // VTK 文件頭：使用 STRUCTURED_GRID
    out << "# vtk DataFile Version 3.0\n";
    out << "13.TVD_MUSCL scheme\n";
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
        n = (NX)*(NX);
        // 重新調整向量大小
        a.assign(n+2, vector<double>(n+2, 0.0));
        b.assign(n+2, 0.0);
        x.assign(n+2, 0.0);
        x_old.assign(n+2, 0.0); //向量大小必須唯一
        T.assign(NX, vector<double>(NX, 0.0));
        r_w.assign(NX, vector<double>(NX, 1.0));
        r_e.assign(NX, vector<double>(NX, 1.0));
        r_s.assign(NX, vector<double>(NX, 1.0));
        r_n.assign(NX, vector<double>(NX, 1.0)); //四個矩陣初始值為0
        cout << "========================================" << endl;
        cout << "c++ programming is coming...." << endl;
        cout << "Mesh : " << NX+1 << " x " << NX+1 << endl;
        cout << "Grid spacing: dx = dy = " << dx << endl;
        cout << "Boundary conditions:" << endl;
        cout << "  left boundary(Dirichlet): T = " << T_left << endl;
        cout << "  bottom boundary(Dirichlet): T = " << T_bottom << endl;
        cout << "  top & right  boundary(Neumann):  T_[i][NY] = T_[i][NY+1](ghost point)" << endl;
        cout << "Relaxation factor: λ = " << lamda << endl;
        cout << "Convergence criterion: " << tolerance << endl;
        cout << "Peclet Number = " << Pelect << endl ;
        cout << "========================================" << endl;
        steadystate = false;
        initial(a , b , x) ; //初始化，係數矩陣，非齊性項，虛擬節點
        for(G = 0; G < max_G; G++) {
            Gamma0(a, b ,x);         //設定係數矩陣（含邊界條件）
            Jacobi(a, b, x, n);  //執行一次Jacobi迭代
            if(G % 1000 == 0) {
                cout << "itteration = " << G;
                if(G > 0) {
                    cout << ", the Max error  = " << scientific << maxerror;
                }
                cout << endl;
                output(G);
            }
           
            if(G > 100 && maxerror < tolerance) {
                steadystate = true;
                cout << "\nhave already converged!" << endl;
                cout << "Final iteration: " << G << endl;
                cout << "Final error: " << scientific << maxerror << endl;
                break;
            }
        }
       
        if(!steadystate) {
            cout << "\nWarning: Reached maximum iterations, but not converged!" << endl;
        }
       
        output(G);  //輸出最終結果
        cout << "Mesh " << NX+1 << "x" << NX+1 << " computation completed\n" << endl;
    }
    cout << "All computations completed!" << endl;
    return 0;
}
