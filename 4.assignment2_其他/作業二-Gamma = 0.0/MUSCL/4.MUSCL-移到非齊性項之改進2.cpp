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
int Nx[] = {80,80};
//有限體積法採用Wet node處理，故計算點數目 = 網格數目, 節點數目 = 網格數目+1
int NX; //NX = NY => dx = dy
int n;

vector<vector<double> > T;
vector<vector<double> > T_old ;
double dx; //一個網格間距
//邊界條件定義
const double T_left = 1.0;    //左邊界Dirichlet條件
const double T_right = 0.0;   //右邊界Dirichlet條件  
const double T_bottom = 0.0;  //下邊界Dirichlet條件
//迭代參數
const double lamda = 0.1; //逐次超鬆弛迭代法
double TotalVarition ;
int G, max_G = 10000000;
const double w = 1e-6 ;
double maxerror ;
const float tolerance = 1e-10;
bool steadystate;
double MUSCL_w(int i,int j){
	double T_D = T[i][j] ;
	double T_U = T[i-1][j] ;
	double T_UU = T[i-2][j];
	double x = (T_U-T_UU)/(T_D-T_U) ;
    return 0.3 ;
}
double MUSCL_e(int i,int j){
	double T_D = T[i+1][j] ;
	double T_U = T[i][j] ; 
	double T_UU = T[i-1][j];
	double x = (T_U-T_UU)/(T_D-T_U) ;
	return 0.3 ;
}
double MUSCL_s(int i,int j){
	double T_D = T[i][j] ;
	double T_U = T[i][j-1] ;
	double T_UU = T[i][j-2];
	double x = (T_U-T_UU)/(T_D-T_U) ;
	return 0.3 ;
}
double MUSCL_n(int i,int j){
	double T_D = T[i][j+1] ;
	double T_U = T[i][j] ;
	double T_UU = T[i][j-1];
	double x = (T_U-T_UU)/(T_D-T_U) ;
    return 0.3 ;
}
//TVD-MUSCL格式為在東西方向為三點格式,南北方向為三點格式

//顯示迭代法（考慮Neumann邊界條件）
void Gamma0() { //因為係數式時變的,每一次迭代都必須要做.....
    if(G == 0){
    	for(int j = 1 ; j <= NX ; j++){
    		T[0][j]=  2*T_left-T[1][j] ;
		}
		for(int i = 1 ; i <= NX ; i++){
			T[i][0] = 2*T_bottom - T[i][1] ;
		}
	}
	for(int i = 0 ; i <= NX+1 ; i++){
		for(int j = 0 ; j <= NX+1 ; j++){
			T_old[i][j] = T[i][j] ;
		}
	}
    double T_new1 =  ((dx*T_left + dx*T_bottom) - T[2][1]*dx*(MUSCL_e(1,1))/(2.0) - T[1][2]*dx*(MUSCL_n(1,1))/2.0)/ (dx*(1-(MUSCL_e(1,1))) + dx*(1-(MUSCL_n(1,1))/2.0)) ;
    T[1][1] = T_old[1][1] + lamda * (T_new1 - T_old[1][1]);
    T[0][1] = 2*T_left-T[1][1] ;
    T[1][0] = 2*T_bottom-T[1][1] ;
       

    double T_new2 =  ((-dx*T_right + dx*T_bottom) + T[NX-1][1]*dx*(1-(MUSCL_w(NX,1))/2.0) - T[NX][2]*dx*(MUSCL_n(NX,1))/2.0 )/((-dx*(MUSCL_w(NX,1))/2.0) + dx*(1-(MUSCL_n(NX,1))/2.0));
    T[NX][1] = T_old[NX][1] + lamda * (T_new2 - T_old[NX][1]);
    T[NX+1][1] = 0.0 ;
    T[NX][0] = 2*T_bottom-T[NX][1] ;
    

    double T_new3 =  (dx*T_left - T[2][NX]*dx*(MUSCL_e(1,NX))/2.0 + T[1][NX-1]*(dx*MUSCL_s(1,NX)/2.0) )/(dx*(1-(MUSCL_e(1,NX)/2.0)) + dx*(1-MUSCL_s(1,NX)/2.0)) ;
    T[1][NX] = T_old[1][NX] + lamda * (T_new3 - T_old[1][NX]);
	T[0][NX] = 2*T_left-T[1][NX] ; 
	T[1][NX+1] = 0.0 ;   
    int index4 = (NX-1)*NX+NX ;  
	
     
    double T_new4 =  ((-dx*T_right) + T[NX-1][NX]*(dx*(1-MUSCL_w(NX,NX)/2.0)) + T[NX][NX-1]*(dx*(1-MUSCL_s(NX,NX)/2.0)))/(dx*(-MUSCL_w(NX,NX)/2.0) + dx*(1-MUSCL_s(NX,NX)/2.0)) ;
    T[NX][NX] = T_old[NX][NX] + lamda * (T_new4 - T_old[NX][NX]);
	T[NX][NX+1] = 0.0 ; 
	T[NX+1][NX] = 0.0 ;          
    // 邊界處理（不含角點）
    // 下邊界 (除角點外)(i = 2~NX-1 ; j = 1)
    for(int i = 2 ; i <= NX-1; i++) {

        double T_new = (dx*T_bottom - T[i+1][1]*dx*MUSCL_e(i,1)/2.0 + T[i-1][1]*(dx*(1-MUSCL_w(i,1)/2.0)) -T[i][2]*dx*MUSCL_n(i,1)/2.0 )/(dx*(1-MUSCL_e(i,1)/2.0-MUSCL_w(i,1)/2.0)+dx*(1-MUSCL_n(i,1)/2.0));
        T[i][1] = T_old[i][1] + lamda * (T_new - T_old[i][1]);
        T[i][0] = 2*T_bottom-T[i][1] ;
    }
    //上邊界(i = 2~NX-1 ; j = NX)
    for(int i = 2 ; i <= NX-1; i++) {

        double T_new = (0.0 - T[i+1][NX]*dx*MUSCL_e(i,NX)/2.0 + T[i-1][NX]*(dx*(1-MUSCL_w(i,NX)/2.0)) + T[i][NX-1]*(dx*(1-MUSCL_s(i,NX)/2.0)))/(dx*(1-MUSCL_e(i,NX)/2.0-MUSCL_w(i,NX)/2.0)+dx*(1-MUSCL_s(i,NX)/2.0)) ;
        T[i][NX] = T_old[i][NX] + lamda * (T_new - T_old[i][NX]);
        T[i][NX+1] = 0.0 ;
    }

    // 左邊界 (除角點外)
    for(int j = 2; j <= NX-1; j++) {

        double T_new = (dx*T_left - T[2][j]*dx*MUSCL_e(1,j)/2.0 - T[1][j+1]*dx*MUSCL_n(1,j)/2.0 +T[1][j-1]*(dx*(1-MUSCL_s(1,j)/2.0)) )/(dx*(1-MUSCL_e(1,j)/2.0)+dx*(1-MUSCL_n(1,j)/2.0-MUSCL_s(1,j)/2.0)) ;
        T[1][j] = T_old[1][j] + lamda * (T_new - T_old[1][j]);
		T[0][j] = 2*T_left-T[1][j] ; 
    }
   
    // 右邊界 (除角點外)
    for(int j = 2; j <= NX-1; j++) {
    	
  
        double T_new = (-dx*T_right +T[NX-1][j]*(dx*(1-MUSCL_w(NX,j)/2.0)) -T[NX][j+1]*dx*MUSCL_n(NX,j)/2.0 +T[NX][j-1]*(dx*(1-MUSCL_s(NX,j)/2.0)) )/(dx*(-MUSCL_w(NX,j)/2.0)+dx*(1-MUSCL_n(NX,j)/2.0-MUSCL_s(NX,j)/2.0)) ;
        T[NX][j] = T_old[NX][j] + lamda * (T_new - T_old[NX][j]);
        T[NX+1][j] = 0.0 ;
    }
    // 內點
    for(int i = 2; i <= NX-1; i++) {
        for(int j = 2; j <= NX-1; j++) {

            double T_new = (0.0 - T[i+1][j]*dx*MUSCL_e(i,j)/2.0 + T[i-1][j]*(dx*(1-MUSCL_w(i,j)/2.0)) - T[i+1][j+1]*dx*MUSCL_n(i,j)/2.0 + T[i+1][j-1]*(dx*(1-MUSCL_s(i,j)/2.0)))/(dx*(1-MUSCL_e(i,j)/2.0-MUSCL_w(i,j)/2.0)+dx*(1-MUSCL_n(i,j)/2.0-MUSCL_s(i,j)/2.0)) ;
            T[i][j] = T_old[i][j] + lamda * (T_new - T_old[i][j]);
        }
    }
}

void MaxError(int NX) {
    //每一次迭代完成都更新溫度場
    //計算當前步最大誤差
    maxerror = 0;
    for(int i = 1; i <= NX; i++) {
    	for(int j = 1 ; j <= NX ; j++){
    		double error = fabs(T[i][j] - T_old[i][j]);
            if(maxerror < error) {
                maxerror = error;
            }
		}
    }//計算i = 1 : NX ; j = 1 : NX-1 的最大誤差
    double TotalVariationx = 0 ;
    double TotalVariationy = 0 ;
    for(int i = 1 ; i <= NX-1 ; i++){
    	for(int j = 1 ; j<=NX ; j++){
    		TotalVariationx += fabs(T[i][j]-T[i+1][j]) ;
		}
	}
	for(int i = 1 ; i <= NX ; i++){
		for(int j = 1 ; j <= NX ; j++){
			TotalVariationy += fabs(T[i][j]-T[i][j+1]) ;
		}
	}
	TotalVarition = TotalVariationx + TotalVariationy ;
	if(G % 500 == 0){
		cout << "迭代步數:" << G << "(TV) = " <<  TotalVarition << endl ; 
	}
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
        T.assign(NX+2, vector<double>(NX+2, 0.0));
        T_old.assign(NX+2, vector<double>(NX+2, 0.0));
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
        for(G = 0; G < max_G; G++) {
            Gamma0();  //設定係數矩陣（含邊界條件）
            MaxError(NX);  //計算最大誤差
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