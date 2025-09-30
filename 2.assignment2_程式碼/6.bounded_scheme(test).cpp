//利用有界格式(MINMOD)求解二維穩態熱擴散方程
//以矩陣描述通量限制器
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <iomanip>
#include <vector>
using namespace std ;
/*float A(float ) ;
float B(float ) ;
float C(float ) ;
float rf_w(int ,int ) ;
float rf_e(int ,int ) ;
float rf_s(int ,int ) ;
float rf_n(int ,int ) ;
void initial(float ,float ,int ) ;
void GaussSeidel(double , double , double , int) ;
void output(int ) ;*/
//會用到的參數
const int NX = 80;
const int NY = 80;
const int n = NX * NY;
double a[n+1][n+1], b[n+1]; //a二維矩陣為等式左側係數矩陣;b一維矩陣為等式右側之外源項，方程組之非齊性項
double x[n+1], x_old[n+1], T[NX][NY];
double rf_w[NX+1][NY+1] ,rf_e[NX+1][NY+1] ,rf_s[NX+1][NY+1] ,rf_n[NX+1][NY+1] ;//各個點在特定邊界的通量限制器函數值; 
const double Pelect = 10000 ;
const double  Gamma = 0.1 ;//熱擴散係數 alpha
const double dx = 1.0/double(NX-1) ;  // 修正：網格間距應為1/(NX-1)
const double dy = 1.0/double(NY-1) ;  // 修正：網格間距應為1/(NY-1)
const double u = 1.0;  // x方向速度
const double v = 1.0;  // y方向速度

int G, max_G = 3000;
double maxerror; 
const float tolerance = 1e-6;
bool steadystate; 
/////////////////////////////////////////////////////////////////////////////////////////
double A(double p){
	if(p>1){
		return 0.5 ;//通量限制器參數設定(p要輸入rf_{w,e,s,n}) 
	}else if(p<1 && (p>0 || p == 0)){
		return 0.0 ;
	}else{
		return 0.0 ;
	}
} 
double B(double p){
	if(p>1){
		return -0.5 ;//通量限制器參數設定(p要輸入rf_{w,e,s,n})
	}else if(p<1 && (p>0 || p == 0)){
		return 0.5 ;
	}else{
		return 0 ;
	}
} 
double C(double p){
	if(p>1){
		return 0.0 ;//通量限制器參數設定(p要輸入rf_{w,e,s,n})
	}else if(p<1 && (p>0 || p == 0)){
		return -0.5 ;
	}else{
		return 0 ;
	}
} 
/* 
//西邊界通量限制參數 
double rf_w(int i ,int j){//只討論i = 2 : NX , j = 1 : Ny //當i = 2時，rf_w需要特殊討論
    if(i == 2){
    	return 2*(x[(j-1)*NX+i-1] - 1)/(x[(j-1)*NX+i] - x[(j-1)*NX+i-1]) ;
	}else{
		return(x[(j-1)*NX+i-1] -x[(j-1)*NX+i-2])/(x[(j-1)*NX+i] - x[(j-1)*NX+i-1]) ;
	}//i = 1 ; 用不到 
}
//東邊界通量限制參數 
double rf_e(int i ,int j){//只討論i = 1 : NX-1 , j = 1 : Ny //當i = 1時，rf_e需要特殊討論
    if(i == 1){
    	return 2*(x[(j-1)*NX+i] - 1)/(x[(j-1)*NX+i+1] - x[(j-1)*NX+i]) ;
	}else{
		return(x[(j-1)*NX+i] -x[(j-1)*NX+i-1])/(x[(j-1)*NX+i+1] - x[(j-1)*NX+i]) ;
	}//i = NX ; 用不到 
}
//南邊界通量限制器函數
double rf_s(int i,int j){//只討論i = 1 : NX , j = 2 : Ny //當j = 2時，rf_s需要特殊討論 
	if(j== 2){
		return 2*(x[(j-2)*NX+i] - 1)/(x[(j-1)*NX+i] - x[(j-2)*NX+i]) ;
	}else{
		return (x[(j-2)*NX+i] - x[(j-3)*NX+i])/(x[(j-1)*NX+i] - x[(j-2)*NX+i]) ;
	}//j = 1 ; 用不到 
} 
//北邊界通量限制器函數
double rf_n(int i , int j){//只討論i = 1 : NX , j = 1 : NY-1 //當j=1時，rf_n需要特殊討論 
	if(j == 1){
		return 2*(x[(j-1)*NX+i] - 1)/(x[(j)*NX+i] - x[(j-1)*NX+i]) ;
	}else{
		return (x[(j-1)*NX+i] - x[(j-2)*NX+i])/(x[(j)*NX+i] - x[(j-1)*NX+i]) ;
	}//j = NY ; 用不到 
}*/
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//A對應到下游，B對應到上游，C對應到上上游之參數。
void initial(double a[][n+1] , double b[] , int n){
	// 初始化所有元素為0
    for(int i = 1; i <= n; i++) {
        for(int j = 1; j <= n; j++) {
            a[i][j] = 0.0;
        }
        b[i] = 0.0;//聯立方程式非齊性項 
        x[i] = 0.0;//傳輸變量，作為聯立方程式的解 
        x_old[i] = 0.0;//傳輸變量，作為聯立方程式的迭代前的解 
    }
    
    //初始化各點之特定邊的通量限制器函數值
	//初始化各點在西邊的通量限制器函數值 
	for(int i = 2 ; i <= NX ; i++){
		for(int j = 1 ; j <= NY ; j++){
			if(i == 2){
				rf_w[i][j] = 2*(x[(j-1)*NX+i-1] - 1)/(x[(j-1)*NX+i] - x[(j-1)*NX+i-1]) ;
			}else{
				rf_w[i][j] = (x[(j-1)*NX+i-1] -x[(j-1)*NX+i-2])/(x[(j-1)*NX+i] - x[(j-1)*NX+i-1]) ;
			}
		}
	} //i == 1 ; 用不到 
	//初始化各點在東邊的通量限制器函數值
	for(int i = 1 ; i <= NX-1 ; i++){
		for(int j = 1 ; j <= NY ; j++){
			if(i == 1){
				rf_e[i][j] = 2*(x[(j-1)*NX+i] - 1)/(x[(j-1)*NX+i+1] - x[(j-1)*NX+i]) ;
			}else{
				rf_e[i][j] = (x[(j-1)*NX+i] -x[(j-1)*NX+i-1])/(x[(j-1)*NX+i+1] - x[(j-1)*NX+i]) ;
			}
		}
	} //i == NX ; 用不到 
	//初始化各點在南邊的通量限制器函數值
	for(int i = 1 ; i <= NX ; i++){
		for(int j = 2 ; j <= NY ; j++){
			if(j == 2){
				rf_s[i][j] = 2*(x[(j-2)*NX+i] - 1)/(x[(j-1)*NX+i] - x[(j-2)*NX+i]) ;
			}else{
				rf_s[i][j] = (x[(j-2)*NX+i] - x[(j-3)*NX+i])/(x[(j-1)*NX+i] - x[(j-2)*NX+i]) ;
			}
		}
	} //j == 1 ; 用不到 
		//初始化各點在南邊的通量限制器函數值
	for(int i = 1 ; i <= NX ; i++){
		for(int j = 1 ; j <= NY-1 ; j++){
			if(j == 1){
				rf_n[i][j] = 2*(x[(j-1)*NX+i] - 1)/(x[(j)*NX+i] - x[(j-1)*NX+i]) ;
			}else{
				rf_n[i][j] = (x[(j-1)*NX+i] - x[(j-2)*NX+i])/(x[(j)*NX+i] - x[(j-1)*NX+i]) ;
			}
		}
	} //j == NY ; 用不到 
	
	
	
    //內點初始化(i = 3 : NX-1 ; j = 3 : NY-1 ) 
    //左右內點X上下內點 
    for(int i = 3 ; i <= NX-1 ; i++){
    	for(int j = 3 ; j <= NY - 1 ; j++){
    		a[(j-1)*NX+i][(j-1)*NX+i] = (2*Gamma + dy*(1+B(rf_e[i][j]-A(rf_w[i][j])))) + (2*Gamma + dx*(1+B(rf_n[i][j])-A(rf_s[i][j])));//本點
			a[(j-1)*NX+i][(j-1)*NX+i-1] = -Gamma + dy*(C(rf_e[i][j])-1-B(rf_w[i][j])) ;//W西鄰界計算點 
			a[(j-1)*NX+i][(j-1)*NX+i+1] = -Gamma + dy * A(rf_e[i][j]);//E東鄰界計算點
			a[(j-1)*NX+i][(j-2)*NX+i] =-Gamma + dx*(C(rf_n[i][j])-1-B(rf_s[i][j]));//S南鄰界計算點
			a[(j-1)*NX+i][(j)*NX+i] = -Gamma + dx* A(rf_n[i][j]);//N北鄰界計算點
			a[(j-1)*NX+i][(j-1)*NX+i-2] = -dy*C(rf_w[i][j]);//WW西西鄰界計算點
			a[(j-1)*NX+i][(j-3)*NX+i] = -dx*C(rf_s[i][j]);//SS南南鄰界計算點
			b[(j-1)*NX+i] = 0.0 ;//非齊性項 
			 
		}
	} 
    //六個邊界初始化(left\right\up\bottom\left2\bottom2) 
    //左邊界計算點
    //左右左邊界X上下內點 
	for(int j = 3 ; j <= NY-1 ; j++){//範圍:i = 1 ; j = 3 : NY-1 
		a[(j-1)*NX+1][(j-1)*NX+1] = (3*Gamma+dy*(1+B(rf_e[1][j])-C(rf_e[1][j])))+(2*Gamma + dx*(1+B(rf_n[1][j])-A(rf_s[1][j]))) ; //本點
		//a[(j-1)*NX+1][(j-1)*NX] = ? ;//不存在
		a[(j-1)*NX+1][(j-1)*NX+2] = (-Gamma+dy*A(rf_e[1][j])); //E東鄰界計算點 
		a[(j-1)*NX+1][(j-2)*NX+1] = -Gamma + dx*(C(rf_n[1][j])-1-B(rf_s[1][j]));//S南邊界計算點
		a[(j-1)*NX+1][(j)*NX+1] = -Gamma + dx* A(rf_n[1][j]);//N北邊界計算點
		a[(j-1)*NX+1][(j-3)*NX+1] = -dx*C(rf_s[1][j]) ;//SS南南邊界計算點
		b[(j-1)*NX+1] = -(-Gamma+dy*(2*C(rf_e[1][j])-1));
	} 
	//右邊界計算點
    //左右右邊界X上下內點
	for(int j = 3 ; j <= NY-1 ; j++){//範圍:i = NX ; j = 3 : NY-1 
		a[(j-1)*NX + NX][(j-1)*NX + NX] = (3*Gamma-dy*A(rf_w[NX][j]))+(2*Gamma + dx*(1+B(rf_n[NX][j])-A(rf_s[NX][j]))) ; //本點
		a[(j-1)*NX + NX][(j-1)*NX + NX-1]  = (-Gamma -dy*(B(rf_w[NX][j])+1));//西點
		a[(j-1)*NX + NX][(j-1)*NX] = -Gamma + dx*(C(rf_n[NX][j])-1-B(rf_s[NX][j]));//南點
		a[(j-1)*NX + NX][(j-1)*NX + 2*NX] = -Gamma + dx* A(rf_n[NX][j]);//北點
		a[(j-1)*NX + NX][(j-1)*NX + NX-2] = -dy * C(rf_w[NX][j]); //西西點 
	    a[(j-1)*NX + NX][(j-1)*NX -NX] = -dx*C(rf_s[NX][j]) ;//南南點 
	    b[(j-1)*NX + NX] = 0.0 ;//非齊性項 
	} 
	//南邊界計算點
	//左右內點X上下下邊界 
	for(int i = 3 ; i <= NX-1 ; i++){//範圍:i = 3 : NX-1 ; j = 1 
		a[i][i] = (2*Gamma + dy*(1+B(rf_e[i][1]))-A(rf_w[i][1]))+(3*Gamma+dx*(1+B(rf_n[i][1])-C(rf_n[i][1]))) ;//本點
		a[i][i-1] = -Gamma + dy*(C(rf_e[i][1])-1-B(rf_w[i][1])) ;//西點
		a[i][i+1] = -Gamma + dy * A(rf_e[i][1]);//東點
		a[i][i+NX] = (-Gamma+dx*A(rf_n[i][1]));//北點
		a[i][i-2] = -dy*C(rf_w[i][1]); //西西點
		b[i]  = 0.0 ; //非齊性項 
	}
	//北邊界計算點
	//左右內點X上下上邊界
	for(int i = 3 ;i <= NX-1 ; i++){//範圍:i = 3 : NX-1 ; j = NY 
		a[(NY-1)*NX+i][(NY-1)*NX+i] = (2*Gamma + dy*(1+B(rf_e[i][NY])-A(rf_w[i][NY])))+(3*Gamma-dx*A(rf_s[i][NY]));//本點
		a[(NY-1)*NX+i][(NY-1)*NX+i-1] =-Gamma + dy*(C(rf_e[i][NY])-1-B(rf_w[i][NY])) ;//西點
		a[(NY-1)*NX+i][(NY-1)*NX+i+1] =-Gamma + dy * A(rf_e[i][NY]);//東點
		a[(NY-1)*NX+i][(NY-2)*NX+i] = (-Gamma-dx*(B(rf_s[i][NY])+1)); //南點
		a[(NY-1)*NX+i][(NY-1)*NX+i-2] = -dy*C(rf_w[i][NY]); //西西點
		a[(NY-1)*NX+i][(NY-3)*NX+i]  = -dx*(C(rf_s[i][NY])); //南南點
		b[(NY-1)*NX+i] = -(-2*Gamma+dx);//非齊性項 
	} 
	//左邊界第二計算點
	//左右左第二排X上下內點
	for(int j = 3 ; j <= NY-1 ; j++){//範圍:i = 2 ; j = 3 : NY-1 
		a[(j-1)*NX+2][(j-1)*NX+2] = (2*Gamma+dy*(B(rf_e[2][j])+1-A(rf_w[2][j])))+(2*Gamma + dx*(1+B(rf_n[2][j])-A(rf_s[2][j]))); //本點
		a[(j-1)*NX+2][(j-1)*NX+1] = (-Gamma+dy*(C(rf_e[2][j])-1-B(rf_w[2][j])+C(rf_w[2][j]))); //西點
		a[(j-1)*NX+2][(j-1)*NX+3] = (-Gamma+dy*A(rf_e[2][j])); //東點
		a[(j-1)*NX+2][(j-2)*NX+2] = -Gamma + dx*(C(rf_n[2][j])-1-B(rf_s[2][j]));//南點
		a[(j-1)*NX+2][(j)*NX+2] = -Gamma + dx* A(rf_n[2][j]); //北點
		a[(j-1)*NX+2][(j-3)*NX+2] = -dx*C(rf_s[2][j]) ; //南南點
		b[(j-1)*NX+2] = 2*dy ;//非齊性項 
	} 
	//下邊界第二計算點
	//左右內點X上下下第二排 
	for(int i = 3 ; i <= NX-1 ; i++){//範圍: i = 3 : NX-1 ; j = 2
		a[NX + i][NX + i] = (2*Gamma + dy*(1+B(rf_e[i][2])-A(rf_w[i][2]))) + (2*Gamma+dx*(B(rf_n[i][2])+1-A(rf_s[i][2])));//本點
		a[NX + i][NX + i-1] = -Gamma + dy*(C(rf_e[i][2])-1-B(rf_w[i][2])); //西點 
		a[NX + i][NX + i+1] =-Gamma + dy * A(rf_e[i][2]) ; //東點
		a[NX + i][i] = (-Gamma+dx*(C(rf_n[i][2])-1-B(rf_s[i][2])+C(rf_s[i][2]))); //南點
		a[NX + i][2*NX + i]  = (-Gamma+dx*A(rf_n[i][2])); //北點
		a[NX + i][NX + i-2]  =-dy*C(rf_w[i][2]) ; //西西點
		b[NX+i] = 0.0 ; //非齊性項 
	}
	//九個角點初始化
	//第一個角點
	//左右左邊界X上下下邊界
	//範圍:i = 1 ; j = 1 ;
	a[1][1]  = (3*Gamma+dy*(1+B(rf_e[1][1])-C(rf_e[1][1])))+(3*Gamma+dx*(1+B(rf_n[1][1])-C(rf_n[1][1])))  ; //本點 
	a[1][2] = (-Gamma+dy*A(rf_e[1][1]));//東點 
	a[1][1+NX]= (-Gamma+dx*A(rf_n[1][1]));//北點
	b[1] = -(-Gamma+dy*(2*C(rf_e[1][1])-1)); //非齊性項 
	//第二個角點
	//左右左第二排邊界X上下下邊界
	//範圍:i = 2 ; j = 1 ;
	a[2][2]  = (2*Gamma+dy*(B(rf_e[2][1])+1-A(rf_w[2][1])))+(3*Gamma+dx*(1+B(rf_n[2][1])-C(rf_n[2][1])))  ; //本點 
	a[2][1] = (-Gamma+dy*(C(rf_e[2][1])-1-B(rf_w[2][1])+C(rf_w[2][1]))); //西點 
	a[2][3] = (-Gamma+dy*A(rf_e[2][1]));//東點 
	a[1][1+NX] = (-Gamma+dx*A(rf_n[2][1]));//北點
	b[2] = 2*dy ;//非齊性項
	//第三角點
	//左右右邊界X上下下邊界
	//範圍i = NX ; j = 1 ;
	a[NX][NX]  = (3*Gamma-dy*(A(rf_w[NX][1])))+(3*Gamma+dx*(1+B(rf_n[NX][1])-C(rf_n[NX][1])))  ; //本點 
	a[NX][NX-1] = (-Gamma-dy*(B(rf_w[NX][1])+1));//西點
	a[NX][2*NX] = (-Gamma+dx*A(rf_n[NX][1]));//北點
	a[NX][NX-2] = (-dy*C(rf_w[NX][1])); //西西點
	b[NX] = 0.0 ;//非齊性項  
	//第四角點
	//左右左邊界X上下下第二排
	//範圍:i = 1 ; j = 2 ;
	a[NX+1][NX+1]  = (3*Gamma+dy*(1+B(rf_e[1][2])-C(rf_e[1][2]))) + (2*Gamma+dx*(B(rf_n[1][2])+1-A(rf_s[1][2])));//本點 
	a[NX+1][NX+2] = (-Gamma+dy*A(rf_e[1][2])) ;//東點 
	a[NX+1][1] = (-Gamma+dx*(C(rf_n[1][2])-1-B(rf_s[1][2])+C(rf_s[1][2])));//南點 
	a[NX+1][2*NX+1] = (-Gamma+dx*A(rf_n[1][2]));//北點 
	b[NX+1] = -(-Gamma+dy*(2*C(rf_e[1][2])-1)); //非齊性項
	//第五角點
	//左右左第二排X上下下第二排
	//範圍:i = 2 ; j = 2 ;
	a[NX+2][NX+2] = (2*Gamma+dy*(B(rf_e[2][2])+1-A(rf_w[2][2]))) + (2*Gamma+dx*(B(rf_n[2][2])+1-A(rf_s[2][2])));//本點  
	a[NX+2][NX+1] = (-Gamma+dy*(C(rf_e[2][2])-1-B(rf_w[2][2])+C(rf_w[2][2]))); //西點 
	a[NX+2][NX+3] = (-Gamma+dy*A(rf_e[2][2]));//東點
	a[NX+2][2] = (-Gamma+dx*(C(rf_n[2][2])-1-B(rf_s[2][2])+C(rf_s[2][2])));//南點 
	a[NX+2][2*NX+2] = (-Gamma+dx*A(rf_n[2][2]));//北點
	b[NX+2] = 2*dy ;//非齊性項 
	//第六角點
	//左右右邊界X上下下第二排
	//範圍:i = NX ; j = 2 ;
	a[2*NX][2*NX] = (3*Gamma-dy*(A(rf_w[NX][2])))+ (2*Gamma+dx*(B(rf_n[NX][2])+1-A(rf_s[NX][2])));//本點
	a[2*NX][2*NX-1] = (-Gamma-dy*(B(rf_w[NX][2])+1));//西點
	a[2*NX][NX] = (-Gamma+dx*(C(rf_n[NX][2])-1-B(rf_s[NX][2])+C(rf_s[NX][2])));//南點 
	a[2*NX][3*NX] = (-Gamma+dx*A(rf_n[NX][2]));//北點
	a[2*NX][2*NX-2] = (-dy*C(rf_w[NX][2])); //西西點
	b[2*NX] = 0.0 ;//非齊性項 
    //第七角點
	//左右左邊界X上下上邊界
	//範圍:i = 1 ; j = NY ; 
	a[NX*(NY-1)+1][NX*(NY-1)+1] = (3*Gamma+dy*(1+B(rf_e[1][NY])-C(rf_e[1][NY])))+(3*Gamma-dx*A(rf_s[1][NY])) ;//本點
	a[NX*(NY-1)+1][NX*(NY-1)+2] = (-Gamma+dy*A(rf_e[1][NY])) ;//東點 
	a[NX*(NY-1)+1][NX*(NY-2)+1] = (-Gamma-dx*(B(rf_s[1][NY])+1)); //南點 
	a[NX*(NY-1)+1][NX*(NY-3)+1] = -dx*C(rf_s[1][NY]); //南南點 
	b[NX*(NY-1)+1] = -(-Gamma+dy*(2*C(rf_e[1][NY])-1))-(-2*Gamma+dx); //非齊性項
	//第八角點
	//左右左第二排X上下上邊界
	//範圍:i = 2 ; j = NY ;
	a[NX*(NY-1)+2][NX*(NY-1)+2] = (2*Gamma+dy*(B(rf_e[2][NY])+1-A(rf_w[2][NY])))+(3*Gamma-dx*A(rf_s[2][NY])) ;//本點
	a[NX*(NY-1)+2][NX*(NY-1)+1] = (-Gamma+dy*(C(rf_e[2][NY])-1-B(rf_w[2][NY])+C(rf_w[2][NY]))); //西點  
	a[NX*(NY-1)+2][NX*(NY-1)+3] = (-Gamma+dy*A(rf_e[2][NY]));//東點 
	a[NX*(NY-1)+2][NX*(NY-2)+2] = (-Gamma-dx*(B(rf_s[2][NY])+1)); //南點  
	a[NX*(NY-1)+2][NX*(NY-3)+2] = -dx*C(rf_s[2][NY]); //南南點  
	b[NX*(NY-1)+2] = 2*dy -(-2*Gamma+dx); //非齊性項
	//第九角點
	//左右右邊界X上下上邊界
	//範圍:i = NX ; j = NY ;
	a[NX*NY][NX*NY] = (3*Gamma-dy*(A(rf_w[NX][NY])))+(3*Gamma-dx*A(rf_s[NX][NY])) ;//本點
	a[NX*NY][NX*NY-1] = (-Gamma-dy*(B(rf_w[NX][NY])+1));//西點 
	a[NX*NY][NX*(NY-1)] = (-Gamma-dx*(B(rf_s[NX][NY])+1)); //南點  
	a[NX*NY][NX*NY-2] = (-dy*C(rf_w[NX][NY])); //西西點
	a[NX*NY][NX*(NY-2)] = -dx*C(rf_s[NX][NY]); //南南點  
	b[NX*NY] = 0.0 -(-2*Gamma+dx); //非齊性項
}
/////////////////////////////////////////////////////////////////////////////
void newcoefficient(double a[][n+1] , double b[] , int n){
    //更新各點之特定邊的通量限制器函數值
	//更新各點在西邊的通量限制器函數值 
	for(int i = 2 ; i <= NX ; i++){
		for(int j = 1 ; j <= NY ; j++){
			if(i == 2){
				rf_w[i][j] = 2*(x[(j-1)*NX+i-1] - 1)/(x[(j-1)*NX+i] - x[(j-1)*NX+i-1]) ;
			}else{
				rf_w[i][j] = (x[(j-1)*NX+i-1] -x[(j-1)*NX+i-2])/(x[(j-1)*NX+i] - x[(j-1)*NX+i-1]) ;
			}
		}
	} //i == 1 ; 用不到 
	//更新各點在東邊的通量限制器函數值
	for(int i = 1 ; i <= NX-1 ; i++){
		for(int j = 1 ; j <= NY ; j++){
			if(i == 1){
				rf_e[i][j] = 2*(x[(j-1)*NX+i] - 1)/(x[(j-1)*NX+i+1] - x[(j-1)*NX+i]) ;
			}else{
				rf_e[i][j] = (x[(j-1)*NX+i] -x[(j-1)*NX+i-1])/(x[(j-1)*NX+i+1] - x[(j-1)*NX+i]) ;
			}
		}
	} //i == NX ; 用不到 
	//更新各點在南邊的通量限制器函數值
	for(int i = 1 ; i <= NX ; i++){
		for(int j = 2 ; j <= NY ; j++){
			if(j == 2){
				rf_s[i][j] = 2*(x[(j-2)*NX+i] - 1)/(x[(j-1)*NX+i] - x[(j-2)*NX+i]) ;
			}else{
				rf_s[i][j] = (x[(j-2)*NX+i] - x[(j-3)*NX+i])/(x[(j-1)*NX+i] - x[(j-2)*NX+i]) ;
			}
		}
	} //j == 1 ; 用不到 
	//更新各點在南邊的通量限制器函數值
	for(int i = 1 ; i <= NX ; i++){
		for(int j = 1 ; j <= NY-1 ; j++){
			if(j == 1){
				rf_n[i][j] = 2*(x[(j-1)*NX+i] - 1)/(x[(j)*NX+i] - x[(j-1)*NX+i]) ;
			}else{
				rf_n[i][j] = (x[(j-1)*NX+i] - x[(j-2)*NX+i])/(x[(j)*NX+i] - x[(j-1)*NX+i]) ;
			}
		}
	} //j == NY ; 用不到 
	
	
	
    //內點更新(i = 3 : NX-1 ; j = 3 : NY-1 ) 
    //左右內點X上下內點 
    for(int i = 3 ; i <= NX-1 ; i++){
    	for(int j = 3 ; j <= NY - 1 ; j++){
    		a[(j-1)*NX+i][(j-1)*NX+i] = (2*Gamma + dy*(1+B(rf_e[i][j]-A(rf_w[i][j])))) + (2*Gamma + dx*(1+B(rf_n[i][j])-A(rf_s[i][j])));//本點
			a[(j-1)*NX+i][(j-1)*NX+i-1] = -Gamma + dy*(C(rf_e[i][j])-1-B(rf_w[i][j])) ;//W西鄰界計算點 
			a[(j-1)*NX+i][(j-1)*NX+i+1] = -Gamma + dy * A(rf_e[i][j]);//E東鄰界計算點
			a[(j-1)*NX+i][(j-2)*NX+i] =-Gamma + dx*(C(rf_n[i][j])-1-B(rf_s[i][j]));//S南鄰界計算點
			a[(j-1)*NX+i][(j)*NX+i] = -Gamma + dx* A(rf_n[i][j]);//N北鄰界計算點
			a[(j-1)*NX+i][(j-1)*NX+i-2] = -dy*C(rf_w[i][j]);//WW西西鄰界計算點
			a[(j-1)*NX+i][(j-3)*NX+i] = -dx*C(rf_s[i][j]);//SS南南鄰界計算點
			b[(j-1)*NX+i] = 0.0 ;//非齊性項 
			 
		}
	} 
    //六個邊界更新(left\right\up\bottom\left2\bottom2) 
    //左邊界計算點
    //左右左邊界X上下內點 
	for(int j = 3 ; j <= NY-1 ; j++){//範圍:i = 1 ; j = 3 : NY-1 
		a[(j-1)*NX+1][(j-1)*NX+1] = (3*Gamma+dy*(1+B(rf_e[1][j])-C(rf_e[1][j])))+(2*Gamma + dx*(1+B(rf_n[1][j])-A(rf_s[1][j]))) ; //本點
		//a[(j-1)*NX+1][(j-1)*NX] = ? ;//不存在
		a[(j-1)*NX+1][(j-1)*NX+2] = (-Gamma+dy*A(rf_e[1][j])); //E東鄰界計算點 
		a[(j-1)*NX+1][(j-2)*NX+1] = -Gamma + dx*(C(rf_n[1][j])-1-B(rf_s[1][j]));//S南邊界計算點
		a[(j-1)*NX+1][(j)*NX+1] = -Gamma + dx* A(rf_n[1][j]);//N北邊界計算點
		a[(j-1)*NX+1][(j-3)*NX+1] = -dx*C(rf_s[1][j]) ;//SS南南邊界計算點
		b[(j-1)*NX+1] = -(-Gamma+dy*(2*C(rf_e[1][j])-1));
	} 
	//右邊界計算點
    //左右右邊界X上下內點
	for(int j = 3 ; j <= NY-1 ; j++){//範圍:i = NX ; j = 3 : NY-1 
		a[(j-1)*NX + NX][(j-1)*NX + NX] = (3*Gamma-dy*A(rf_w[NX][j]))+(2*Gamma + dx*(1+B(rf_n[NX][j])-A(rf_s[NX][j]))) ; //本點
		a[(j-1)*NX + NX][(j-1)*NX + NX-1]  = (-Gamma -dy*(B(rf_w[NX][j])+1));//西點
		a[(j-1)*NX + NX][(j-1)*NX] = -Gamma + dx*(C(rf_n[NX][j])-1-B(rf_s[NX][j]));//南點
		a[(j-1)*NX + NX][(j-1)*NX + 2*NX] = -Gamma + dx* A(rf_n[NX][j]);//北點
		a[(j-1)*NX + NX][(j-1)*NX + NX-2] = -dy * C(rf_w[NX][j]); //西西點 
	    a[(j-1)*NX + NX][(j-1)*NX -NX] = -dx*C(rf_s[NX][j]) ;//南南點 
	    b[(j-1)*NX + NX] = 0.0 ;//非齊性項 
	} 
	//南邊界計算點
	//左右內點X上下下邊界 
	for(int i = 3 ; i <= NX-1 ; i++){//範圍:i = 3 : NX-1 ; j = 1 
		a[i][i] = (2*Gamma + dy*(1+B(rf_e[i][1]))-A(rf_w[i][1]))+(3*Gamma+dx*(1+B(rf_n[i][1])-C(rf_n[i][1]))) ;//本點
		a[i][i-1] = -Gamma + dy*(C(rf_e[i][1])-1-B(rf_w[i][1])) ;//西點
		a[i][i+1] = -Gamma + dy * A(rf_e[i][1]);//東點
		a[i][i+NX] = (-Gamma+dx*A(rf_n[i][1]));//北點
		a[i][i-2] = -dy*C(rf_w[i][1]); //西西點
		b[i]  = 0.0 ; //非齊性項 
	}
	//北邊界計算點
	//左右內點X上下上邊界
	for(int i = 3 ;i <= NX-1 ; i++){//範圍:i = 3 : NX-1 ; j = NY 
		a[(NY-1)*NX+i][(NY-1)*NX+i] = (2*Gamma + dy*(1+B(rf_e[i][NY])-A(rf_w[i][NY])))+(3*Gamma-dx*A(rf_s[i][NY]));//本點
		a[(NY-1)*NX+i][(NY-1)*NX+i-1] =-Gamma + dy*(C(rf_e[i][NY])-1-B(rf_w[i][NY])) ;//西點
		a[(NY-1)*NX+i][(NY-1)*NX+i+1] =-Gamma + dy * A(rf_e[i][NY]);//東點
		a[(NY-1)*NX+i][(NY-2)*NX+i] = (-Gamma-dx*(B(rf_s[i][NY])+1)); //南點
		a[(NY-1)*NX+i][(NY-1)*NX+i-2] = -dy*C(rf_w[i][NY]); //西西點
		a[(NY-1)*NX+i][(NY-3)*NX+i]  = -dx*(C(rf_s[i][NY])); //南南點
		b[(NY-1)*NX+i] = -(-2*Gamma+dx);//非齊性項 
	} 
	//左邊界第二計算點
	//左右左第二排X上下內點
	for(int j = 3 ; j <= NY-1 ; j++){//範圍:i = 2 ; j = 3 : NY-1 
		a[(j-1)*NX+2][(j-1)*NX+2] = (2*Gamma+dy*(B(rf_e[2][j])+1-A(rf_w[2][j])))+(2*Gamma + dx*(1+B(rf_n[2][j])-A(rf_s[2][j]))); //本點
		a[(j-1)*NX+2][(j-1)*NX+1] = (-Gamma+dy*(C(rf_e[2][j])-1-B(rf_w[2][j])+C(rf_w[2][j]))); //西點
		a[(j-1)*NX+2][(j-1)*NX+3] = (-Gamma+dy*A(rf_e[2][j])); //東點
		a[(j-1)*NX+2][(j-2)*NX+2] = -Gamma + dx*(C(rf_n[2][j])-1-B(rf_s[2][j]));//南點
		a[(j-1)*NX+2][(j)*NX+2] = -Gamma + dx* A(rf_n[2][j]); //北點
		a[(j-1)*NX+2][(j-3)*NX+2] = -dx*C(rf_s[2][j]) ; //南南點
		b[(j-1)*NX+2] = 2*dy ;//非齊性項 
	} 
	//下邊界第二計算點
	//左右內點X上下下第二排 
	for(int i = 3 ; i <= NX-1 ; i++){//範圍: i = 3 : NX-1 ; j = 2
		a[NX + i][NX + i] = (2*Gamma + dy*(1+B(rf_e[i][2])-A(rf_w[i][2]))) + (2*Gamma+dx*(B(rf_n[i][2])+1-A(rf_s[i][2])));//本點
		a[NX + i][NX + i-1] = -Gamma + dy*(C(rf_e[i][2])-1-B(rf_w[i][2])); //西點 
		a[NX + i][NX + i+1] =-Gamma + dy * A(rf_e[i][2]) ; //東點
		a[NX + i][i] = (-Gamma+dx*(C(rf_n[i][2])-1-B(rf_s[i][2])+C(rf_s[i][2]))); //南點
		a[NX + i][2*NX + i]  = (-Gamma+dx*A(rf_n[i][2])); //北點
		a[NX + i][NX + i-2]  =-dy*C(rf_w[i][2]) ; //西西點
		b[NX+i] = 0.0 ; //非齊性項 
	}
	//九個角點更新
	//第一個角點
	//左右左邊界X上下下邊界
	//範圍:i = 1 ; j = 1 ;
	a[1][1]  = (3*Gamma+dy*(1+B(rf_e[1][1])-C(rf_e[1][1])))+(3*Gamma+dx*(1+B(rf_n[1][1])-C(rf_n[1][1])))  ; //本點 
	a[1][2] = (-Gamma+dy*A(rf_e[1][1]));//東點 
	a[1][1+NX]= (-Gamma+dx*A(rf_n[1][1]));//北點
	b[1] = -(-Gamma+dy*(2*C(rf_e[1][1])-1)); //非齊性項 
	//第二個角點
	//左右左第二排邊界X上下下邊界
	//範圍:i = 2 ; j = 1 ;
	a[2][2]  = (2*Gamma+dy*(B(rf_e[2][1])+1-A(rf_w[2][1])))+(3*Gamma+dx*(1+B(rf_n[2][1])-C(rf_n[2][1])))  ; //本點 
	a[2][1] = (-Gamma+dy*(C(rf_e[2][1])-1-B(rf_w[2][1])+C(rf_w[2][1]))); //西點 
	a[2][3] = (-Gamma+dy*A(rf_e[2][1]));//東點 
	a[1][1+NX] = (-Gamma+dx*A(rf_n[2][1]));//北點
	b[2] = 2*dy ;//非齊性項
	//第三角點
	//左右右邊界X上下下邊界
	//範圍i = NX ; j = 1 ;
	a[NX][NX]  = (3*Gamma-dy*(A(rf_w[NX][1])))+(3*Gamma+dx*(1+B(rf_n[NX][1])-C(rf_n[NX][1])))  ; //本點 
	a[NX][NX-1] = (-Gamma-dy*(B(rf_w[NX][1])+1));//西點
	a[NX][2*NX] = (-Gamma+dx*A(rf_n[NX][1]));//北點
	a[NX][NX-2] = (-dy*C(rf_w[NX][1])); //西西點
	b[NX] = 0.0 ;//非齊性項  
	//第四角點
	//左右左邊界X上下下第二排
	//範圍:i = 1 ; j = 2 ;
	a[NX+1][NX+1]  = (3*Gamma+dy*(1+B(rf_e[1][2])-C(rf_e[1][2]))) + (2*Gamma+dx*(B(rf_n[1][2])+1-A(rf_s[1][2])));//本點 
	a[NX+1][NX+2] = (-Gamma+dy*A(rf_e[1][2])) ;//東點 
	a[NX+1][1] = (-Gamma+dx*(C(rf_n[1][2])-1-B(rf_s[1][2])+C(rf_s[1][2])));//南點 
	a[NX+1][2*NX+1] = (-Gamma+dx*A(rf_n[1][2]));//北點 
	b[NX+1] = -(-Gamma+dy*(2*C(rf_e[1][2])-1)); //非齊性項
	//第五角點
	//左右左第二排X上下下第二排
	//範圍:i = 2 ; j = 2 ;
	a[NX+2][NX+2] = (2*Gamma+dy*(B(rf_e[2][2])+1-A(rf_w[2][2]))) + (2*Gamma+dx*(B(rf_n[2][2])+1-A(rf_s[2][2])));//本點  
	a[NX+2][NX+1] = (-Gamma+dy*(C(rf_e[2][2])-1-B(rf_w[2][2])+C(rf_w[2][2]))); //西點 
	a[NX+2][NX+3] = (-Gamma+dy*A(rf_e[2][2]));//東點
	a[NX+2][2] = (-Gamma+dx*(C(rf_n[2][2])-1-B(rf_s[2][2])+C(rf_s[2][2])));//南點 
	a[NX+2][2*NX+2] = (-Gamma+dx*A(rf_n[2][2]));//北點
	b[NX+2] = 2*dy ;//非齊性項 
	//第六角點
	//左右右邊界X上下下第二排
	//範圍:i = NX ; j = 2 ;
	a[2*NX][2*NX] = (3*Gamma-dy*(A(rf_w[NX][2])))+ (2*Gamma+dx*(B(rf_n[NX][2])+1-A(rf_s[NX][2])));//本點
	a[2*NX][2*NX-1] = (-Gamma-dy*(B(rf_w[NX][2])+1));//西點
	a[2*NX][NX] = (-Gamma+dx*(C(rf_n[NX][2])-1-B(rf_s[NX][2])+C(rf_s[NX][2])));//南點 
	a[2*NX][3*NX] = (-Gamma+dx*A(rf_n[NX][2]));//北點
	a[2*NX][2*NX-2] = (-dy*C(rf_w[NX][2])); //西西點
	b[2*NX] = 0.0 ;//非齊性項 
    //第七角點
	//左右左邊界X上下上邊界
	//範圍:i = 1 ; j = NY ; 
	a[NX*(NY-1)+1][NX*(NY-1)+1] = (3*Gamma+dy*(1+B(rf_e[1][NY])-C(rf_e[1][NY])))+(3*Gamma-dx*A(rf_s[1][NY])) ;//本點
	a[NX*(NY-1)+1][NX*(NY-1)+2] = (-Gamma+dy*A(rf_e[1][NY])) ;//東點 
	a[NX*(NY-1)+1][NX*(NY-2)+1] = (-Gamma-dx*(B(rf_s[1][NY])+1)); //南點 
	a[NX*(NY-1)+1][NX*(NY-3)+1] = -dx*C(rf_s[1][NY]); //南南點 
	b[NX*(NY-1)+1] = -(-Gamma+dy*(2*C(rf_e[1][NY])-1))-(-2*Gamma+dx); //非齊性項
	//第八角點
	//左右左第二排X上下上邊界
	//範圍:i = 2 ; j = NY ;
	a[NX*(NY-1)+2][NX*(NY-1)+2] = (2*Gamma+dy*(B(rf_e[2][NY])+1-A(rf_w[2][NY])))+(3*Gamma-dx*A(rf_s[2][NY])) ;//本點
	a[NX*(NY-1)+2][NX*(NY-1)+1] = (-Gamma+dy*(C(rf_e[2][NY])-1-B(rf_w[2][NY])+C(rf_w[2][NY]))); //西點  
	a[NX*(NY-1)+2][NX*(NY-1)+3] = (-Gamma+dy*A(rf_e[2][NY]));//東點 
	a[NX*(NY-1)+2][NX*(NY-2)+2] = (-Gamma-dx*(B(rf_s[2][NY])+1)); //南點  
	a[NX*(NY-1)+2][NX*(NY-3)+2] = -dx*C(rf_s[2][NY]); //南南點  
	b[NX*(NY-1)+2] = 2*dy -(-2*Gamma+dx); //非齊性項
	//第九角點
	//左右右邊界X上下上邊界
	//範圍:i = NX ; j = NY ;
	a[NX*NY][NX*NY] = (3*Gamma-dy*(A(rf_w[NX][NY])))+(3*Gamma-dx*A(rf_s[NX][NY])) ;//本點
	a[NX*NY][NX*NY-1] = (-Gamma-dy*(B(rf_w[NX][NY])+1));//西點 
	a[NX*NY][NX*(NY-1)] = (-Gamma-dx*(B(rf_s[NX][NY])+1)); //南點  
	a[NX*NY][NX*NY-2] = (-dy*C(rf_w[NX][NY])); //西西點
	a[NX*NY][NX*(NY-2)] = -dx*C(rf_s[NX][NY]); //南南點  
	b[NX*NY] = 0.0 -(-2*Gamma+dx); //非齊性項
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
        //從1~N得到新的一組溫度場 
    }
    //重新賦值各點各邊之A(rf_w)\A(rf_e)\A(rf_s)\A(rf_n)\B(rf_w)\B(rf_e)\B(rf_s)\B(rf_n)\C(rf_w)\C(rf_e)\C(rf_s)\C(rf_n)
    //更新各點各邊之通量限制器函數
	 
} 
//////////////////////////////////////////////////////////////////////////////////////////////////
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
    name << "bounded_scheme_Pe = " << Pelect <<","<<"Grid "<< nodes_i << "x" << nodes_j << ","<< setfill('0') << setw(6) << m << ".vtk";
    ofstream out(name.str().c_str());
    
    // VTK 文件頭：使用 STRUCTURED_GRID
    out << "# vtk DataFile Version 3.0\n";
    out << "bounded_scheme\n";
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
//主程式 
int main() {
    cout << "MINMOD程式碼開始執行...." << endl;
    cout << "網格大小: " << NX << " x " << NY << endl;
    cout << "網格間距: dx=" << dx << ", dy=" << dy << endl;
    cout << "邊界條件: 左邊界和上邊界=1.0, 右邊界和下邊界=0.0" << endl;
    cout << "目前Pelect number = " << Pelect << ",Gamma = " << Gamma << endl ; 
    steadystate = false;
    initial(a, b, n);
    
    for(G = 0; G < max_G; G++) {
        if(G > 0) {
            GaussSeidel(a, b, x, n);
            newcoefficient(a ,b ,n) ;
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

///////////////////////////////////////////////////////////////////////////////////
//問題:要如何邊迭代邊修正係數 
//在每一輪迭代中，A\B\C係數需要重新賦值 