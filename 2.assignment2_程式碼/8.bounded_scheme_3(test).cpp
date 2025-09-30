//利用有界格式(MINMOD)求解二維穩態熱擴散方程
//bounded scheme 3 加入初始cd模擬 
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <iomanip>
#include <vector>
#define W 1e-10 
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
const int NX = 40;
const int NY = 40;
const int n = NX * NY;
double a[n+1][n+1], b[n+1]; //a二維矩陣為等式左側係數矩陣;b一維矩陣為等式右側之外源項，方程組之非齊性項
double x[n+1], x_old[n+1], T[NX][NY];
const double Pelect = 1.0 ;
const double  Gamma = sqrt(2.0)/Pelect + W ;//熱擴散係數 alpha
const double dx = 1.0/double(NX-1) ;  // 修正：網格間距應為1/(NX-1)
const double dy = 1.0/double(NY-1) ;  // 修正：網格間距應為1/(NY-1)
const double T_Bl = 1.0 ;//邊界溫度left  = 1.0 ;
const double T_Bu = 1.0 ;//邊界溫度up = 1.0 ; 
const double a_W = (Gamma + dy/2.0);
const double a_E = (Gamma - dy/2.0);
const double a_S = (Gamma + dx/2.0);
const double a_N = (Gamma - dx/2.0);
int G, max_G = 8000;
double maxerror ;//第0次迭代開始計算...... 
const double tolerance = 1e-10;
bool steadystate; 
/////////////////////////////////////////////////////////////////////////////////////////

double A(double p){
	if(p>1){
		return 0.5 ;//通量限制器參數設定(p要輸入rf_{w,e,s,n}) 
	}else if(p<1 && p>= 0.0 ){
		return 0.0 ;
	}else{
		return 0.0 ;
	}
} 
double B(double p){
	if(p>1){
		return -0.5 ;//通量限制器參數設定(p要輸入rf_{w,e,s,n})
	}else if(p<1 && p>= 0.0 ){
		return 0.5 ;
	}else{
		return 0.0 ;
	}
} 
double C(double p){
	if(p>1){
		return 0.0 ;//通量限制器參數設定(p要輸入rf_{w,e,s,n})
	}else if(p<1 && p>= 0.0 ){
		return -0.5 ;
	}else{
		return 0.0 ;
	}
} 
//西邊界通量限制參數 
double rf_w(int i ,int j){//只討論i = 2 : NX , j = 1 : Ny //當i = 2時，rf_w需要特殊討論
    if(i == 2){
    	return 2*(x[(j-1)*NX+i-1] - T_Bl)/((x[(j-1)*NX+i] - x[(j-1)*NX+i-1])+W) ;
	}else{
		return(x[(j-1)*NX+i-1] -x[(j-1)*NX+i-2])/((x[(j-1)*NX+i] - x[(j-1)*NX+i-1])+W) ;
	}//i = 1 ; 用不到 
}
//東邊界通量限制參數 
double rf_e(int i ,int j){//只討論i = 1 : NX-1 , j = 1 : Ny //當i = 1時，rf_e需要特殊討論
    if(i == 1){
    	return 2*(x[(j-1)*NX+i] - T_Bl)/((x[(j-1)*NX+i+1] - x[(j-1)*NX+i])+W) ;
	}else{
		return(x[(j-1)*NX+i] -x[(j-1)*NX+i-1])/((x[(j-1)*NX+i+1] - x[(j-1)*NX+i])+W) ;
	}//i = NX ; 用不到 
}
//南邊界通量限制器函數
double rf_s(int i,int j){//只討論i = 1 : NX , j = 2 : Ny //當j = 2時，rf_s需要特殊討論 
	if(j== 2){
		return 2*(x[(j-2)*NX+i] )/((x[(j-1)*NX+i] - x[(j-2)*NX+i])+W) ;
	}else{
		return (x[(j-2)*NX+i] - x[(j-3)*NX+i])/((x[(j-1)*NX+i] - x[(j-2)*NX+i])+W) ;
	}//j = 1 ; 用不到 
} 
//北邊界通量限制器函數
double rf_n(int i , int j){//只討論i = 1 : NX , j = 1 : NY-1 //當j=1時，rf_n需要特殊討論 
	if(j == 1){
		return 2*(x[(j-1)*NX+i] )/((x[(j)*NX+i] - x[(j-1)*NX+i])+W) ;
	}else{
		return (x[(j-1)*NX+i] - x[(j-2)*NX+i])/((x[(j)*NX+i] - x[(j-1)*NX+i])+W) ;
	}//j = NY ; 用不到 
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//A對應到下游，B對應到上游，C對應到上上游之參數。
void initial(double a[][n+1] , double b[] , int n){
	// 初始化所有元素為0
    // 初始化所有元素為0
    for(int i = 0; i <= n; i++) {
        for(int j = 0; j <= n; j++) {
            a[i][j] = 0.0;
        }
        b[i] = 0.0;
        x[i] = 0.0; // 初始化解向量
    }
}
/////////////////////////////////////////////////////////////////////////////
void orderCD(double a[][n+1]  ,double b[] , int n){
	// 設定邊界條件和係數矩陣
    // 四個角點
    a[1][1] = (6.0*Gamma) + (dx/2.0) + (dy/2.0) ; 
    a[1][2] = -a_E;
    a[1][NX+1] = -a_N;
    b[1] = 2.0*a_W; // 邊界溫度
    
    a[NX][NX] = (6.0*Gamma) - (dx/2.0) + (dy/2.0);
    a[NX][NX-1] = -a_W;
    a[NX][2*NX] = -a_N;
    b[NX] = 0.0;// 邊界溫度
    
    a[n-NX+1][n-NX+1] = (6.0*Gamma) + (dx/2.0) - (dy/2.0);
    a[n-NX+1][n-NX+2] =  -a_E;
    a[n-NX+1][n-2*NX+1] = -a_S;
    b[n-NX+1] = 2.0*(a_W+a_N); //邊界效應引起之外源項
    
    a[n][n] =(6.0*Gamma) - (dx/2.0) - (dy/2.0);
    a[n][n-1] = -a_W;
    a[n][n-NX] = -a_S;
    b[n] = 2.0*a_N;// 邊界溫度
    
    // 下邊界 (除角點外)
    for(int i = 2; i <= NX-1; i++) {
        a[i][i] = (5*Gamma)+(dx/2.0);
        a[i][i-1] =-(Gamma + (dy/2.0));
        a[i][i+1] = -(Gamma - (dy/2.0));
        a[i][i+NX] = - (Gamma - (dx/2.0));
        b[i] = 0.0; // 邊界溫度
    }
    
    // 上邊界 (除角點外)
    for(int i = n-NX+2; i < n; i++) {
        a[i][i] = (5*Gamma)-(dx/2.0);
        a[i][i+1] = -(Gamma + (dy/2.0));
        a[i][i-1] = -(Gamma - (dy/2.0));
        a[i][i-NX] = -(Gamma + (dx/2.0));
        b[i] =2.0*Gamma- dx;
    }
    
    // 左邊界 (除角點外)
    for(int i = 1; i <= NY-2; i++) {
        int idx = NX*i+1;
        a[idx][idx] = (5*Gamma)+(dy/2.0);
        a[idx][idx+1] = -(-(dy/2.0) + Gamma) ;
        a[idx][idx-NX] = -((dx/2.0) + Gamma);
        a[idx][idx+NX] = -(-(dx/2.0) + Gamma);
        b[idx] = dy + 2.0 * Gamma ; //邊界溫度
    }
    
    // 右邊界 (除角點外)
    for(int i = 2; i <= NY-1; i++) {
        int idx = NX*i;
        a[idx][idx] = (5*Gamma)-(dy/2.0);
        a[idx][idx-1] = -((dy/2.0)+Gamma);
        a[idx][idx-NX] = -((dx/2.0)+Gamma);
        a[idx][idx+NX] = -(-(dx/2.0)+Gamma);
        b[idx] = 0.0;
    }
    
    // 內點
    for(int i = 2; i <= NX-1; i++) {
        for(int j = 1; j <= NY-2; j++) {
            int idx = j*NX+i;
            a[idx][idx] = (4*Gamma);
            a[idx][idx+1] = -(Gamma+(dy/2.0));
            a[idx][idx-1] = -(Gamma-(dy/2.0));
            a[idx][idx-NX] = -(Gamma+(dx/2.0));
            a[idx][idx+NX] = -(Gamma-(dx/2.0));
            b[idx] = 0.0; // 內點無熱源
        }
    }
}

/////////////////////////////////////////////////////////////////////////////
void coefficient(double a[][n+1] , double b[] , int n){
	//係數矩陣需要重新賦值
    //內點更新(i = 3 : NX-1 ; j = 3 : NY-1 ) 
    //左右內點X上下內點 
    for(int i = 3 ; i <= NX-1 ; i++){
    	for(int j = 3 ; j <= NY - 1 ; j++){
    		a[(j-1)*NX+i][(j-1)*NX+i] = 2*Gamma + dy*(1+B(rf_e(i,j))-A(rf_w(i,j))) + 2*Gamma + dx*(1+B(rf_n(i,j))-A(rf_s(i,j))) ;//本點
			a[(j-1)*NX+i][(j-1)*NX+i-1] = -Gamma + dy*(C(rf_e(i,j))-1-B(rf_w(i,j))) ;//W西鄰界計算點 
			a[(j-1)*NX+i][(j-1)*NX+i+1] = -Gamma + dy * A(rf_e(i,j));//E東鄰界計算點
			a[(j-1)*NX+i][(j-2)*NX+i] =-Gamma + dx*(C(rf_n(i,j))-1-B(rf_s(i,j)));//S南鄰界計算點
			a[(j-1)*NX+i][(j)*NX+i] = -Gamma + dx* A(rf_n(i,j));//N北鄰界計算點
			a[(j-1)*NX+i][(j-1)*NX+i-2] = -dy*C(rf_w(i,j));//WW西西鄰界計算點
			a[(j-1)*NX+i][(j-3)*NX+i] = -dx*C(rf_s(i,j));//SS南南鄰界計算點
			b[(j-1)*NX+i] = 0.0 ;//非齊性項 
			 
		}
	} 
    //六個邊界初始化(left\right\up\bottom\left2\bottom2) 
    //左邊界計算點
    //左右左邊界X上下內點 
	for(int j = 3 ; j <= NY-1 ; j++){//範圍:i = 1 ; j = 3 : NY-1 
		a[(j-1)*NX+1][(j-1)*NX+1] = (3*Gamma+dy*(1+B(rf_e(1,j))-C(rf_e(1,j))))+(2*Gamma + dx*(1+B(rf_n(1,j))-A(rf_s(1,j)))) ; //本點
		//a[(j-1)*NX+1][(j-1)*NX] = ? ;//不存在
		a[(j-1)*NX+1][(j-1)*NX+2] = (-Gamma+dy*A(rf_e(1,j))); //E東鄰界計算點 
		a[(j-1)*NX+1][(j-2)*NX+1] = -Gamma + dx*(C(rf_n(1,j))-1-B(rf_s(1,j)));//S南邊界計算點
		a[(j-1)*NX+1][(j)*NX+1] = -Gamma + dx* A(rf_n(1,j));//N北邊界計算點
		a[(j-1)*NX+1][(j-3)*NX+1] = -dx*C(rf_s(1,j)) ;//SS南南邊界計算點
		b[(j-1)*NX+1] = -(-Gamma+dy*(2*C(rf_e(1,j))-1))*T_Bl;
	} 
	//右邊界計算點
    //左右右邊界X上下內點
	for(int j = 3 ; j <= NY-1 ; j++){//範圍:i = NX ; j = 3 : NY-1 
		a[(j-1)*NX + NX][(j-1)*NX + NX] = (3*Gamma-dy*A(rf_w(NX,j)))+(2*Gamma + dx*(1+B(rf_n(NX,j))-A(rf_s(NX,j)))) ; //本點
		a[(j-1)*NX + NX][(j-1)*NX + NX-1]  = (-Gamma -dy*(B(rf_w(NX,j))+1));//西點
		a[(j-1)*NX + NX][(j-1)*NX] = -Gamma + dx*(C(rf_n(NX,j))-1-B(rf_s(NX,j)));//南點
		a[(j-1)*NX + NX][(j-1)*NX + 2*NX] = -Gamma + dx* A(rf_n(NX,j));//北點
		a[(j-1)*NX + NX][(j-1)*NX + NX-2] = -dy * C(rf_w(NX,j)); //西西點 
	    a[(j-1)*NX + NX][(j-1)*NX -NX] = -dx*C(rf_s(NX,j)) ;//南南點 
	    b[(j-1)*NX + NX] = 0.0 ;//非齊性項 
	} 
	//南邊界計算點
	//左右內點X上下下邊界 
	for(int i = 3 ; i <= NX-1 ; i++){//範圍:i = 3 : NX-1 ; j = 1 
		a[i][i] = (2*Gamma + dy*(1+B(rf_e(i,1)))-A(rf_w(i,1)))+(3*Gamma+dx*(1+B(rf_n(i,1))-C(rf_n(i,1)))) ;//本點
		a[i][i-1] = -Gamma + dy*(C(rf_e(i,1))-1-B(rf_w(i,1))) ;//西點
		a[i][i+1] = -Gamma + dy * A(rf_e(i,1));//東點
		a[i][i+NX] = (-Gamma+dx*A(rf_n(i,1)));//北點
		a[i][i-2] = -dy*C(rf_w(i,1)); //西西點
		b[i]  = 0.0 ; //非齊性項 
	}
	//北邊界計算點
	//左右內點X上下上邊界
	for(int i = 3 ;i <= NX-1 ; i++){//範圍:i = 3 : NX-1 ; j = NY 
		a[(NY-1)*NX+i][(NY-1)*NX+i] = (2*Gamma + dy*(1+B(rf_e(i,NY))-A(rf_w(i,NY))))+(3*Gamma-dx*A(rf_s(i,NY)));//本點
		a[(NY-1)*NX+i][(NY-1)*NX+i-1] =-Gamma + dy*(C(rf_e(i,NY))-1-B(rf_w(i,NY))) ;//西點
		a[(NY-1)*NX+i][(NY-1)*NX+i+1] =-Gamma + dy * A(rf_e(i,NY));//東點
		a[(NY-1)*NX+i][(NY-2)*NX+i] = (-Gamma-dx*(B(rf_s(i,NY))+1)); //南點
		a[(NY-1)*NX+i][(NY-1)*NX+i-2] = -dy*C(rf_w(i,NY)); //西西點
		a[(NY-1)*NX+i][(NY-3)*NX+i]  = -dx*(C(rf_s(i,NY))); //南南點
		b[(NY-1)*NX+i] = -(-2*Gamma+dx)*T_Bu;//非齊性項 
	} 
	//左邊界第二計算點
	//左右左第二排X上下內點
	for(int j = 3 ; j <= NY-1 ; j++){//範圍:i = 2 ; j = 3 : NY-1 
		a[(j-1)*NX+2][(j-1)*NX+2] = (2*Gamma+dy*(B(rf_e(2,j))+1-A(rf_w(2,j))))+(2*Gamma + dx*(1+B(rf_n(2,j))-A(rf_s(2,j)))); //本點
		a[(j-1)*NX+2][(j-1)*NX+1] = (-Gamma+dy*(C(rf_e(2,j))-1-B(rf_w(2,j))+C(rf_w(2,j)))); //西點
		a[(j-1)*NX+2][(j-1)*NX+3] = (-Gamma+dy*A(rf_e(2,j))); //東點
		a[(j-1)*NX+2][(j-2)*NX+2] = -Gamma + dx*(C(rf_n(2,j))-1-B(rf_s(2,j)));//南點
		a[(j-1)*NX+2][(j)*NX+2] = -Gamma + dx* A(rf_n(2,j)); //北點
		a[(j-1)*NX+2][(j-3)*NX+2] = -dx*C(rf_s(2,j)) ; //南南點
		b[(j-1)*NX+2] = 2*dy*T_Bl;//非齊性項 
	} 
	//下邊界第二計算點
	//左右內點X上下下第二排 
	for(int i = 3 ; i <= NX-1 ; i++){//範圍: i = 3 : NX-1 ; j = 2
		a[NX + i][NX + i] = (2*Gamma + dy*(1+B(rf_e(i,2))-A(rf_w(i,2)))) + (2*Gamma+dx*(B(rf_n(i,2))+1-A(rf_s(i,2))));//本點
		a[NX + i][NX + i-1] = -Gamma + dy*(C(rf_e(i,2))-1-B(rf_w(i,2))); //西點 
		a[NX + i][NX + i+1] =-Gamma + dy * A(rf_e(i,2)) ; //東點
		a[NX + i][i] = (-Gamma+dx*(C(rf_n(i,2))-1-B(rf_s(i,2))+C(rf_s(i,2)))); //南點
		a[NX + i][2*NX + i]  = (-Gamma+dx*A(rf_n(i,2))); //北點
		a[NX + i][NX + i-2]  =-dy*C(rf_w(i,2)) ; //西西點
		b[NX+i] = 0.0 ; //非齊性項 
	}
	//九個角點初始化
	//第一個角點
	//左右左邊界X上下下邊界
	//範圍:i = 1 ; j = 1 ;
	a[1][1]  = (3*Gamma+dy*(1+B(rf_e(1,1))-C(rf_e(1,1))))+(3*Gamma+dx*(1+B(rf_n(1,1))-C(rf_n(1,1))))  ; //本點 
	a[1][2] = (-Gamma+dy*A(rf_e(1,1)));//東點 
	a[1][1+NX]= (-Gamma+dx*A(rf_n(1,1)));//北點
	b[1] = -(-Gamma+dy*(2*C(rf_e(1,1))-1))*T_Bl; //非齊性項 
	//第二個角點
	//左右左第二排邊界X上下下邊界
	//範圍:i = 2 ; j = 1 ;
	a[2][2]  = (2*Gamma+dy*(B(rf_e(2,1))+1-A(rf_w(2,1))))+(3*Gamma+dx*(1+B(rf_n(2,1))-C(rf_n(2,1))))  ; //本點 
	a[2][1] = (-Gamma+dy*(C(rf_e(2,1))-1-B(rf_w(2,1))+C(rf_w(2,1)))); //西點 
	a[2][3] = (-Gamma+dy*A(rf_e(2,1)));//東點 
	a[2][2+NX] = (-Gamma+dx*A(rf_n(2,1)));//北點
	b[2] = 2*dy*T_Bl;//非齊性項
	//第三角點
	//左右右邊界X上下下邊界
	//範圍i = NX ; j = 1 ;
	a[NX][NX]  = (3*Gamma-dy*(A(rf_w(NX,1))))+(3*Gamma+dx*(1+B(rf_n(NX,1))-C(rf_n(NX,1))))  ; //本點 
	a[NX][NX-1] = (-Gamma-dy*(B(rf_w(NX,1))+1));//西點
	a[NX][2*NX] = (-Gamma+dx*A(rf_n(NX,1)));//北點
	a[NX][NX-2] = (-dy*C(rf_w(NX,1))); //西西點
	b[NX] = 0.0 ;//非齊性項  
	//第四角點
	//左右左邊界X上下下第二排
	//範圍:i = 1 ; j = 2 ;
	a[NX+1][NX+1]  = (3*Gamma+dy*(1+B(rf_e(1,2))-C(rf_e(1,2)))) + (2*Gamma+dx*(B(rf_n(1,2))+1-A(rf_s(1,2))));//本點 
	a[NX+1][NX+2] = (-Gamma+dy*A(rf_e(1,2))) ;//東點 
	a[NX+1][1] = (-Gamma+dx*(C(rf_n(1,2))-1-B(rf_s(1,2))+C(rf_s(1,2))));//南點 
	a[NX+1][2*NX+1] = (-Gamma+dx*A(rf_n(1,2)));//北點 
	b[NX+1] = -(-Gamma+dy*(2*C(rf_e(1,2))-1))*T_Bl; //非齊性項
	//第五角點
	//左右左第二排X上下下第二排
	//範圍:i = 2 ; j = 2 ;
	a[NX+2][NX+2] = (2*Gamma+dy*(B(rf_e(2,2))+1-A(rf_w(2,2)))) + (2*Gamma+dx*(B(rf_n(2,2))+1-A(rf_s(2,2))));//本點  
	a[NX+2][NX+1] = (-Gamma+dy*(C(rf_e(2,2))-1-B(rf_w(2,2))+C(rf_w(2,2)))); //西點 
	a[NX+2][NX+3] = (-Gamma+dy*A(rf_e(2,2)));//東點
	a[NX+2][2] = (-Gamma+dx*(C(rf_n(2,2))-1-B(rf_s(2,2))+C(rf_s(2,2))));//南點 
	a[NX+2][2*NX+2] = (-Gamma+dx*A(rf_n(2,2)));//北點
	b[NX+2] = 2*dy*T_Bl;//非齊性項 
	//第六角點
	//左右右邊界X上下下第二排
	//範圍:i = NX ; j = 2 ;
	a[2*NX][2*NX] = (3*Gamma-dy*(A(rf_w(NX,2))))+ (2*Gamma+dx*(B(rf_n(NX,2))+1-A(rf_s(NX,2))));//本點
	a[2*NX][2*NX-1] = (-Gamma-dy*(B(rf_w(NX,2))+1));//西點
	a[2*NX][NX] = (-Gamma+dx*(C(rf_n(NX,2))-1-B(rf_s(NX,2))+C(rf_s(NX,2))));//南點 
	a[2*NX][3*NX] = (-Gamma+dx*A(rf_n(NX,2)));//北點
	a[2*NX][2*NX-2] = (-dy*C(rf_w(NX,2))); //西西點
	b[2*NX] = 0.0 ;//非齊性項 
    //第七角點
	//左右左邊界X上下上邊界
	//範圍:i = 1 ; j = NY ; 
	a[NX*(NY-1)+1][NX*(NY-1)+1] = (3*Gamma+dy*(1+B(rf_e(1,NY))-C(rf_e(1,NY))))+(3*Gamma-dx*A(rf_s(1,NY))) ;//本點
	a[NX*(NY-1)+1][NX*(NY-1)+2] = (-Gamma+dy*A(rf_e(1,NY))) ;//東點 
	a[NX*(NY-1)+1][NX*(NY-2)+1] = (-Gamma-dx*(B(rf_s(1,NY))+1)); //南點 
	a[NX*(NY-1)+1][NX*(NY-3)+1] = -dx*C(rf_s(1,NY)); //南南點 
	b[NX*(NY-1)+1] = -(-Gamma+dy*(2*C(rf_e(1,NY))-1))*T_Bl-(-2*Gamma+dx)*T_Bu; //非齊性項
	//第八角點
	//左右左第二排X上下上邊界
	//範圍:i = 2 ; j = NY ;
	a[NX*(NY-1)+2][NX*(NY-1)+2] = (2*Gamma+dy*(B(rf_e(2,NY))+1-A(rf_w(2,NY))))+(3*Gamma-dx*A(rf_s(2,NY))) ;//本點
	a[NX*(NY-1)+2][NX*(NY-1)+1] = (-Gamma+dy*(C(rf_e(2,NY))-1-B(rf_w(2,NY))+C(rf_w(2,NY)))); //西點  
	a[NX*(NY-1)+2][NX*(NY-1)+3] = (-Gamma+dy*A(rf_e(2,NY)));//東點 
	a[NX*(NY-1)+2][NX*(NY-2)+2] = (-Gamma-dx*(B(rf_s(2,NY))+1)); //南點  
	a[NX*(NY-1)+2][NX*(NY-3)+2] = -dx*C(rf_s(2,NY)); //南南點  
	b[NX*(NY-1)+2] = 2*dy*T_Bl-(-2*Gamma+dx)*T_Bu; //非齊性項
	//第九角點
	//左右右邊界X上下上邊界
	//範圍:i = NX ; j = NY ;
	a[NX*NY][NX*NY] = (3*Gamma-dy*(A(rf_w(NX,NY))))+(3*Gamma-dx*A(rf_s(NX,NY))) ;//本點
	a[NX*NY][NX*NY-1] = (-Gamma-dy*(B(rf_w(NX,NY))+1));//西點 
	a[NX*NY][NX*(NY-1)] = (-Gamma-dx*(B(rf_s(NX,NY))+1)); //南點  
	a[NX*NY][NX*NY-2] = (-dy*C(rf_w(NX,NY))); //西西點
	a[NX*NY][NX*(NY-2)] = -dx*C(rf_s(NX,NY)); //南南點  
	b[NX*NY] = 0.0 -(-2*Gamma+dx)*T_Bu; //非齊性項
	
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Jacobi(double a[][n+1], double b[], double x[], int n) {
    // 先複製當前解到x_old
    for(int k = 1; k <= n; k++) {
        x_old[k] = x[k];
    }
    
    // 計算新的解
    for(int k = 1; k <= n; k++) {
        double sum = 0;
        for(int p = 1; p <= n; p++) {
            if(p != k) {
                sum += a[k][p] * x_old[p];
            }
        }
        x[k] = (b[k] - sum) / a[k][k];
    }
    
    // 計算最大誤差
    maxerror = 0;
    for(int k = 1; k <= n; k++) {
        double error = fabs(x[k] - x_old[k]);
        if(maxerror < error) {
            maxerror = error;
        }
    }
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
    name << "bounded_scheme_3_Pe = " << Pelect <<","<<"Grid "<< nodes_i << "x" << nodes_j << ","<< setfill('0') << setw(6) << m << ".vtk";
    ofstream out(name.str().c_str());
    
    // VTK 文件頭：使用 STRUCTURED_GRID
    out << "# vtk DataFile Version 3.0\n";
    out << "bounded_scheme_3\n";
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
    initial(a, b, n);//首先\，溫度，係數矩陣，非齊性項先進行初始化 
    orderCD(a,b,n) ;//初始狀態先以 
   for(G = 0; G < max_G; G++) {
        if(G % 1000== 0) {
            cout << "迭代次數 = " << G << endl;
            if(G > 0) {
                cout << "最大變化量max_change = " << scientific << maxerror << endl;
            }
            output(G);
        }
        
        if(G > 0 && maxerror < tolerance) {  // 修正：確保第一次迭代後才檢查收斂
            steadystate = true;
            cout << "已經達到迭代終點，溫度場收斂!!" << endl;
            break;
        }
        
        Jacobi(a,b,x,n) ;//求解溫度計算誤差 
        coefficient(a,b,n) ;//更新溫度利馬更新係數 
    }
    
    if(!steadystate) {
        cout << "達到最大迭代次數，但未達到穩態!" << endl;
    }
    
    output(G);
    cout << "格子大小 " << NX << "x" << NY << " 計算完成\n" << endl;
    return 0;
}

///////////////////////////////////////////////////////////////////////////////////