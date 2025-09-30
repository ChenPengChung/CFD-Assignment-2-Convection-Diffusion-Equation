//利用TVD求解二維穩態熱擴散對流方程
//引入鬆弛技術
//增加對角係數 
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <iomanip>
#include <vector>
#define W 1e-10 

using namespace std;
//會用到的參數
//會用到的參數
int Nx[] = {40,80};
//有限體積法採用Wet node處理，故計算點數目 = 網格數目, 節點數目 = 網格數目+1
int NX , NY; //NX = NY => dx = dy
int n;
vector<vector<double> > a; //a二維矩陣為等式左側係數矩陣
vector<double> b; //b一維矩陣為等式右側之外源項，方程組之非齊性項
vector<double> x;
vector<double> x_old;
vector<double> x_old2;
vector<double> x_nuse ;
vector<vector<double> > T;
const double Pelect = 10.0 ;
const double  Gamma = (1.0/Pelect) ;//熱擴散係數 alpha
double dx ,dy  ;  // 修正：網格間距應為1/(NX-1)
const double T_Bl = 1.0 ;//邊界溫度left  = 1.0 ;
const double T_Bb = 0.0 ;//邊界溫度bottom  0.0 ;
const double u = 1.0 ;
const double v = 1.0 ;
double lamda ; //使計算結果顯示欠鬆弛under relaxtion 
double S_u1l = 2*Gamma+dy*u ;
double S_u1r = 2*Gamma-dy*u ;
double S_u2b = 2*Gamma+dx*v ;
double S_u2u = 2*Gamma-dx*v ;
int G, max_G = 100000000;
double maxerror ;
const double tolerance = 1e-10;
bool steadystate; 
/////////////////////////////////////////////////////////////////////////////////////////
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
		return 2*(x[(j-2)*NX+i]-T_Bb )/((x[(j-1)*NX+i] - x[(j-2)*NX+i])+W) ;
	}else{
		return (x[(j-2)*NX+i] - x[(j-3)*NX+i])/((x[(j-1)*NX+i] - x[(j-2)*NX+i])+W) ;
	}//j = 1 ; 用不到 
} 
//北邊界通量限制器函數
double rf_n(int i , int j){//只討論i = 1 : NX , j = 1 : NY-1 //當j=1時，rf_n需要特殊討論 
	if(j == 1){
		return 2*(x[(j-1)*NX+i]-T_Bb )/((x[(j)*NX+i] - x[(j-1)*NX+i])+W) ;
	}else{
		return (x[(j-1)*NX+i] - x[(j-2)*NX+i])/((x[(j)*NX+i] - x[(j-1)*NX+i])+W) ;
	}//j = NY ; 用不到 
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
double psi(double p){
	return max(0.0, min(2.0*p, min((p+1.0)/2.0, 2.0))) ;  // 修正多重比較問題//muscl格式各個點各個邊界的通量限制器函數 
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double a_W(double r){
	return Gamma + dy*u*(1-(r/2.0)) ;
} 
double a_E(double r){
	return Gamma - dy*u*((r)/2.0);
}
double a_S(double r){
	return Gamma + dx*v*(1-(r/2.0)) ;
} 
double a_N(double r){
	return Gamma - dx*v*((r)/2.0) ;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void initial(vector<vector<double> >& a, vector<double>& b , vector<double>& x){
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
/////////////////////////////////////////////////////////////////////////////
void coefficient(vector<vector<double> >& a, vector<double>& b , vector<double>& x){
	//係數矩陣需要重新賦值
    //內點更新(i = 3 : NX-1 ; j = 3 : NY-1 ) 
    //左右內點X上下內點 
    for(int i = 2 ; i <= NX-1 ; i++){
    	for(int j = 2 ; j <= NY - 1 ; j++){
    		a[(j-1)*NX+i][(j-1)*NX+i] = (a_W(psi(rf_w(i,j)))+a_E(psi(rf_e(i,j))))+(a_S(psi(rf_s(i,j)))+a_N(psi(rf_n(i,j)))) ;//本點
			a[(j-1)*NX+i][(j-1)*NX+i-1] = -a_W(psi(rf_w(i,j))) ; //W西鄰界計算點 
			a[(j-1)*NX+i][(j-1)*NX+i+1] = -a_E(psi(rf_e(i,j))) ; //E東鄰界計算點
			a[(j-1)*NX+i][(j-2)*NX+i]   = -a_S(psi(rf_s(i,j))) ; //S南鄰界計算點
			a[(j-1)*NX+i][(j)*NX+i]     = -a_N(psi(rf_n(i,j))) ; //N北鄰界計算點
			b[(j-1)*NX+i] = 0.0 ;//非齊性項 
			 
		}
	} 
    //四個邊界初始化(left\right\up\bottom) 
    //左邊界計算點
    //左右左邊界X上下內點 
	for(int j = 2 ; j <= NY-1 ; j++){//範圍:i = 1 ; j = 2 : NY-1 
		a[(j-1)*NX+1][(j-1)*NX+1] =  (a_E(psi(rf_e(1,j))) + S_u1l) + (a_S(psi(rf_s(1,j))) +a_N(psi(rf_n(1,j)))); //本點
		a[(j-1)*NX+1][(j-1)*NX+2] =  -a_E(psi(rf_e(1,j))) ; //E東鄰界計算點 
		a[(j-1)*NX+1][(j-2)*NX+1] =  -a_S(psi(rf_s(1,j))) ; //S南邊界計算點
		a[(j-1)*NX+1][(j)*NX+1]   =  -a_N(psi(rf_n(1,j))) ; //N北邊界計算點
		b[(j-1)*NX+1] = S_u1l*T_Bl ;
	} 
	//右邊界計算點
    //左右右邊界X上下內點
	for(int j = 2 ; j <= NY-1 ; j++){//範圍:i = NX ; j = 2 : NY-1 
		a[(j-1)*NX + NX][(j-1)*NX + NX]    = (a_W(psi(rf_w(NX,j))) - dx + dx ) + (a_S(psi(rf_s(NX,j))) + a_N(psi(rf_n(NX,j)))); //本點
		a[(j-1)*NX + NX][(j-1)*NX + NX-1]  = -a_W(psi(rf_w(NX,j))) ; //西點
		a[(j-1)*NX + NX][(j-1)*NX]         = -a_S(psi(rf_s(NX,j))) ; //南點
		a[(j-1)*NX + NX][(j-1)*NX + 2*NX]  = -a_N(psi(rf_n(NX,j))) ; //北點
	    b[(j-1)*NX + NX]                   = 0.0  ;//非齊性項 
	} 
	//下邊界計算點
	//左右內點X上下下邊界 
	for(int i = 2 ; i <= NX-1 ; i++){//範圍:i = 2: NX-1 ; j = 1 
		a[i][i]    = (a_W(psi(rf_w(i,1))) + a_E(psi(rf_e(i,1))))+(a_N(psi(rf_n(i,1))) + S_u2b); //本點
		a[i][i-1]  = -a_W(psi(rf_w(i,1))) ; //西點
		a[i][i+1]  = -a_E(psi(rf_e(i,1))) ; //東點
		a[i][i+NX] = -a_N(psi(rf_n(i,1))) ; //北點
		b[i]       = S_u2b*T_Bb; //非齊性項 
	}
	//上邊界計算點
	//左右內點X上下上邊界
	for(int i = 2 ;i <= NX-1 ; i++){//範圍:i = 2 : NX-1 ; j = NY 
		a[(NY-1)*NX+i][(NY-1)*NX+i]   = (a_W(psi(rf_w(i,NY))) + a_E(psi(rf_e(i,NY))))+(a_S(psi(rf_s(i,NY))) -dx + dx ); //本點
		a[(NY-1)*NX+i][(NY-1)*NX+i-1] = -a_W(psi(rf_w(i,NY))) ; //西點
		a[(NY-1)*NX+i][(NY-1)*NX+i+1] = -a_E(psi(rf_e(i,NY))) ; //東點
		a[(NY-1)*NX+i][(NY-2)*NX+i]   = -a_S(psi(rf_s(i,NY))) ; //南點
		b[(NY-1)*NX+i]                =  0.0            ; //非齊性項 
	} 
	//四個角點初始化
	//第一個角點
	//左右左邊界X上下下邊界
	//範圍:i = 1 ; j = 1 ;
	a[1][1]    = (a_E(psi(rf_e(1,1))) + S_u1l) + (a_N(psi(rf_n(1,1))) + S_u2b ); //本點 
	a[1][2]    = -a_E(psi(rf_e(1,1))) ; //東點 
	a[1][1+NX] = -a_N(psi(rf_n(1,1)))  ; //北點
	b[1]       = (T_Bl*S_u1l) + (T_Bb*S_u2b); //非齊性項 
	
	//第二角點
	//左右右邊界X上下下邊界
	//範圍i = NX ; j = 1 ;
	a[NX][NX]   = (a_W(psi(rf_w(NX,1))) - dx + dx ) + (a_N(psi(rf_n(NX,1))) + S_u2b ); //本點 
	a[NX][NX-1] = -a_W(psi(rf_w(NX,1))) ; //西點
	a[NX][2*NX] = -a_N(psi(rf_n(NX,1))) ; //北點
	b[NX]       =  (T_Bb*S_u2b); //非齊性項  
	
    //第三角點
	//左右左邊界X上下上邊界
	//範圍:i = 1 ; j = NY ; 
	a[NX*(NY-1)+1][NX*(NY-1)+1] =  (a_E(psi(rf_e(1,NY)))+S_u1l) + (a_S(psi(rf_s(1,NY)))  - dx + dx ) ; //本點
	a[NX*(NY-1)+1][NX*(NY-1)+2] = -a_E(psi(rf_e(1,NY))) ; //東點 
	a[NX*(NY-1)+1][NX*(NY-2)+1] = -a_S(psi(rf_s(1,NY))) ; //南點 
	b[NX*(NY-1)+1] = (T_Bl*S_u1l)+0.0 ; //非齊性項

	//第四角點
	//左右右邊界X上下上邊界
	//範圍:i = NX ; j = NY ;
	a[NX*NY][NX*NY]     = (a_W(psi(rf_w(NX,NY)))  -dx + dx ) + (a_S(psi(rf_s(NX,NY))) -dx + dx  ) ;//本點
	a[NX*NY][NX*NY-1]   = -a_W(psi(rf_w(NX,NY))) ; //西點 
	a[NX*NY][NX*(NY-1)] = -a_S(psi(rf_s(NX,NY))) ; //南點   
	b[NX*NY]            = 0.0  ; //非齊性項

}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Jacobi(vector<vector<double> >& a, vector<double>& b , vector<double>& x , int n) {
	if(NX ==  40){
        lamda = 0.5 ;
    }else{
        lamda = 0.8 ;
    } //自適應鬆弛因子
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
    name << "10.TVD-MUSCLprimescheme,Pe = " << Pelect <<","<<"Grid "<< nodes_i << "x" << nodes_j << ","<< setfill('0') << setw(6) << m << ".vtk";
    ofstream out(name.str().c_str());
    
    // VTK 文件頭：使用 STRUCTURED_GRID
    out << "# vtk DataFile Version 3.0\n";
    out << "10.TVD-MUSCLprimescheme,Pe = " << Pelect <<","<<"Grid "<< nodes_i << "x" << nodes_j << ","<< setfill('0') << setw(6) << m << ".vtk\n";
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
        NY=NX ;
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
        cout << "Boundary condition : up & right: Neumann n dot Gradient of T = 0 " << endl;
        cout << "Boundary condition : left = 1.0,  and bottom = 0.0" << endl;


        steadystate = false;
        initial(a, b, x);
        coefficient(a,b,x) ;//根據初始溫度進行賦值係數 
        for(G = 0; G < max_G; G++) {
        if(G % 1000== 0) {
            cout << "itteration number = " << G << endl;
            if(G > 0) {
                cout << "max_change = " << scientific << maxerror << endl;
            }
            output(G);
        }
        
        if(G > 0 && maxerror < tolerance) {  // 修正：確保第一次迭代後才檢查收斂
            steadystate = true;
            cout << "Have already converged!!" << endl;
            break;
        }
        
        Jacobi(a,b,x,n) ;//跟新溫度解併計算迭代誤差 (計算 G==1下的迭代誤差) 
        coefficient(a,b,x) ;//更新溫度利馬更新係數 
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


