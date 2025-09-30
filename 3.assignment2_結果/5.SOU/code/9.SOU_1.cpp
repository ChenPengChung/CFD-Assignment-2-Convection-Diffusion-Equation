//利用空間二階精度upwind格式(SOU)求解熱擴散對流問題 
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
//會用到的參數
int Nx[] = {40,80};
//有限體積法採用Wet node處理，故計算點數目 = 網格數目, 節點數目 = 網格數目+1
int NX; //NX = NY => dx = dy
int NY = NX ;
int n;
vector<vector<double> > a; //a二維矩陣為等式左側係數矩陣
vector<double> b; //b一維矩陣為等式右側之外源項，方程組之非齊性項
vector<double> x;
vector<double> x_old;
vector<double> x_old2;
vector<double> x_nuse ;
vector<vector<double> > T;
const double Pelect = 1.0 ;
const double  Gamma = (sqrt(2)/Pelect) ;//熱擴散係數 alpha
double dx ,dy  ;  // 修正：網格間距應為1/(NX-1)
const double T_Bl = 1.0 ;
const double T_Bb = 0.0 ;
double lamda = 0.8 ;
int G, max_G = 10000000000;
double maxerror; 
const float tolerance = 1e-10;
bool steadystate;


//初始化矩陣
void initial(vector<vector<double> >& a, vector<double>& b , vector<double>& x) {
    // 初始化所有元素為0
    for(int i = 1; i <= n; i++) {
        for(int j = 1; j <= n; j++) {
            a[i][j] = 0.0;
        }
        b[i] = 0.0;
        x[i] = 0.0;
        x_old[i] = 0.0;
    }
    a[0][0] = 0.0 ;
    a[n+1][n+1] = 0.0 ;
    //係數矩陣賦值 
    //內點更新(i = 3 : NX-1 ; j = 3 : NY-1 ) 
    //左右內點X上下內點 
    for(int i = 3 ; i <= NX-1 ; i++){
    	for(int j = 3 ; j <= NY - 1 ; j++){
    		a[(j-1)*NX+i][(j-1)*NX+i]   =  4*Gamma + (3.0/2.0)*dx+ (3.0/2.0)*dx ;//本點
			a[(j-1)*NX+i][(j-1)*NX+i-1] =  -(Gamma+2.0*dx);//W西鄰界計算點 
			a[(j-1)*NX+i][(j-1)*NX+i+1] =  -(Gamma);//E東鄰界計算點
			a[(j-1)*NX+i][(j-2)*NX+i]   =  -(Gamma+2.0*dx);//S南鄰界計算點
			a[(j-1)*NX+i][(j)*NX+i]     =  -(Gamma);//N北鄰界計算點
			a[(j-1)*NX+i][(j-1)*NX+i-2] = -(-(1.0/2.0)*dx);//WW西西鄰界計算點
			a[(j-1)*NX+i][(j-3)*NX+i]   = -(-(1.0/2.0)*dx);//SS南南鄰界計算點
			b[(j-1)*NX+i] = 0.0 ;//非齊性項 
			 
		}
	} 
    //六個邊界初始化(left\right\up\bottom\left2\bottom2) 
    //左邊界計算點
	for(int j = 3 ; j <= NY-1 ; j++){//範圍:i = 1 ; j = 3 : NY-1 
		a[(j-1)*NX+1][(j-1)*NX+1] = 5*Gamma +2*dx+(3.0/2.0)*dx; //本點
		//a[(j-1)*NX+1][(j-1)*NX] = ? ;//不存在
		a[(j-1)*NX+1][(j-1)*NX+2] = -Gamma ; //E東鄰界計算點 
		a[(j-1)*NX+1][(j-2)*NX+1] = -(Gamma+2*dx);//S南邊界計算點
		a[(j-1)*NX+1][(j)*NX+1]   = -Gamma;//N北邊界計算點
		a[(j-1)*NX+1][(j-3)*NX+1] = (1.0/2.0)*dx;//SS南南邊界計算點
		b[(j-1)*NX+1]             = T_Bl*(2*Gamma+2.0*dx ) ;
	} 
	//右邊界計算點
    //左右右邊界X上下內點
	for(int j = 3 ; j <= NY-1 ; j++){//範圍:i = NX ; j = 3 : NY-1 
		a[(j-1)*NX + NX][(j-1)*NX + NX]    = 3*Gamma + dx + (3.0/2.0)*dx ; //本點
		a[(j-1)*NX + NX][(j-1)*NX + NX-1]  = -Gamma - dx * (3.0/2.0); //西點
		a[(j-1)*NX + NX][(j-1)*NX]         = -Gamma - dx * 2.0 ; //南點
		a[(j-1)*NX + NX][(j-1)*NX + 2*NX]  = -Gamma ; //北點
		a[(j-1)*NX + NX][(j-1)*NX + NX-2]  = dx*(1.0/2.0) ; //西西點 
	    a[(j-1)*NX + NX][(j-1)*NX -NX]     = dx*(1.0/2.0) ; //南南點 
	    b[(j-1)*NX + NX] = 0.0 ; //非齊性項 
	} 
	//下邊界計算點
	//左右內點X上下下邊界 
	for(int i = 3 ; i <= NX-1 ; i++){//範圍:i = 3 : NX-1 ; j = 1 
		a[i][i] = 5*Gamma + (3.0/2.0)*dx + 2*dx;//本點
		a[i][i-1] = -Gamma -2*dx; //西點
		a[i][i+1] = -Gamma; //東點
		a[i][i+NX]= -Gamma;//北點
		a[i][i-2] = (1.0/2.0)*dx; //西西點
		b[i]  =  T_Bb*(2*Gamma +2*dx); //非齊性項 
	}
	//上邊界計算點
	//左右內點X上下上邊界
	for(int i = 3 ;i <= NX-1 ; i++){//範圍:i = 3 : NX-1 ; j = NY 
		a[(NY-1)*NX+i][(NY-1)*NX+i]   = 3*Gamma + dx*(3.0/2.0) + dx ; //本點
		a[(NY-1)*NX+i][(NY-1)*NX+i-1] =  -Gamma -2*dx; //西點
		a[(NY-1)*NX+i][(NY-1)*NX+i+1] =  -Gamma; //東點
		a[(NY-1)*NX+i][(NY-2)*NX+i]   =  -Gamma -dx*(3.0/2.0); //南點
		a[(NY-1)*NX+i][(NY-1)*NX+i-2] = (1.0/2.0)*dx ; //西西點
		a[(NY-1)*NX+i][(NY-3)*NX+i]   = (1.0/2.0)*dx ; //南南點
		b[(NY-1)*NX+i]                = 0.0  ; //非齊性項 
	} 
	//左邊界第二計算點
	//左右左第二排X上下內點
	for(int j = 3 ; j <= NY-1 ; j++){//範圍:i = 2 ; j = 3 : NY-1 
		a[(j-1)*NX+2][(j-1)*NX+2] = 4*Gamma +dx*(3.0/2.0)+dx*(3.0/2.0) ; //本點
		a[(j-1)*NX+2][(j-1)*NX+1] = -Gamma -dx*(5.0/2.0) ; //西點
		a[(j-1)*NX+2][(j-1)*NX+3] = -Gamma ; //東點
		a[(j-1)*NX+2][(j-2)*NX+2] = -Gamma -2*dx;//南點
		a[(j-1)*NX+2][(j)*NX+2] =   -Gamma ; //北點
		a[(j-1)*NX+2][(j-3)*NX+2] = (1.0/2.0)*dx; //南南點
		b[(j-1)*NX+2]             = T_Bl * (-dx);//非齊性項 
	} 
	//下邊界第二計算點
	//左右內點X上下下第二排 
	for(int i = 3 ; i <= NX-1 ; i++){//範圍: i = 3 : NX-1 ; j = 2
		a[NX + i][NX + i]   = 4*Gamma + (3.0/2.0)*dx + (3.0/2.0)*dx; //本點
		a[NX + i][NX + i-1] = -Gamma -2*dx ; //西點 
		a[NX + i][NX + i+1] = -Gamma ; //東點
		a[NX + i][i]        = -Gamma -(5.0/2.0)*dx ; //南點
		a[NX + i][2*NX + i] = -Gamma ; //北點
		a[NX + i][NX + i-2] = (1.0/2.0)*dx ; //西西點
		b[NX+i]             = T_Bb * (-dx) ; //非齊性項 
	}
	//九個角點初始化
	//第一個角點
	//左右左邊界X上下下邊界
	//範圍:i = 1 ; j = 1 ;
	a[1][1]   = 6*Gamma + 4*dx; //本點 
	a[1][2]   = -Gamma ; //東點 
	a[1][1+NX]= -Gamma ; //北點
	b[1]      = T_Bl*(2.0*Gamma+2.0*dx) + T_Bb*(2.0*Gamma+2.0*dx) ; //非齊性項 
	//第二個角點
	//左右左第二排邊界X上下下邊界
	//範圍:i = 2 ; j = 1 ;
	a[2][2]    = 5*Gamma + (3.0/2.0)*dx + 2*dx ; //本點 
	a[2][1]    = -Gamma-(5.0/2.0)*dx ;//西點 
	a[2][3]    = -Gamma; //東點 
	a[2][2+NX] = -Gamma; //北點
	b[2]       =  T_Bb*(2*Gamma-dx+2*dx ); //非齊性項
	//第三角點
	//左右右邊界X上下下邊界
	//範圍i = NX ; j = 1 ;
	a[NX][NX]   = 4*Gamma + dx + 2*dx; //本點 
	a[NX][NX-1] = -Gamma -(3.0/2.0)*dx ; //西點
	a[NX][2*NX] = -Gamma ; //北點
	a[NX][NX-2] = dx*(1.0/2.0); //西西點
	b[NX]       = T_Bb*(2*Gamma+dx+dx); //非齊性項  
	//第四角點
	//左右左邊界X上下下第二排
	//範圍:i = 1 ; j = 2 ;
	a[NX+1][NX+1]  = 5*Gamma + 2.0*dx + (3.0/2.0)*dx ; //本點 
	a[NX+1][NX+2]  = -Gamma ; //東點 
	a[NX+1][1]     = -Gamma -(5.0/2.0)*dx ; //南點 
	a[NX+1][2*NX+1]= -Gamma ; //北點 
	b[NX+1]        = T_Bl*(2*Gamma+2.0*dx-dx ); //非齊性項
	//第五角點
	//左右左第二排X上下下第二排
	//範圍:i = 2 ; j = 2 ;
	a[NX+2][NX+2]   = 4*Gamma + (3.0/2.0)*dx + (3.0/2.0)*dx ; //本點  
	a[NX+2][NX+1]   = -Gamma - (5.0/2.0)*dx ; //西點 
	a[NX+2][NX+3]   = -Gamma ; //東點
	a[NX+2][2]      = -Gamma - (5.0/2.0)*dx ; //南點 
	a[NX+2][2*NX+2] = -Gamma ; //北點
	b[NX+2]         = -dx*T_Bl -dx*T_Bb ; //非齊性項 
	//第六角點
	//左右右邊界X上下下第二排
	//範圍:i = NX ; j = 2 ;
	a[2*NX][2*NX]   = 3*Gamma + (5.0/2.0)*dx ; //本點
	a[2*NX][2*NX-1] = -Gamma-(3.0/2.0)*dx ; //西點
	a[2*NX][NX]     = -Gamma-(5.0/2.0)*dx ; //南點 
	a[2*NX][3*NX]   = -Gamma; //北點
	a[2*NX][2*NX-2] = (1.0/2.0)*dx ; //西西點
	b[2*NX]         = -dx*T_Bb ; //非齊性項 
    //第七角點
	//左右左邊界X上下上邊界
	//範圍:i = 1 ; j = NY ; 
	a[NX*(NY-1)+1][NX*(NY-1)+1] = 4*Gamma + 3.0*dx ; //本點
	a[NX*(NY-1)+1][NX*(NY-1)+2] = -Gamma ; //東點 
	a[NX*(NY-1)+1][NX*(NY-2)+1] = -Gamma-(3.0/2.0)*dx ; //南點 
	a[NX*(NY-1)+1][NX*(NY-3)+1] = (1.0/2.0)*dx ; //南南點 
	b[NX*(NY-1)+1]              = T_Bl*(2.0*dx+2*Gamma ); //非齊性項
	//第八角點
	//左右左第二排X上下上邊界
	//範圍:i = 2 ; j = NY ;
	a[NX*(NY-1)+2][NX*(NY-1)+2] = 3*Gamma + dx + (3.0/2.0)*dx ; //本點
	a[NX*(NY-1)+2][NX*(NY-1)+1] = -Gamma-(5.0/2.0)*dx; //西點  
	a[NX*(NY-1)+2][NX*(NY-1)+3] = -Gamma; //東點 
	a[NX*(NY-1)+2][NX*(NY-2)+2] = -Gamma - (3.0/2.0)*dx ; //南點  
	a[NX*(NY-1)+2][NX*(NY-3)+2] = (1.0/2.0)*dx ; //南南點 
	b[NX*(NY-1)+2]              = -dx*T_Bl ; //非齊性項
	//第九角點
	//左右右邊界X上下上邊界
	//範圍:i = NX ; j = NY ;
	a[NX*NY][NX*NY]     =  2*Gamma+2*dx; //本點
	a[NX*NY][NX*NY-1]   =  -Gamma-(3.0/2.0)*dx ; //西點 
	a[NX*NY][NX*(NY-1)] =  -Gamma-(3.0/2.0)*dx ; //南點  
	a[NX*NY][NX*NY-2]   =  (1.0/2.0)*dx ; //西西點
	a[NX*NY][NX*(NY-2)] =  (1.0/2.0)*dx ; //南南點  
	b[NX*NY]            =  0.0 ; //非齊性項
	
}

void Jacobi(vector<vector<double> >& a, vector<double>& b , vector<double>& x , int n) {
	if(NX ==  40){
        lamda = 0.05 ;
    }else{
        lamda = 0.08 ;
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
    name << "9.SOU_1Pe = " << Pelect <<","<<"Grid "<< nodes_i << "x" << nodes_j << ","<< setfill('0') << setw(6) << m << ".vtk";
    ofstream out(name.str().c_str());
    
    // VTK 文件頭：使用 STRUCTURED_GRID
    out << "# vtk DataFile Version 3.0\n";
    out << "9.SOU_1Pe = " << Pelect <<","<<"Grid "<< nodes_i << "x" << nodes_j << ","<< setfill('0') << setw(6) << m << ".vtk\n";
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
		NY = NX ;
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
        cout << "Boundary condition : up&right : Neumann n dot Gradient of T = 0 " << endl;
        cout << "Boundary condition : left = 1.0,  and bottom = 0.0" << endl;


        steadystate = false;
        initial(a, b, x);
    
        for(G = 0; G < max_G; G++) {
            if(G > 0) {
                Jacobi(a, b, x, n);
            }
        if(G % 500 == 0) {
            cout << "itteration = " << G << endl;
            if(G > 0) {
                cout << "the max_change = " << scientific << maxerror << endl;
            }
            output(G);
        }
        
        if(G > 100 && maxerror < tolerance) {
            steadystate = true;
            cout << "will be steady state and temperature field converged!!" << endl;
            break;
        }
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