//利用QUICK格式求解二維穩態純對流方程式
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
vector<vector<double> > A; //a二維矩陣為等式左側係數矩陣
vector<double> B; //b一維矩陣為等式右側之外源項，方程組之非齊性項
vector<double> X;
vector<double> X_old;
double dx; //一個網格間距
//邊界條件定義
const double T_left = 1.0;    //左邊界Dirichlet條件
const double T_right = 0.0;   //右邊界Dirichlet條件  
const double T_bottom = 0.0;  //下邊界Dirichlet條件
//上邊界為Neumann條件: ?T/?n = 0

//QUICK格式係數
double a1,a2,a3,a4;

//迭代參數
const double lamda = 0.3; //逐次超鬆弛迭代法
int G, max_G = 10000000;
double maxerror , maxerror2 ; 
const float tolerance = 1e-10;
bool Temp , steadystate;


//為了應對 Neumann條件，先求解最上層邊界作為一維溫度場的求解
//求解完之後，上層邊界的溫度場即作為已知的溫度再來迭代求解下面的點 
//上層邊界點的計算 執行一次就好 
//利用最上層邊界計算點作為一維場點求解溫度場
void cofficient(vector<vector<double> >& A, vector<double>& B , vector<double>& X){
	// 初始化所有元素為0
    for(int i = 1; i <= NX ; i++) {
            A[i][i] = 0.0 ;
            B[i] = 0.0;
            X[i] = 0.0; 
    }
    A[0][0] = 0.0 ;
    A[NX+1][NX+1] = 1.0 ; //ghost points
    //三個點需討論x[1]x[2]...x[NX]
    for(int i = 3 ; i <= NX-1 ; i++){
    	A[i][i] = a2 ;//本點
		A[i][i+1] = a2 ;//東側計算點
		A[i][i-1] = -a4 ;//西側計算點
		A[i][i-2] = a1 ;//西西側計算點 
		B[i] = 0.0 ; 
	}
	A[1][1] = a4 ;//最左側計算點 (7/8)*dx
	A[1][2] = a2 ;//(3/8)*dx
	B[1] = dx*(10.0/8.0)*T_left ;
	A[2][3] = a2;//最左側第二排計算點 //(3/8)*dx
	A[2][2] = a2;//(3/8)*dx 
	A[2][1] = -dx;
	B[2] = dx*(-2.0/8.0)*T_left ;
	A[NX][NX] = a2;//最右側計算點 //(-3.0/8.0)*dx
	A[NX][NX-1] = a3;
	A[NX][NX-2] = -a1;
	B[NX] = dx*T_right ;
	cout << "賦值完畢!" << endl ; 
}

void Jacobiup(vector<vector<double> >& A, vector<double>& B, vector<double>& X, int NX) {
	for(int k = 1; k <= NX ; k++) {
        X_old[k] = X[k] ;
    }
	for(int k = 1; k <= NX ; k++) {
        double sum = 0;
        for(int p = 1; p <= NX; p++) {
            if(p != k) {
                sum += A[k][p] * X[p];
            }
        }
        double X_new = (B[k] - sum) / A[k][k];
        X[k] = X_old[k] + 0.1 * (X_new - X_old[k]);
    }
    
    //計算當前步最大誤差
    maxerror2 = 0;
    for(int k = 1; k <= NX ; k++) {
        double error1 = fabs(X[k] - X_old[k]);
        if(maxerror2 < error1) {
            maxerror2 = error1 ;
        }
    }
}





//係數矩陣賦值（考慮Neumann邊界條件）
void Gamma0(vector<vector<double> >& a, vector<double>& b , vector<double>& x) {
	// 初始化所有元素為0
    for(int i = 1; i <= n-NX ; i++) {
        for(int j = 1; j <= n-NX ; j++) {
            a[i][j] = 0.0 ;
        }  
        b[i] = 0.0 ;
        x[i] = 0.0 ;
    }
	a[0][0] = 0.0 ;
	a[n-NX+1][n-NX+1] = 1.0 ;
	//記得不要每一次迭代都把 解賦值為0 
    ///////////////////////////////////////////////////////////////  
    //第一角點 (1,1) - 左下角
    int p1 = (1-1)*NX + 1;
    a[p1][p1] = a4+a4 ;     //本點
    a[p1][p1+1] = a2;      //東側計算點
    a[p1][p1+NX] = a2;     //北側計算點
    b[p1] = dx*(10.0/8.0)*(T_bottom+T_left);
    
    //第二角點 (1,2)
    int p2 = (1-1)*NX + 2;
    a[p2][p2] = a2+a4;     //本點
    a[p2][p2+1] = a2;      //東側計算點
    a[p2][p2-1] = -dx;     //西側計算點
    a[p2][p2+NX] = a2;     //北側計算點
    b[p2] = dx*(-2.0/8.0)*T_left + dx*(10.0/8.0)*T_bottom;
    
    //第三角點 (1,NX) - 右下角
    int p3 = (1-1)*NX + NX;
    a[p3][p3] = -a2+a4;    //本點
    a[p3][p3-1] = -a3;     //西側計算點
    a[p3][p3-2] = a1;      //西西側計算點
    a[p3][p3+NX] = a2;     //北側計算點
    b[p3] = -dx*T_right + dx*(10.0/8.0)*T_bottom;
    ///////////////////////////////////////////////////////////////
    //第二排角點處理
    ///////////////////////////////////////////////////////////////
    
    //第四角點 (2,1)
    int p4 = (2-1)*NX + 1;
    a[p4][p4] = a4+a2;     //本點
    a[p4][p4+1] = a2;      //東側計算點
    a[p4][p4+NX] = a2;     //北側計算點
    a[p4][p4-NX] = -dx;    //南側計算點
    b[p4] = dx*(10.0/8.0)*T_left + dx*(-2.0/8.0)*T_bottom;
    
    //第五角點 (2,2)
    int p5 = (2-1)*NX + 2;
    a[p5][p5] = a2+a2;     //本點
    a[p5][p5+1] = a2;      //東側計算點
    a[p5][p5-1] = -dx;     //西側計算點
    a[p5][p5+NX] = a2;     //北側計算點
    a[p5][p5-NX] = -dx;    //南側計算點
    b[p5] = dx*(-2.0/8.0)*T_left + dx*(-2.0/8.0)*T_bottom;
    
    //第六角點 (2,NX)
    int p6 = (2-1)*NX + NX;
    a[p6][p6] = -a2+a2 + dx ; //本點 (對角係數 = 0)
    a[p6][p6-1] = -a3;     //西側計算點
    a[p6][p6-2] = a1;      //西西側計算點
    a[p6][p6+NX] = a2;     //北側計算點
    a[p6][p6-NX] = -dx;    //南側計算點
    b[p6] = -dx*T_right + dx*(-2.0/8.0)*T_bottom;
    
    ///////////////////////////////////////////////////////////////
    //上邊界角點處理（最後一排）- Neumann條件
    ///////////////////////////////////////////////////////////////
    
    //第七個角點 (1,NX-1)
    int p7 = (NX-2)*NX + 1;
    a[p7][p7] = a4 - a2 ;        //本點（原始係數組合）
    a[p7][p7+1] = a2 ;            //東側計算點
    a[p7][p7+NX] = 0.0;
    a[p7][p7-NX] = -a3;          //南側計算點
    a[p7][p7-2*NX] = a1;         //南南側計算點
    b[p7] = dx*(10.0/8.0)*T_left - dx*T[1-1][NX-1] ;
    
    ///第八個角點 (2,NX-1)
    int p8 = (NX-2)*NX + 2;
    // 原本北側的a2項，因為T_N = T_{N-1}，所以移到對角線
    a[p8][p8] = a2-a2+dx ;         //本點（a2原有 + a2從北側移來）
    a[p8][p8+1] = a2;            //東側計算點
    a[p8][p8-1] = -dx;           //西側計算點
    a[p8][p8+NX] = 0.0;
    a[p8][p8-NX] = -a3;          //南側計算點
    a[p8][p8-2*NX] = a1;         //南南側計算點
    b[p8] = dx*(-2.0/8.0)*T_left- dx*T[2-1][NX-1] ;
    
    //第九個角點 (NX, NX-1)
    int p9 = (NX-2)*NX + NX;
    a[p9][p9] = -a2 -a2 ;         //本點（a2原有 + a2從北側移來）
    a[p9][p9-1] = -a3 ;           //西側計算點
    a[p9][p9-2] = a1 ; 
    a[p9][p9+NX] = 0.0;
    a[p9][p9-NX] = -a3;          //南側計算點
    a[p9][p9-2*NX] = a1;         //南南側計算點
    b[p9] = -dx*T_right -dx*T[NX-1][NX-1] ;
    
    /////////////////////////////////////////////////////////////////////
    // 邊界處理（不含角點）
    /////////////////////////////////////////////////////////////////////
    
    // 下邊界 (除角點外)
    for(int i = 3; i <= NX-1; i++) {
        int index = (1-1)*NX + i;
        a[index][index] = a2+a4;   //本點
        a[index][index+1] = a2;     //東側計算點
        a[index][index-1] = -a4;    //西側計算點
        a[index][index-2] = a1;     //西西側計算點
        a[index][index+NX] = a2;    //北側計算點
        b[index] = dx*(10.0/8.0)*T_bottom;
    }
    
    //上第二排邊界(i = 3~NX-1 ; j = NX-1) 
    for(int i = 3; i <= NX-1; i++) {
        int index = (NX-2)*NX + i;
        a[index][index] = a2 -a2 + dx ;       //本點（注意：無上側項）
        a[index][index+1] = a2;     //東側計算點
        a[index][index-1] = -a4;    //西側計算點
        a[index][index-2] = a1;     //西西側計算點
        a[index][index+NX] = 0.0 ;
        a[index][index-NX] = -a3 ;
        a[index][index-2*NX] = a1 ;
        b[index] = -dx*T[i-1][NX-1];          //Neumann邊界無額外源項
    }
    // 左邊界 (除角點外)
    for(int j = 3; j <= NX-2; j++) {
        int index = (j-1)*NX + 1;
        a[index][index] = a4+a2;    //本點
        a[index][index+1] = a2;      //東側計算點
        a[index][index-NX] = -a4;    //南側計算點
        a[index][index+NX] = a2;     //北側計算點
        a[index][index-2*NX] = a1;   //南南側計算點
        b[index] = dx*(10.0/8.0)*T_left;
    }
    
    // 右邊界 (除角點外)
    for(int j = 3; j <= NX-2; j++) {
        int index = (j-1)*NX + NX;
        a[index][index] = -a2+a2 + dx ;   //本點
        a[index][index-1] = -a3;     //西側計算點
        a[index][index-2] = a1;      //西西側計算點
        a[index][index+NX] = a2;     //北側計算點
        a[index][index-NX] = -a4;    //南側計算點
        a[index][index-2*NX] = a1;   //南南側計算點
        b[index] = -dx*T_right;
    }
    
    // 左邊第二層邊界
    for(int j = 3; j <= NX-2; j++) {
        int index = (j-1)*NX + 2;
        a[index][index] = a2+a2;     //本點
        a[index][index+1] = a2;       //東側計算點
        a[index][index-1] = -dx;      //西側計算點
        a[index][index+NX] = a2;      //北側計算點
        a[index][index-NX] = -a4;     //南側計算點
        a[index][index-2*NX] = a1;    //南南側計算點
        b[index] = dx*(-2.0/8.0)*T_left;
    }
    
    // 下面第二層邊界
    for(int i = 3; i <= NX-1; i++){
        int index = (2-1)*NX + i;
        a[index][index] = a2+a2;     //本點
        a[index][index+1] = a2;       //東側計算點
        a[index][index-1] = -a4 ;      //西側計算點
        a[index][index-2] = a1;       //西西側計算點
        a[index][index+NX] = a2;      //北側計算點
        a[index][index-NX] = -dx ;     //南側計算點
        b[index] = dx*(-2.0/8.0)*T_bottom;
    }
    
    // 內點
    for(int i = 3; i <= NX-1; i++) {
        for(int j = 3; j <= NX-2; j++) {
            int index = NX*(j-1)+i;
            a[index][index] = a2+a2;     //本點
            a[index][index-1] = -a4;      //西側計算點
            a[index][index+1] = a2;       //東側計算點
            a[index][index-2] = a1;       //西西側計算點
            a[index][index-NX] = -a4;     //南側計算點
            a[index][index+NX] = a2;      //北側計算點
            a[index][index-2*NX] = a1;    //南南側計算點
            b[index] = 0.0;
        }
    }
}
//Jacobi迭代求解
void Jacobi(vector<vector<double> >& a, vector<double>& b, vector<double>& x, int n) {
    for(int i = 1 ; i <= n-NX ; i++){
    	x_old[i] = x[i] ;
	}
    for(int k = 1; k <= n-NX ; k++) {
        if(fabs(a[k][k]) < 1e-15) continue; // 跳過奇異矩陣
        
        double sum = 0;
        for(int p = 1; p <= n-NX ; p++) {
            if(p != k) {
                sum += a[k][p] * x[p];
            }
        }
        double x_new = (b[k] - sum) / a[k][k];
        x[k] = x_old[k] + lamda * (x_new - x_old[k]);
    }
    
    //計算當前步最大誤差
    maxerror = 0;
    for(int k = 1; k <= n-NX; k++) {
        double error = fabs(x[k] - x_old[k]);
        if(maxerror < error) {
            maxerror = error;
        }
    }//計算i = 1 : NX ; j = 1 : NX-1 的最大誤差 
}

//輸出VTK檔案
void output(int m) {
    // 將一維解轉換為二維溫度場
    for(int j = 1; j <= NX-1; j++) {
        for(int i = 1; i <= NX; i++) {
            T[i-1][j-1] = x[(j-1)*NX + i];
        }
    }
    
    ostringstream name;
    name << "QUICK上邊界計算點先算"<< NX+1 << "x" << NX+1 << "步數 = "<< setfill('0') << setw(6)  <<  m << ".vtk";
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
        
        // QUICK格式係數定義
        a1 = dx*(1.0/8.0);  //係數1/8
        a2 = dx*(3.0/8.0);  //係數3/8
        a3 = dx*(6.0/8.0);  //係數6/8
        a4 = dx*(7.0/8.0);  //係數7/8
        
        n = (NX)*(NX);
        
        // 重新調整向量大小
        a.assign(n+2, vector<double>(n+2, 0.0));
        b.assign(n+2, 0.0);
        x.assign(n+2, 0.0);
        x_old.assign(n+2, 0.0);
        T.assign(NX, vector<double>(NX, 0.0));
        A.assign(NX+2, vector<double>(NX+2, 0.0));
        B.assign(NX+2, 0.0);
        X.assign(NX+2, 0.0);
        X_old.assign(NX+2, 0.0);
        cout << "========================================" << endl;
        cout << "程式碼開始執行...." << endl;
        cout << "網格大小: " << NX+1 << " x " << NX+1 << endl;
        cout << "網格間距: dx = dy = " << dx << endl;
        cout << "邊界條件設定：" << endl;
        cout << "左邊界(Dirichlet): T = " << T_left << endl;
        cout << "右邊界(Dirichlet): T = " << T_right << endl;
        cout << "下邊界(Dirichlet): T = " << T_bottom << endl;
        cout << "上邊界(Neumann): \partial T/\partial n = 0" << endl;
        cout << "鬆弛因子: λ = " << lamda << endl;
        cout << "收斂準則: " << tolerance << endl;
        cout << "========================================" << endl;
        steadystate = false;
        cofficient(A,B,X) ; //一維場點係數賦值 
	    int K_max = 1000000 ;
	    Temp = 0 ;
	    for(int K= 0 ; K<=K_max ; K++){
		    Jacobiup(A,B,X,NX) ;
		    if(K>0 && maxerror2 < 0.1){
			    cout << "得到上層溫度場!"<< endl ;
			    for(int i = 1 ; i <= NX ; i++){
				    T[i-1][NX-1] = X[i] ; //賦與最上層溫度值 
			    }
			Temp = 1 ;
			break ; 
	    	}
     	}
	    if(!Temp){
		    cout << "最上層溫度場不收斂!"<< endl ;
	    }
        Gamma0(a, b ,x);         //設定係數矩陣（含邊界條件）
        
        for(G = 0; G < max_G; G++) {
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
        cout << "網格 " << NX << "x" << NX << " 計算完成\n" << endl;
    }
    
    cout << "所有計算完成！" << endl;
    return 0;
}
