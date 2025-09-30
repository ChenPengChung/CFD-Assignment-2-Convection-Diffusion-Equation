//�Q��QUICK�榡�D�ѤG��í�A�¹�y��{��
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>
#include <iomanip>
using namespace std;

//�|�Ψ쪺�Ѽ�
int Nx[] = {40,80};
//������n�k�ĥ�Wet node�B�z�A�G�p���I�ƥ� = ����ƥ�, �`�I�ƥ� = ����ƥ�+1
int NX; //NX = NY => dx = dy
int n;
vector<vector<double> > a; //a�G���x�}�����������Y�Ưx�}
vector<double> b; //b�@���x�}�������k�����~�����A��{�դ��D���ʶ�
vector<double> x;
vector<double> x_old;
vector<vector<double> > T;
vector<vector<double> > A; //a�G���x�}�����������Y�Ưx�}
vector<double> B; //b�@���x�}�������k�����~�����A��{�դ��D���ʶ�
vector<double> X;
vector<double> X_old;
double dx; //�@�Ӻ��涡�Z
//��ɱ���w�q
const double T_left = 1.0;    //�����Dirichlet����
const double T_right = 0.0;   //�k���Dirichlet����  
const double T_bottom = 0.0;  //�U���Dirichlet����
//�W��ɬ�Neumann����: ?T/?n = 0

//QUICK�榡�Y��
double a1,a2,a3,a4;

//���N�Ѽ�
const double lamda = 0.3; //�v���W�P�����N�k
int G, max_G = 10000000;
double maxerror , maxerror2 ; 
const float tolerance = 1e-10;
bool Temp , steadystate;


//���F���� Neumann����A���D�ѳ̤W�h��ɧ@���@���ū׳����D��
//�D�ѧ�����A�W�h��ɪ��ū׳��Y�@���w�����ūצA�ӭ��N�D�ѤU�����I 
//�W�h����I���p�� ����@���N�n 
//�Q�γ̤W�h��ɭp���I�@���@�����I�D�ѷū׳�
void cofficient(vector<vector<double> >& A, vector<double>& B , vector<double>& X){
	// ��l�ƩҦ�������0
    for(int i = 1; i <= NX ; i++) {
            A[i][i] = 0.0 ;
            B[i] = 0.0;
            X[i] = 0.0; 
    }
    A[0][0] = 0.0 ;
    A[NX+1][NX+1] = 1.0 ; //ghost points
    //�T���I�ݰQ��x[1]x[2]...x[NX]
    for(int i = 3 ; i <= NX-1 ; i++){
    	A[i][i] = a2 ;//���I
		A[i][i+1] = a2 ;//�F���p���I
		A[i][i-1] = -a4 ;//�谼�p���I
		A[i][i-2] = a1 ;//��谼�p���I 
		B[i] = 0.0 ; 
	}
	A[1][1] = a4 ;//�̥����p���I (7/8)*dx
	A[1][2] = a2 ;//(3/8)*dx
	B[1] = dx*(10.0/8.0)*T_left ;
	A[2][3] = a2;//�̥����ĤG�ƭp���I //(3/8)*dx
	A[2][2] = a2;//(3/8)*dx 
	A[2][1] = -dx;
	B[2] = dx*(-2.0/8.0)*T_left ;
	A[NX][NX] = a2;//�̥k���p���I //(-3.0/8.0)*dx
	A[NX][NX-1] = a3;
	A[NX][NX-2] = -a1;
	B[NX] = dx*T_right ;
	cout << "��ȧ���!" << endl ; 
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
    
    //�p���e�B�̤j�~�t
    maxerror2 = 0;
    for(int k = 1; k <= NX ; k++) {
        double error1 = fabs(X[k] - X_old[k]);
        if(maxerror2 < error1) {
            maxerror2 = error1 ;
        }
    }
}





//�Y�Ưx�}��ȡ]�Ҽ{Neumann��ɱ���^
void Gamma0(vector<vector<double> >& a, vector<double>& b , vector<double>& x) {
	// ��l�ƩҦ�������0
    for(int i = 1; i <= n-NX ; i++) {
        for(int j = 1; j <= n-NX ; j++) {
            a[i][j] = 0.0 ;
        }  
        b[i] = 0.0 ;
        x[i] = 0.0 ;
    }
	a[0][0] = 0.0 ;
	a[n-NX+1][n-NX+1] = 1.0 ;
	//�O�o���n�C�@�����N���� �ѽ�Ȭ�0 
    ///////////////////////////////////////////////////////////////  
    //�Ĥ@���I (1,1) - ���U��
    int p1 = (1-1)*NX + 1;
    a[p1][p1] = a4+a4 ;     //���I
    a[p1][p1+1] = a2;      //�F���p���I
    a[p1][p1+NX] = a2;     //�_���p���I
    b[p1] = dx*(10.0/8.0)*(T_bottom+T_left);
    
    //�ĤG���I (1,2)
    int p2 = (1-1)*NX + 2;
    a[p2][p2] = a2+a4;     //���I
    a[p2][p2+1] = a2;      //�F���p���I
    a[p2][p2-1] = -dx;     //�谼�p���I
    a[p2][p2+NX] = a2;     //�_���p���I
    b[p2] = dx*(-2.0/8.0)*T_left + dx*(10.0/8.0)*T_bottom;
    
    //�ĤT���I (1,NX) - �k�U��
    int p3 = (1-1)*NX + NX;
    a[p3][p3] = -a2+a4;    //���I
    a[p3][p3-1] = -a3;     //�谼�p���I
    a[p3][p3-2] = a1;      //��谼�p���I
    a[p3][p3+NX] = a2;     //�_���p���I
    b[p3] = -dx*T_right + dx*(10.0/8.0)*T_bottom;
    ///////////////////////////////////////////////////////////////
    //�ĤG�ƨ��I�B�z
    ///////////////////////////////////////////////////////////////
    
    //�ĥ|���I (2,1)
    int p4 = (2-1)*NX + 1;
    a[p4][p4] = a4+a2;     //���I
    a[p4][p4+1] = a2;      //�F���p���I
    a[p4][p4+NX] = a2;     //�_���p���I
    a[p4][p4-NX] = -dx;    //�n���p���I
    b[p4] = dx*(10.0/8.0)*T_left + dx*(-2.0/8.0)*T_bottom;
    
    //�Ĥ����I (2,2)
    int p5 = (2-1)*NX + 2;
    a[p5][p5] = a2+a2;     //���I
    a[p5][p5+1] = a2;      //�F���p���I
    a[p5][p5-1] = -dx;     //�谼�p���I
    a[p5][p5+NX] = a2;     //�_���p���I
    a[p5][p5-NX] = -dx;    //�n���p���I
    b[p5] = dx*(-2.0/8.0)*T_left + dx*(-2.0/8.0)*T_bottom;
    
    //�Ĥ����I (2,NX)
    int p6 = (2-1)*NX + NX;
    a[p6][p6] = -a2+a2 + dx ; //���I (�﨤�Y�� = 0)
    a[p6][p6-1] = -a3;     //�谼�p���I
    a[p6][p6-2] = a1;      //��谼�p���I
    a[p6][p6+NX] = a2;     //�_���p���I
    a[p6][p6-NX] = -dx;    //�n���p���I
    b[p6] = -dx*T_right + dx*(-2.0/8.0)*T_bottom;
    
    ///////////////////////////////////////////////////////////////
    //�W��ɨ��I�B�z�]�̫�@�ơ^- Neumann����
    ///////////////////////////////////////////////////////////////
    
    //�ĤC�Ө��I (1,NX-1)
    int p7 = (NX-2)*NX + 1;
    a[p7][p7] = a4 - a2 ;        //���I�]��l�Y�ƲզX�^
    a[p7][p7+1] = a2 ;            //�F���p���I
    a[p7][p7+NX] = 0.0;
    a[p7][p7-NX] = -a3;          //�n���p���I
    a[p7][p7-2*NX] = a1;         //�n�n���p���I
    b[p7] = dx*(10.0/8.0)*T_left - dx*T[1-1][NX-1] ;
    
    ///�ĤK�Ө��I (2,NX-1)
    int p8 = (NX-2)*NX + 2;
    // �쥻�_����a2���A�]��T_N = T_{N-1}�A�ҥH����﨤�u
    a[p8][p8] = a2-a2+dx ;         //���I�]a2�즳 + a2�q�_�����ӡ^
    a[p8][p8+1] = a2;            //�F���p���I
    a[p8][p8-1] = -dx;           //�谼�p���I
    a[p8][p8+NX] = 0.0;
    a[p8][p8-NX] = -a3;          //�n���p���I
    a[p8][p8-2*NX] = a1;         //�n�n���p���I
    b[p8] = dx*(-2.0/8.0)*T_left- dx*T[2-1][NX-1] ;
    
    //�ĤE�Ө��I (NX, NX-1)
    int p9 = (NX-2)*NX + NX;
    a[p9][p9] = -a2 -a2 ;         //���I�]a2�즳 + a2�q�_�����ӡ^
    a[p9][p9-1] = -a3 ;           //�谼�p���I
    a[p9][p9-2] = a1 ; 
    a[p9][p9+NX] = 0.0;
    a[p9][p9-NX] = -a3;          //�n���p���I
    a[p9][p9-2*NX] = a1;         //�n�n���p���I
    b[p9] = -dx*T_right -dx*T[NX-1][NX-1] ;
    
    /////////////////////////////////////////////////////////////////////
    // ��ɳB�z�]���t���I�^
    /////////////////////////////////////////////////////////////////////
    
    // �U��� (�����I�~)
    for(int i = 3; i <= NX-1; i++) {
        int index = (1-1)*NX + i;
        a[index][index] = a2+a4;   //���I
        a[index][index+1] = a2;     //�F���p���I
        a[index][index-1] = -a4;    //�谼�p���I
        a[index][index-2] = a1;     //��谼�p���I
        a[index][index+NX] = a2;    //�_���p���I
        b[index] = dx*(10.0/8.0)*T_bottom;
    }
    
    //�W�ĤG�����(i = 3~NX-1 ; j = NX-1) 
    for(int i = 3; i <= NX-1; i++) {
        int index = (NX-2)*NX + i;
        a[index][index] = a2 -a2 + dx ;       //���I�]�`�N�G�L�W�����^
        a[index][index+1] = a2;     //�F���p���I
        a[index][index-1] = -a4;    //�谼�p���I
        a[index][index-2] = a1;     //��谼�p���I
        a[index][index+NX] = 0.0 ;
        a[index][index-NX] = -a3 ;
        a[index][index-2*NX] = a1 ;
        b[index] = -dx*T[i-1][NX-1];          //Neumann��ɵL�B�~����
    }
    // ����� (�����I�~)
    for(int j = 3; j <= NX-2; j++) {
        int index = (j-1)*NX + 1;
        a[index][index] = a4+a2;    //���I
        a[index][index+1] = a2;      //�F���p���I
        a[index][index-NX] = -a4;    //�n���p���I
        a[index][index+NX] = a2;     //�_���p���I
        a[index][index-2*NX] = a1;   //�n�n���p���I
        b[index] = dx*(10.0/8.0)*T_left;
    }
    
    // �k��� (�����I�~)
    for(int j = 3; j <= NX-2; j++) {
        int index = (j-1)*NX + NX;
        a[index][index] = -a2+a2 + dx ;   //���I
        a[index][index-1] = -a3;     //�谼�p���I
        a[index][index-2] = a1;      //��谼�p���I
        a[index][index+NX] = a2;     //�_���p���I
        a[index][index-NX] = -a4;    //�n���p���I
        a[index][index-2*NX] = a1;   //�n�n���p���I
        b[index] = -dx*T_right;
    }
    
    // ����ĤG�h���
    for(int j = 3; j <= NX-2; j++) {
        int index = (j-1)*NX + 2;
        a[index][index] = a2+a2;     //���I
        a[index][index+1] = a2;       //�F���p���I
        a[index][index-1] = -dx;      //�谼�p���I
        a[index][index+NX] = a2;      //�_���p���I
        a[index][index-NX] = -a4;     //�n���p���I
        a[index][index-2*NX] = a1;    //�n�n���p���I
        b[index] = dx*(-2.0/8.0)*T_left;
    }
    
    // �U���ĤG�h���
    for(int i = 3; i <= NX-1; i++){
        int index = (2-1)*NX + i;
        a[index][index] = a2+a2;     //���I
        a[index][index+1] = a2;       //�F���p���I
        a[index][index-1] = -a4 ;      //�谼�p���I
        a[index][index-2] = a1;       //��谼�p���I
        a[index][index+NX] = a2;      //�_���p���I
        a[index][index-NX] = -dx ;     //�n���p���I
        b[index] = dx*(-2.0/8.0)*T_bottom;
    }
    
    // ���I
    for(int i = 3; i <= NX-1; i++) {
        for(int j = 3; j <= NX-2; j++) {
            int index = NX*(j-1)+i;
            a[index][index] = a2+a2;     //���I
            a[index][index-1] = -a4;      //�谼�p���I
            a[index][index+1] = a2;       //�F���p���I
            a[index][index-2] = a1;       //��谼�p���I
            a[index][index-NX] = -a4;     //�n���p���I
            a[index][index+NX] = a2;      //�_���p���I
            a[index][index-2*NX] = a1;    //�n�n���p���I
            b[index] = 0.0;
        }
    }
}
//Jacobi���N�D��
void Jacobi(vector<vector<double> >& a, vector<double>& b, vector<double>& x, int n) {
    for(int i = 1 ; i <= n-NX ; i++){
    	x_old[i] = x[i] ;
	}
    for(int k = 1; k <= n-NX ; k++) {
        if(fabs(a[k][k]) < 1e-15) continue; // ���L�_���x�}
        
        double sum = 0;
        for(int p = 1; p <= n-NX ; p++) {
            if(p != k) {
                sum += a[k][p] * x[p];
            }
        }
        double x_new = (b[k] - sum) / a[k][k];
        x[k] = x_old[k] + lamda * (x_new - x_old[k]);
    }
    
    //�p���e�B�̤j�~�t
    maxerror = 0;
    for(int k = 1; k <= n-NX; k++) {
        double error = fabs(x[k] - x_old[k]);
        if(maxerror < error) {
            maxerror = error;
        }
    }//�p��i = 1 : NX ; j = 1 : NX-1 ���̤j�~�t 
}

//��XVTK�ɮ�
void output(int m) {
    // �N�@�����ഫ���G���ū׳�
    for(int j = 1; j <= NX-1; j++) {
        for(int i = 1; i <= NX; i++) {
            T[i-1][j-1] = x[(j-1)*NX + i];
        }
    }
    
    ostringstream name;
    name << "QUICK�W��ɭp���I����"<< NX+1 << "x" << NX+1 << "�B�� = "<< setfill('0') << setw(6)  <<  m << ".vtk";
    ofstream out(name.str().c_str());
    
    // VTK ����Y
    out << "# vtk DataFile Version 3.0\n";
    out << "QUICK_Gamma=0.0\n";
    out << "ASCII\n";
    out << "DATASET STRUCTURED_POINTS\n";
    out << "DIMENSIONS " << NX << " " << NX << " 1\n";
    out << "ORIGIN 0 0 0\n";
    out << "SPACING " << dx << " " << dx << " 1\n";
    out << "POINT_DATA " << NX*NX << "\n";
    
    // ��X�ū׳�
    out << "SCALARS Temperature double 1\n";
    out << "LOOKUP_TABLE default\n";
    for(int j = 0; j < NX; j++) {
        for(int i = 0; i < NX; i++) {
            out << scientific << setprecision(6) << T[i][j] << "\n";
        }
    }
    out.close();
    cout << "VTK ���w��X: " << name.str() << endl;
}

int main() {
    for(int grid_idx = 0; grid_idx < 2; grid_idx++){
        NX = Nx[grid_idx];
        dx = 1.0/double(NX);
        
        // QUICK�榡�Y�Ʃw�q
        a1 = dx*(1.0/8.0);  //�Y��1/8
        a2 = dx*(3.0/8.0);  //�Y��3/8
        a3 = dx*(6.0/8.0);  //�Y��6/8
        a4 = dx*(7.0/8.0);  //�Y��7/8
        
        n = (NX)*(NX);
        
        // ���s�վ�V�q�j�p
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
        cout << "�{���X�}�l����...." << endl;
        cout << "����j�p: " << NX+1 << " x " << NX+1 << endl;
        cout << "���涡�Z: dx = dy = " << dx << endl;
        cout << "��ɱ���]�w�G" << endl;
        cout << "�����(Dirichlet): T = " << T_left << endl;
        cout << "�k���(Dirichlet): T = " << T_right << endl;
        cout << "�U���(Dirichlet): T = " << T_bottom << endl;
        cout << "�W���(Neumann): \partial T/\partial n = 0" << endl;
        cout << "�P���]�l: �f = " << lamda << endl;
        cout << "���ķǫh: " << tolerance << endl;
        cout << "========================================" << endl;
        steadystate = false;
        cofficient(A,B,X) ; //�@�����I�Y�ƽ�� 
	    int K_max = 1000000 ;
	    Temp = 0 ;
	    for(int K= 0 ; K<=K_max ; K++){
		    Jacobiup(A,B,X,NX) ;
		    if(K>0 && maxerror2 < 0.1){
			    cout << "�o��W�h�ū׳�!"<< endl ;
			    for(int i = 1 ; i <= NX ; i++){
				    T[i-1][NX-1] = X[i] ; //��P�̤W�h�ū׭� 
			    }
			Temp = 1 ;
			break ; 
	    	}
     	}
	    if(!Temp){
		    cout << "�̤W�h�ū׳�������!"<< endl ;
	    }
        Gamma0(a, b ,x);         //�]�w�Y�Ưx�}�]�t��ɱ���^
        
        for(G = 0; G < max_G; G++) {
            Jacobi(a, b, x, n);  //����@��Jacobi���N
            if(G % 100 == 0) {
                cout << "���N���� = " << G;
                if(G > 0) {
                    cout << ", �̤j�ܤƶq = " << scientific << maxerror;
                }
                cout << endl;
                output(G);
            }
            
            if(G > 100 && maxerror < tolerance) {
                steadystate = true;
                cout << "\n�w�F���ı���I" << endl;
                cout << "�̲׭��N����: " << G << endl;
                cout << "�̲׻~�t: " << scientific << maxerror << endl;
                break;
            }
        }
        
        if(!steadystate) {
            cout << "\nĵ�i�G�F��̤j���N���ơA�����F��í�A�I" << endl;
        }
        
        output(G);  //��X�̲׵��G
        cout << "���� " << NX << "x" << NX << " �p�⧹��\n" << endl;
    }
    
    cout << "�Ҧ��p�⧹���I" << endl;
    return 0;
}
