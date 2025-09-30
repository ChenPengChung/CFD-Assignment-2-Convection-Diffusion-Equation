//�Q�ΪŶ��G�����MINMOD�榡�D�ѤG��í�A�X����y��{
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
vector<vector<double> > r_w;
vector<vector<double> > r_e;
vector<vector<double> > r_s;
vector<vector<double> > r_n;
double dx; //�@�Ӻ��涡�Z
//��ɱ���w�q
const double T_left = 1.0;    //�����Dirichlet����
const double T_right = 0.0;   //�k���Dirichlet����  
const double T_bottom = 0.0;  //�U���Dirichlet����
const double Pelect = 1000.0 ;
const double Gamma = sqrt(2.0)/Pelect ;//���X���Y��
//���N�Ѽ�
const double lamda = 0.003; //�v���W�P�����N�k
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
    // ��l�ƩҦ�������0
    for(int i = 1; i <= NX*NX ; i++) {
        for(int j = 1; j <= NX*NX; j++) {
            a[i][j] = 0.0 ;
        }  
        b[i] = 0.0 ;
        x[i] = 0.0 ;
    }//�O�o���n�C�@�����N���� �ѽ�Ȭ�0 
    a[0][0] = 0.0 ;
    a[NX*NX+1][NX*NX+1] = 1.0 ;
    //��l�Ʒū׳�
    //����ɭp���I��������
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

//�Y�Ưx�}��ȡ]�Ҽ{Neumann��ɱ���^
void Gamma0(vector<vector<double> >& a, vector<double>& b , vector<double>& x) { //�C�@�����N�n�����Ʊ��N�O�Na,b�x�}���s���
    for(int i = 1; i <= NX*NX ; i++) {
        for(int j = 1; j <= NX*NX; j++) {
            a[i][j] = 0.0 ;
        }  
        b[i] = 0.0 ;
    }//�O�o���n�C�@�����N���� �ѽ�Ȭ�0 
    a[0][0] = 0.0 ;
    a[NX*NX+1][NX*NX+1] = 1.0 ;
    ///////////////////////////////////////////////////////////////  
    //�Ĥ@���I (1,1) - ���U��
    int p1 = (1-1)*NX + 1;
    a[p1][p1] =(2*Gamma+dx*(U(r_n[1-1][1-1])))+(2*Gamma + dx*U(r_e[1-1][1-1]));//���I
    //a[p1][p1-1] =- ;//�谼�p���I
    a[p1][p1+1] =-(Gamma - dx*D(r_e[1-1][1-1])) ;//�F���p���I
    //a[p1][p1-NX] =-;//�n���p���I
    a[p1][p1+NX] =-(Gamma-dx*D(r_n[1-1][1-1]));//�_���p���I
    b[p1] = dx*T_bottom + dx*T_left + (Gamma - dx*UU(r_e[1-1][1-1])) * T[0][1] + (Gamma-dx*UU(r_n[1-1][1-1])) * T[1][0];
    //�ĤG���I (NX,1) - �k�U��
    int p2 = (1-1)*NX + NX;
    a[p2][p2] =(2*Gamma+dx*(U(r_n[NX-1][1-1])))+(2*Gamma+dx*(-D(r_w[NX-1][1-1]))) ;//���I
    a[p2][p2-1] =-(Gamma-dx*(-U(r_w[NX-1][1-1]))) ;//�谼�p���I
    //a[p2][p2+1] =- ;//�F���p���I
    a[p2][p2-2] =-(-dx*(-UU(r_w[NX-1][1-1]))) ;//��谼�p���I
    //a[p2][p2-NX] =-;//�n���p���I
    a[p2][p2+NX] =-(Gamma-dx*D(r_n[NX-1][1-1]));//�_���p���I
    b[p2] = dx*T_bottom-dx*T_right + (Gamma) * T[NX+1][1] + (Gamma-dx*UU(r_n[NX-1][1-1])) * T[NX][0];
    //�ĤT�Ө��I (1,NX) - ���W��
    int p3 = (NX-1)*NX + 1;
    a[p3][p3] =(1*Gamma+dx*(U(r_n[1-1][NX-1])-D(r_s[1-1][NX-1]))) + (2*Gamma + dx*U(r_e[1-1][NX-1]));//���I
    //a[p3][p3-1] =- ;//�谼�p���I
    a[p3][p3+1] =-(Gamma - dx*D(r_e[1-1][NX-1])) ;//�F���p���I
    a[p3][p3-NX] =-(Gamma-dx*(UU(r_n[1-1][NX-1])-U(r_s[1-1][NX-1])));//�n���p���I
    //a[p3][p3+NX] =-;//�_���p���I
    a[p3][p3-2*NX] =-(-dx*(-UU(r_s[1-1][NX-1])));//�n�n���p���I
    b[p3] = 0.0 + dx*T_left + (Gamma - dx*UU(r_e[1-1][NX-1])) * T[0][NX] + (-dx*D(r_n[1-1][NX-1]))*T[1][NX+1];
    //�ĥ|�Ө��I (NX, NX) - �k�W��
    int p4 = (NX-1)*NX + NX;
    a[p4][p4] =(1*Gamma+dx*(U(r_n[NX-1][NX-1])-D(r_s[NX-1][NX-1])))+ (2*Gamma+dx*(-D(r_w[NX-1][NX-1]))); //���I
    a[p4][p4-1] =-(Gamma-dx*(-U(r_w[NX-1][NX-1]))) ;//�谼�p���I
    //a[p4][p4+1] =- ;//�F���p���I
    a[p4][p4-2] =-(-dx*(-UU(r_w[NX-1][NX-1]))) ;//��谼�p���I
    a[p4][p4-NX] =-(Gamma-dx*(UU(r_n[NX-1][NX-1])-U(r_s[NX-1][NX-1])));//�n���p���I
    //a[p4][p4+NX] =-;//�_���p���I
    a[p4][p4-2*NX] =-(-dx*(-UU(r_s[NX-1][NX-1])));//�n�n���p���I
    b[p4] = 0.0 -dx*T_right + (Gamma) * T[NX+1][NX] + (-dx*D(r_n[NX-1][NX-1])) * T[NX][NX+1];

    /////////////////////////////////////////////////////////////////////
    // ��ɳB�z�]���t���I�^
    /////////////////////////////////////////////////////////////////////
    
    // �U��� (�����I�~)(i = 3~NX-1 ; j = 1) 
    for(int i = 2; i <= NX-1; i++) {
        int index = (1-1)*NX + i;
        a[index][index] =(4*Gamma+dx*(U(r_e[i-1][1-1])-D(r_w[i-1][1-1]))+dx*(U(r_n[i-1][1-1])));//���I
        a[index][index-1] =-(Gamma-dx*(UU(r_e[i-1][1-1])-U(r_w[i-1][1-1])));//�谼�p���I
        a[index][index+1] =-(Gamma - dx*D(r_e[i-1][1-1]));//�F���p���I
        //a[index][index-NX] =-;//�n���p���I
        a[index][index+NX] =-(Gamma-dx*D(r_n[i-1][1-1]));//�_���p���I
        //a[index][index-2] =-;//��谼�p���I
        b[index] = dx*T_bottom + (Gamma-dx*UU(r_n[i-1][1-1])) * T[i][0] + (-dx*(-UU(r_w[i-1][1-1]))) * T[i-2][1];
    }
    
    //�W���(i = 3~NX-1 ; j = NX) 
    for(int i = 2; i <= NX-1; i++) {
        int index = (NX-1)*NX + i;
        a[index][index] =(3*Gamma+dx*(U(r_e[i-1][NX-1])-D(r_w[i-1][NX-1]))+dx*(U(r_n[i-1][NX-1])-D(r_s[i-1][NX-1]))) ;//���I
        a[index][index-1] =-(Gamma-dx*(UU(r_e[i-1][NX-1])-U(r_w[i-1][NX-1])));//�谼�p���I
        a[index][index+1] =-(Gamma-dx*(D(r_e[i-1][NX-1])));//�F���p���I
        a[index][index-NX] =-(Gamma-dx*(UU(r_n[i-1][NX-1])-U(r_s[i-1][NX-1])));//�n���p���I
        //a[index][index+NX] =-;//�_���p���I
        //a[index][index-2] =-;//��谼�p���I
        a[index][index-2*NX] =-(-dx*(-UU(r_s[i-1][NX-1])));//�n�n���p���I
        b[index] = 0.0 + (-dx*D(r_n[i-1][NX-1])) * T[i][NX+1] + (-dx*(-UU(r_w[i-1][NX-1])))*T[i-2][NX];     
    }
    //Break
    
    // ����� (�����I�~)
    for(int j = 2; j <= NX-1; j++) {
        int index = (j-1)*NX + 1;
        a[index][index] =4*Gamma + dx*U(r_e[1-1][j-1]) + dx*(U(r_n[1-1][j-1])-D(r_s[1-1][j-1])) ;//���I
        //a[index][index-1] =- ;//�谼�p���I
        a[index][index+1] =-(Gamma - dx*D(r_e[1-1][j-1])) ;//�F���p���I
        a[index][index-NX] =-(Gamma - dx*(UU(r_n[1-1][j-1]) - U(r_s[1-1][j-1])) ) ; //�n���p���I
        a[index][index+NX] =-(Gamma - dx*D(r_n[1-1][j-1]) ) ;//�_���p���I
        //a[index][index-2*NX] =- ;//�n�n���p���I
        b[index] = dx*T_left + (Gamma - dx*UU(r_e[1-1][j-1])) * T[0][j] + (-dx*(-UU(r_s[1-1][j-1])) ) * T[1][j-2];
    }
    
    // �k��� (�����I�~)
    for(int j = 2; j <= NX-1; j++) {
        int index = (j-1)*NX + NX;
        a[index][index] =(4*Gamma+dx*(-D(r_w[NX-1][j-1]))+dx*(U(r_n[NX-1][j-1])-D(r_s[NX-1][j-1]))) ; //���I
        a[index][index-1] =-(Gamma-dx*(-U(r_w[NX-1][j-1]))) ;//�谼�p���I
        //a[index][index+1] =- ;//�F���p���I
        a[index][index-NX] =-(Gamma-dx*(UU(r_n[NX-1][j-1])-U(r_s[NX-1][j-1]))) ;//�n���p���I
        a[index][index+NX] =-(Gamma-dx*(D(r_n[NX-1][j-1]))) ;//�_���p���I
        a[index][index-2] =-(-dx*(-UU(r_w[NX-1][j-1]))) ;//��谼�p���I
        //a[index][index-2*NX] =- ;//�n�n���p���I
        b[index] =-dx*T_right + (Gamma)*T[NX+1][j] + (-dx*(-UU(r_s[NX-1][j-1]))) * T[NX][j-2];
    }
    
    
    // ���I
    for(int i = 2; i <= NX-1; i++) {
        for(int j = 2; j <= NX-1; j++) {
            int index = NX*(j-1)+i;
            a[index][index] =4.0*Gamma + dx*(U(r_e[i-1][j-1]) - D(r_w[i-1][j-1])) + dx*(U(r_n[i-1][j-1]) - D(r_s[i-1][j-1]));//���I
            a[index][index+1] = -(Gamma - dx*(D(r_e[i-1][j-1]))) ;//�F���p���I
            a[index][index-1] = -(Gamma - dx*(UU(r_e[i-1][j-1]) - U(r_w[i-1][j-1]))) ;//�谼�p���I
            a[index][index+NX] = -(Gamma - dx*(D(r_n[i-1][j-1]))) ;//�_���p���I
            a[index][index-NX] = -(Gamma - dx*(UU(r_n[i-1][j-1]) - U(r_s[i-1][j-1]))) ;//�n���p���I
            //a[index][index-2] =-;//��谼�p���I
            //a[index][index-2*NX] =-;//�n�n���p���I
            b[index] = 0.0+dx*UU(r_w[i-1][j-1])*T[i-2][j] +dx*UU(r_s[i-1][j-1])*T[i][j-2] ;
        }
    }
}
//Jacobi���N�D��
void Jacobi(vector<vector<double> >& a, vector<double>& b, vector<double>& x, int n) {
    for(int i = 1 ; i <= n ; i++){
    	x_old[i] = x[i] ;
	}
    for(int k = 1; k <= n; k++) {
        if(fabs(a[k][k]) < 1e-15) continue; // ���L�_���x�}
        
        double sum = 0;
        for(int p = 1; p <= n; p++) {
            if(p != k) {
                sum += a[k][p] * x[p];
            }
        }
        double x_new = (b[k] - sum) / a[k][k];
        x[k] = x_old[k] + lamda * (x_new - x_old[k]);
    }
     //��s�ū׳�
    // �N�@�����ഫ���G���ū׳�
    for(int j = 1; j <= NX; j++) {
        for(int i = 1; i <= NX; i++) {
            T[i][j] = x[(j-1)*NX + i];
        }
    }
    //��s�����`�I
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
    //�ū׳��@��s���A�U�x�}�N���W��ۧ�s
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
    //�p���e�B�̤j�~�t
    maxerror = 0;
    for(int k = 1; k <= n; k++) {
        double error = fabs(x[k] - x_old[k]);
        if(maxerror < error) {
            maxerror = error;
        }
    }//�p��i = 1 : NX ; j = 1 : NX-1 ���̤j�~�t 
}

//��XVTK�ɮ�
void output(int m) {
    ostringstream name;
    name << "12.TVD_MINMOD"<< NX+1 << "x" << NX+1 << "�B�� = "<< setfill('0') << setw(6)  <<  m << ".vtk";
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
        n = (NX)*(NX);
        // ���s�վ�V�q�j�p
        a.assign(n+2, vector<double>(n+2, 0.0));
        b.assign(n+2, 0.0);
        x.assign(n+2, 0.0);
        x_old.assign(n+2, 0.0); //�V�q�j�p�����ߤ@ 
        T.assign(NX+2, vector<double>(NX+2, 0.0));
        r_w.assign(NX, vector<double>(NX, 1.0));
        r_e.assign(NX, vector<double>(NX, 1.0));
        r_s.assign(NX, vector<double>(NX, 1.0));
        r_n.assign(NX, vector<double>(NX, 1.0)); //�|�ӯx�}��l�Ȭ�0
        cout << "========================================" << endl;
        cout << "�{���X�}�l����...." << endl;
        cout << "����j�p: " << NX+1 << " x " << NX+1 << endl;
        cout << "���涡�Z: dx = dy = " << dx << endl;
        cout << "��ɱ���]�w�G" << endl;
        cout << "  �����(Dirichlet): T = " << T_left << endl;
        cout << "  �k���(Dirichlet): T = " << T_right << endl;
        cout << "  �U���(Dirichlet): T = " << T_bottom << endl;
        cout << "  �W���(Neumann):  T_[i][NY] = T_[i][NY+1](�����I)" << endl;
        cout << "�P���]�l: �f = " << lamda << endl;
        cout << "���ķǫh: " << tolerance << endl;
        cout << "========================================" << endl;
        steadystate = false;
        initial(a , b , x) ; //��l�ơA�Y�Ưx�}�A�D���ʶ��A�����`�I
        for(G = 0; G < max_G; G++) {
            Gamma0(a, b ,x);         //�]�w�Y�Ưx�}�]�t��ɱ���^
            Jacobi(a, b, x, n);  //����@��Jacobi���N
            if(G % 1000000 == 0) {
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
        cout << "���� " << NX+1 << "x" << NX+1 << " �p�⧹��\n" << endl;
    }
    cout << "�Ҧ��p�⧹���I" << endl;
    return 0;
}
