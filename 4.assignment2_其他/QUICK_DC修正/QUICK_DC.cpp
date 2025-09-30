//�Q�ΪŶ��G�����QUICK��(DC�ץ�)�D�ѤG��í���X����y���D
//�S�I:�ޤJ(DC�ץ�)
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
const int NX = 81;
const int NY = 81;
const int n = NX * NY;
double a[n+2][n+2], b[n+2]; //a�G���x�}�����������Y�Ưx�};b�@���x�}�������k�����~�����A��{�դ��D���ʶ�
double x[n+2], x_old[n+2] , x_old2[n+2] , x_nuse[n+2] , T[NX][NY] ; //x_nuse�����ϥ��P���޳N���s���N�ū׳�
const double Pelect = 1000 ;
const double  Gamma = (sqrt(2)/Pelect) ;//���X���Y�� alpha
const double dx = 1.0/double(NX-1) ;  // �ץ��G���涡�Z����1/(NX-1)
const double dy = 1.0/double(NY-1) ;  // �ץ��G���涡�Z����1/(NY-1)
const double u = 1.0;  // x��V�t��
const double v = 1.0;  // y��V�t��
const double T_l = 1.0 ;
const double T_r = 0.0 ;
const double T_b = 0.0 ;
const double T_u = 1.0 ;
const double lamda = 0.3 ;
int G, max_G = 10000000000;
double maxerror; 
const float tolerance = 1e-10;
bool steadystate;

//��l�ƫY�Ưx�}�P�쥻QUICK����P
//��l�Ưx�}//if(G == 0)
void initial(double a[][n+2], double b[], int n) {
    // ��l�ƩҦ�������0
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

//����x�}���� if(G>=0)
void DC_modified(double a[][n+2], double b[], int n) {
    ///////////////////////////////////////////////////////////////
    //�w�q�Ÿ�:(�F��)
    const double D_w = Gamma ; //for(i = 1 : NX-1)//�X���Ÿ� 
    const double D_wb = 2*Gamma ; //for(i = 0)
    const double D_e  = Gamma ; //for(i = 0 ; Nx-2)
    const double D_eb = 2*Gamma ;//for(i = Nx-1)
    const double F_w = dy*u ; //for(i = 0 : NX-1)//��y�Ÿ� 
    const double F_e = dy*u ; //for(i = 0 : NX-1)
    const double a_W = D_w + F_w ; //for(i = 1 : NX-1)//�u�իY�� 
    const double a_Wb = 0.0; //for(i = 0)
    const double a_E = D_e ; //for(i = 0 : NX-2)
    const double a_Eb = 0.0; //for(i = NX-1)
    const double S_ul = D_wb + F_w ; //for(i = 0)
    const double S_ur = D_eb - F_e ;
    
    //�w�q����ץ����k�~����S_{uw}^{DC}
    vector<double> S_uwdc(n+2);
    S_uwdc[0] = 0.0 ;
    S_uwdc[n+1] = 0.0 ;  // �ץ��G���ӬOn+1�Ӥ��On+2
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
    
    //�w�q����ץ����k�~����S_{ue}^{DC} 
    vector<double> S_uedc(n+2);
    S_uedc[0] = 0.0 ;
    S_uedc[n+1] = 0.0 ;  // �ץ��G���ӬOn+1�Ӥ��On+2
    for(int j = 1 ; j <= NY ; j++){
        for(int i = 1 ; i <= NX ; i++){
            if(i== 1){
                S_uedc[NX*(j-1)+i] = F_e*((3.0/8.0)*x_old[NX*(j-1)+i+1] - (1.0/8.0)*x_old[NX*(j-1)+i] - (2.0/8.0)*T_l) ;
            }else if(i == NX){
                S_uedc[NX*(j-1)+i] = 0.0 ;
            }else{
                S_uedc[NX*(j-1)+i] = F_e*((3.0/8.0)*x_old[NX*(j-1)+i+1] - (2.0/8.0)*x_old[NX*(j-1)+i] - (1.0/8.0)*x_old[NX*(j-1)+i-1]) ;  // �ץ��GF_w�אּF_e
            }
        }
    }
    
    //�w�q�Ÿ�:(�n�_) 
    const double D_s = Gamma ; //for(j = 1 : NY-1)//�X���Ÿ� 
    const double D_sb = 2*Gamma ; //for(j = 0)
    const double D_n  = Gamma ; //for(j = 0 ; NY-2)
    const double D_nb = 2*Gamma ;//for(j = NY-1)
    const double F_s = dx*v ; //for(j = 0 : NY-1)//��y�Ÿ� 
    const double F_n = dx*v ; //for(j = 0 : NY-1)
    const double a_S = D_s + F_s ; //for(j = 1 : NY-1)//�u�իY�� 
    const double a_Sb = 0.0; //for(j = 0)
    const double a_N = D_n ; //for(j = 0 : NY-2)
    const double a_Nb = 0.0; //for(j = NY-1)
    const double S_ub = D_sb + F_s ; //for(j = 0)//���k�~���� 
    const double S_uu = D_nb - F_n;//for(j = NY) 
    
    //�w�q����ץ����k�~����S_{us}^{DC}
    vector<double> S_usdc(n+2);
    S_usdc[0] = 0.0 ;
    S_usdc[n+1] = 0.0 ;  // �ץ��G���ӬOn+1�Ӥ��On+2
    for(int j = 1 ; j <= NY ; j++){
        for(int i = 1 ; i <= NX ; i++){
            if(j == 1){
                S_usdc[NX*(j-1)+i] = 0.0 ;
            }else if(j == 2){
                S_usdc[NX*(j-1)+i] = F_s*((3.0/8.0)*x_old[NX*(j-1)+i] - (1.0/8.0)*x_old[NX*(j-2)+i] - (2.0/8.0)*T_b) ;
            }else{
                S_usdc[NX*(j-1)+i] = F_s*((3.0/8.0)*x_old[NX*(j-1)+i] - (2.0/8.0)*x_old[NX*(j-2)+i] - (1.0/8.0)*x_old[NX*(j-3)+i]) ;  // �ץ��GS_uwdc�אּS_usdc
            }
        }
    }
    
    //�w�q����ץ����k�~����S_{un}^{DC} 
    vector<double> S_undc(n+2);
    S_undc[0] = 0.0 ;
    S_undc[n+1] = 0.0 ;  // �ץ��G���ӬOn+1�Ӥ��On+2
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
    
    // ���I�B�z
    //�Ĥ@���IT_{0,0} -> a[1][1] (���U��)
    a[1][1] =   (a_Wb + a_E + F_e) + S_ul - F_w + (a_Sb + a_N + F_n) + S_ub - F_s ;
    a[1][2] =    -a_E ;  // T_{1,0}  
    a[1][1+NX] = -a_N ;  // T_{0,1}
    b[1] =       -S_uedc[1] + S_uwdc[1] -S_undc[1] + S_usdc[1] + S_ul * T_l + S_ub * T_b ;  // �����=1.0

    //�ĤT���IT_{NX-1,0} -> a[NX][NX] (�k�U��)
    a[NX][NX] =    (a_W + a_Eb - F_w) + S_ur + F_e + (a_Sb + a_N + F_n) + S_ub - F_s;
    a[NX][NX-1] =  -a_W ;  // T_{NX-2,0}
    a[NX][2*NX] =  -a_N ;  // T_{NX-1,1}
    b[NX] =        -S_uedc[NX] + S_uwdc[NX] -S_undc[NX] + S_usdc[NX] + S_ur * T_r + S_ub * T_b;  // �k��ɩM�U��ɳ�=0.0

    //�ĤC���IT_{0,NY-1} -> a[n-NX+1][n-NX+1] (���W��)
    a[n-NX+1][n-NX+1] =   (a_Wb + a_E + F_e) + S_ul - F_w + (a_S + a_Nb - F_s) + S_uu + F_n ;
    a[n-NX+1][n-NX+2] =    -a_E ;  // T_{1,NY-1}
    a[n-NX+1][n-2*NX+1] =  -a_S ;  // T_{0,NY-2}
    b[n-NX+1] =          -S_uedc[n-NX+1] + S_uwdc[n-NX+1] -S_undc[n-NX+1] + S_usdc[n-NX+1] + S_ul * T_l + S_uu * T_u ;  // ����ɩM�W��ɳ�=1.0

    //�ĤE���IT_{NX-1,NY-1} -> a[n][n] (�k�W��)
    a[n][n] =      (a_W + a_Eb - F_w) + S_ur + F_e + (a_S + a_Nb - F_s) + S_uu + F_n  ;
    a[n][n-1] =    -a_W ;  // T_{NX-2,NY-1} 
    a[n][n-NX] =   -a_S ;  // T_{NX-1,NY-2}
    b[n] =        -S_uedc[n] + S_uwdc[n] -S_undc[n] + S_usdc[n] + S_ur * T_r + S_uu * T_u ;   

    /////////////////////////////////////////////////////////////////////
    // ��ɳB�z
    // �U��� (�����I�~)
    for(int i = 2; i <= NX-1; i++) {
        a[i][i] =   (a_W + a_E - F_w + F_e)+(a_Sb + a_N + F_n) + S_ub - F_s ;  // �ץ��G���I���ӥ�a_W�Ӥ��Oa_Wb
        a[i][i-1] = - a_W;
        a[i][i+1] = - a_E;
        a[i][i+NX] = - a_N;
        b[i] =  -S_uedc[i] + S_uwdc[i] -S_undc[i] + S_usdc[i]  + S_ub * T_b ;  // �ץ��G�������ӬOi�Ӥ��ONX * (i-1) + 1
    }
    
    // �W��� (�����I�~)
    for(int i = n-NX+2; i <= n-1 ; i++) {
        a[i][i] =  (a_W + a_E - F_w + F_e) +(a_S + a_Nb - F_s) + S_uu + F_n ;  // �ץ��G���I���ӥ�a_W�Ӥ��Oa_Wb
        a[i][i-1] = - a_W;
        a[i][i+1] = - a_E;
        a[i][i-NX] =  -a_S;
        b[i] = -S_uedc[i] + S_uwdc[i] -S_undc[i] + S_usdc[i]  + S_uu * T_u ;  // �ץ��G�������ӬOi
    }
    
    // ����� (�����I�~)(�W�U���I) 
    for(int i = 2 ; i <= NY-1 ; i++) {
        int index = NX * (i-1) + 1;
        a[index][index] = (a_Wb + a_E + F_e)+(a_S + a_N - F_s + F_n) + S_ul - F_w ;
        a[index][index+1]  = -a_E ;//�F�I 
        a[index][index-NX] = -a_S;  // �n�I
        a[index][index+NX] = -a_N;  // �_�I
        b[index] = -S_uedc[index] + S_uwdc[index] -S_undc[index] + S_usdc[index]  + S_ul * T_l ;
    }
    
    // �k��� (�����I�~)(�W�U���I) 
    for(int i = 2; i <= NY-1; i++) {
        int index = NX * i;
        a[index][index] =   (a_W + a_Eb - F_w) + (a_S + a_N - F_s + F_n) + S_ur + F_e;
        a[index][index-1] =  - a_W ;//���I  
        a[index][index-NX] = -a_S;  // �n�I
        a[index][index+NX] = -a_N;  // �_�I
        b[index] = -S_uedc[index] + S_uwdc[index] -S_undc[index] + S_usdc[index] + S_ur * T_r ;
    }
    
    // ���I
    for(int j = 2; j <= NY-1 ; j++) {
        for(int i = 2; i <= NX-1; i++) {
            int index = NX*(j-1)+i;
            a[index][index] = (a_W + a_E - F_w + F_e) + (a_S + a_N - F_s + F_n)  ;
            a[index][index-1] =-a_W ;  // ���I
            a[index][index+1] =-a_E ;  // �F�I
            a[index][index-NX] =-a_S;  // �n�I
            a[index][index+NX] =-a_N;  // �_�I
            b[index] = -S_uedc[index] + S_uwdc[index] -S_undc[index] + S_usdc[index] ;  // ���I�L����
        }
    }
}

void giveTvalue(int m){
    if(m == 0){//�g�L0�����N,�n�i��Ĥ@�����N 
    // ���ƻs��e�Ѩ�x_old
    //�H�ν�ȸѨ�e�G�� //112
        for(int k = 1; k <= n; k++) {
            x_old2[k] = x[k];
            x_old[k] = x[k];
        }
    }else if (m == 1){//�g�L�@�����N�A�n�i��ĤG�����N 
        // ���ƻs��e�Ѩ�x_old
        for(int k = 1; k <= n; k++) {
            //�e�e�Ѥ��ʡC //123
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
    // �p��s����
    for(int k = 1; k <= n; k++) {
        double sum = 0;
        for(int p = 1; p <= n; p++) {
            if(p != k) {
                sum += a[k][p] * x_old[p];
            }
        }
        x[k] = ((b[k] - sum) / a[k][k]);
        //���ͷs���ū׳� 
        x_nuse[k] = x[k] ;//���ϥ��P���޳N���s���N�ū׳� 
        //under relaxtion(���P�����N) && over relation(�W�P�����N)
        // �p��̤j�~�t
        if(fabs(x_nuse[k] - x_old[k]) > fabs (x_old[k] - x_old2[k]) ){
            //�~��ĥΤ��P���޳N
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
    // �N�@�����ഫ���G���ū׳�
    for(int j = 0; j < NY; j++) {
        for(int i = 0; i < NX; i++) {
            T[i][j] = x[j*NX + i + 1]; //��Jx[] : 1~n 
        }
    }
    
    ostringstream name;
    name << "QUICK_DC" << setfill('0') << setw(6) << m << ".vtk";
    ofstream out(name.str().c_str());
    
    // VTK ����Y
    out << "# vtk DataFile Version 3.0\n";
    out << "QUICK_DC\n";
    out << "ASCII\n";
    out << "DATASET STRUCTURED_POINTS\n";
    out << "DIMENSIONS " << NX << " " << NY << " 1\n";
    out << "ORIGIN 0 0 0\n";
    out << "SPACING " << dx << " " << dy << " 1\n";  // �ϥι�ں��涡�Z
    out << "POINT_DATA " << NX * NY << "\n";
    
    // ��X�ū׳�
    out << "SCALARS Temperature double 1\n";
    out << "LOOKUP_TABLE default\n";
    for(int j = 0; j < NY; j++) {
        for(int i = 0; i < NX; i++) {
            out << scientific << setprecision(6) << T[i][j] << "\n";
        }
    }
    
    out.close();
    cout << "VTK ���w��X: " << name.str() << endl;
}

int main() {
    cout << "�{���X�}�l����...." << endl;
    cout << "����j�p: " << NX << " x " << NY << endl;
    cout << "���涡�Z: dx=" << dx << ", dy=" << dy << endl;
    cout << "��ɱ���: ����ɩM�W���=1.0, �k��ɩM�U���=0.0" << endl;
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
            cout << "���N���� = " << G << endl;
            if(G > 0) {
                cout << "�̤j�ܤƶqmax_change = " << scientific << maxerror << endl;
            }
            output(G);
        }
        
        if(G > 100 && maxerror < tolerance) {
            steadystate = true;
            cout << "�w�g�F�쭡�N���I�A�ū׳�����!!" << endl;
            break;
        }
    }
    
    if(!steadystate) {
        cout << "�F��̤j���N���ơA�����F��í�A!" << endl;
    }
    
    output(G);
    cout << "��l�j�p " << NX << "x" << NY << " �p�⧹��\n" << endl;
    return 0;
}
