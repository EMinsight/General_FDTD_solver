//---------------------------------------------------------------------------
#ifndef TFDTD_Hybrid_mod_with_MOMH
#define TFDTD_Hybrid_mod_with_MOMH
//---------------------------------------------------------------------------
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <complex>
#include "TFunctions.h"

// using namespace std;
typedef struct struct_Hybrid_source_info
{
    double f0;
    double ds;
    int SPX;
    int SPY;
    int SPZ;
    int Datalen;

    char  Fname_Xminus[256];
    char  Fname_Xplus[256];
    char  Fname_Yminus[256];
    char  Fname_Yplus[256];
    char  Fname_Zminus[256];
    char  Fname_Zplus[256];
} struct_Hybrid_source_info;

typedef struct struct_Surface_Current
{
	int PX;
	int PY;
	int PZ;

    double  abs_Jx;
    double  abs_Jy;
    double  abs_Jz;
    double  abs_Mx;
    double  abs_My;
    double  abs_Mz;

    double  phase_Jx;
    double  phase_Jy;
    double  phase_Jz;
    double  phase_Mx;
    double  phase_My;
    double  phase_Mz;

} struct_Surface_Current;

class TFDTD_Hybrid_mod_with_MOM
{
	public:
        TFDTD_Hybrid_mod_with_MOM(double  ***Ex, double  ***Ey, double  ***Ez, double  ***Hx, double  ***Hy, double  ***Hz);
        ~TFDTD_Hybrid_mod_with_MOM();

	public:
         void insert_parameters(struct_Hybrid_source_info* pSsource);
         void set_surface_current_data();
         void set_Timestep(int ss){this->ss = ss;};
         void set_dt(double  dt){this->dt = dt;};
         void Apply_sourceE();
         void Apply_sourceH();         
         void Set_In_zeroE();
         void Set_In_zeroH();
         void Calculate_Incident_Fiend();
	private:
        double  ***Ex, ***Ey, ***Ez;
        double  ***Hx, ***Hy, ***Hz;
        int Datalen;
        int ss;
        int PML;
        int SPX;
        int SPY;
        int SPZ;
        double  f0;
        double  dt;
        double pi;
        double ds;
        double  mu0;
        double  eps0;
        double alpha;
        struct_Surface_Current *px_plus, *px_minus, *py_plus, *py_minus, *pz_plus, *pz_minus;
        FILE *pfXplus, *pfXminus, *pfYplus, *pfYminus, *pfZplus, *pfZminus;

        int ii_min, ii_max, jj_min, jj_max, kk_min, kk_max;
        int ii, jj, kk;
        int xi, xt, yi, yt, zi, zt;
        int B_Tx_ini, B_Tx_ter;
        int B_Ty_ini, B_Ty_ter;
        int B_Tz_ini, B_Tz_ter;

        double  **Ex_xyN_correction;
        double  **Ex_xyP_correction;
        double  **Ex_xzN_correction;
        double  **Ex_xzP_correction;

        double  **Ey_xyN_correction;
        double  **Ey_xyP_correction;
        double  **Ey_yzN_correction;
        double  **Ey_yzP_correction;

        double  **Ez_xzN_correction;
        double  **Ez_xzP_correction;
        double  **Ez_yzN_correction;
        double  **Ez_yzP_correction;

        double  **Hx_xyN_correction;
        double  **Hx_xyP_correction;
        double  **Hx_xzN_correction;
        double  **Hx_xzP_correction;

        double  **Hy_xyN_correction;
        double  **Hy_xyP_correction;
        double  **Hy_yzN_correction;
        double  **Hy_yzP_correction;

        double  **Hz_xzN_correction;
        double  **Hz_xzP_correction;
        double  **Hz_yzN_correction;
        double  **Hz_yzP_correction;

    private:
        double  sinrE(double  phase);
        double  sinrH(double  phase);

        void Initial_condition();

        void Apply_xE();
        void Apply_yE();
        void Apply_zE();

        void Apply_xH();
        void Apply_yH();
        void Apply_zH();


        double  get_phase(double re, double im);
        double  get_abs(double re, double im);
        int round(double nominator, double denominator);
};

#endif
