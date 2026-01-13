//---------------------------------------------------------------------------
#ifndef TSaveSnapshotH
#define TSaveSnapshotH
//---------------------------------------------------------------------------
#include <math.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
//#include <complex>
#include <direct.h>

#include "TFunctions.h"
#include "TFourierTransform.h"
// using namespace std;

typedef struct struct_Obs_Composit
{
	int Oz;
	double Roz;
	double *SumEx;
	double *SumEy;
} struct_Obs_Composit;

typedef struct struct_Obs_Point
{
	int Ox, Oy, Oz;
	double Rox, Roy, Roz;
	double *Ex_save_array;
	double *Ey_save_array;
	double *Ez_save_array;

	double *Hx_save_array;
	double *Hy_save_array;
	double *Hz_save_array;
}  struct_Obs_Point;

typedef struct struct_Obs_Plane
{
	int s_plane; // 1: xy, 2: yz, 3: zx
	double r_th_real;
	int r_th;	// n th ���
}  struct_Obs_Plane;

class TSaveSnapshot
{
	public:
		__fastcall TSaveSnapshot(double ***Ex, double ***Ey, double ***Ez, double ***Hx, double ***Hy, double ***Hz);
		__fastcall TSaveSnapshot() { strcpy(JsonFileName, "./IEM_solver_input_format.json"); };
		__fastcall~TSaveSnapshot();
		void get_mesh_length(int txsize, int tysize, int tzsize);
		void set_spatial_and_dt(double dx, double dy, double dz, double dt);
		void Field_3D_All();
		void Field_Obs_Plane(int ss);
		void Field_Obs_Points();
		void Field_XY_2D_Each(int kk, int ss);
		void Field_YZ_2D_Each(int ii, int ss);
		void Field_XZ_2D_Each(int jj, int ss);
		void Read_save_info(); // ���Ϸκ��� ���� �б�.
		//void set_TimeDiff_with_iterationNum(double time, int Nststep_max, int Nststep_iter);
		void set_TimeDiff_NstepMax(double dt, int Nststep_max);
		void Set_Current_iter(int ss) { Nststep_iter = ss; };
		void Set_max_iter();
		void Save_Composition();
    public:
        void get_Obs_point(int *xx1, int *yy1, int *zz1, int *xx2, int *yy2, int *zz2);
		
	protected:
		int isChild;
		char JsonFileName[256];
		double dx, dy, dz, dt;
		double ***EOx, ***EOy, ***EOz;
		double ***HOx, ***HOy, ***HOz;

		int txsize, tysize, tzsize;
		double *time_obs, *time_comp;
		int Nststep_Max;
		int Nststep_iter;
		int Nststep_save;
		//
		FILE *pf_Obs, *pf_Comp, *pf_Plane;
		// ������ ���� �̸�
		char **Obs_save_name;
		char **Plane_save_name;
		char **Comp_save_name;
		char format_size[256];
	protected:
		// Save fields
		void XY_plane_plot(int kk);
		void YZ_plane_plot(int ii);
		void XZ_plane_plot(int jj);
		// Number of Obs points
		int no_points;
		int no_composit;
		struct_Obs_Composit *pSOC;
		struct_Obs_Point *sop;
		struct_Obs_Plane *sopd;

		// Number of Planes
		int no_planes;

		// Field value
		double ExV, EyV, EzV;
		double HxV, HyV, HzV;
		double sumEx;
		double sumEy;
		// Save Field Name
		FILE *pf1, *pfField;

		// ObsPlane
		// E intensity save
		double ***E_value;
		double **E_value_2D;

		// optical chirality
	public:
		// void store_2D_field_frequency();
		// void save_optical_chirality(int ss);
		// // void Calculate_optical_chirality();
		// void Initiate_optical_option();
		
	public:
		void set_Media(int ***Media);
	protected:
		int ***Media;
		// 2���� ����, optical chirality�� ������ ��� ���� : 3���� ���� �Ҵ�.
		// struct_optical_chirality ***pOptical_plane;
		// void Assign_optical_chirality_plane();
		// 
		// void Transform_xy_plane_dft(int xx, int yy);
		// void Transform_chiral_surface_dft(int xx, int yy);
	protected:
		char Fname_Ex1_re[256], Fname_Ey1_re[256], Fname_Ez1_re[256], Fname_Hx1_re[256], Fname_Hy1_re[256], Fname_Hz1_re[256];
		char Fname_Ex2_re[256], Fname_Ey2_re[256], Fname_Ez2_re[256], Fname_Hx2_re[256], Fname_Hy2_re[256], Fname_Hz2_re[256];
		char Fname_Ex1_im[256], Fname_Ey1_im[256], Fname_Ez1_im[256], Fname_Hx1_im[256], Fname_Hy1_im[256], Fname_Hz1_im[256];
		char Fname_Ex2_im[256], Fname_Ey2_im[256], Fname_Ez2_im[256], Fname_Hx2_im[256], Fname_Hy2_im[256], Fname_Hz2_im[256];

		FILE *pF_Ex1_re, *pF_Ey1_re, *pF_Ez1_re, *pF_Hx1_re, *pF_Hy1_re, *pF_Hz1_re;
		FILE *pF_Ex2_re, *pF_Ey2_re, *pF_Ez2_re, *pF_Hx2_re, *pF_Hy2_re, *pF_Hz2_re;
		FILE *pF_Ex1_im, *pF_Ey1_im, *pF_Ez1_im, *pF_Hx1_im, *pF_Hy1_im, *pF_Hz1_im;
		FILE *pF_Ex2_im, *pF_Ey2_im, *pF_Ez2_im, *pF_Hx2_im, *pF_Hy2_im, *pF_Hz2_im;

		int **Extract_Optical_Chirality_plane;
		int sampling;
		TFourierTransform *pFourier;
		double Ex1_of, Ey1_of, Ez1_of, Hx1_of, Hy1_of, Hz1_of;
		double Ex2_of, Ey2_of, Ez2_of, Hx2_of, Hy2_of, Hz2_of;
		double *pSave_Temp_re, *pSave_Temp_im;

	protected:
		int ss_loc;

};
#endif
