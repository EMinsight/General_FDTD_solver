#pragma once
#include "TSaveSnapshot.h"
#include "TFourierTransform.h"
#include "TFunctions.h"

 typedef struct struct_Current_density
 {
 	double Amp[3], Phase[3];
 	double J1[3], J2[3];
 
 	// f1
 	double Hx1, Hy1, Hz1;
 	double Hxdy1, Hxdz1;
 	double Hydx1, Hydz1;
 	double Hzdx1, Hzdy1;
 
 	//f2
 	double Hx2, Hy2, Hz2;
 	double Hxdy2, Hxdz2;
 	double Hydx2, Hydz2;
 	double Hzdx2, Hzdy2;
 }struct_Current_density; // 2차원 곡면, 

//typedef struct struct_Current_density
//{
//	double** Hx1, ** Hy1, ** Hz1;
//	double** Hx2, ** Hy2, ** Hz2;
//
//	double** Jx_save_plane1, ** Jx_save_plane2;
//	double** Jy_save_plane1, ** Jy_save_plane2;
//	double** Jz_save_plane1, ** Jz_save_plane2;
//
//	double** AmpJx, ** PhaseJx;
//	double** AmpJy, ** PhaseJy;
//	double** AmpJz, ** PhaseJz;
//}struct_Current_density; // 2차원 곡면, 

typedef struct struct_Obs_Plane_datasave
{
	double** Ex_save_plane1, ** Ex_save_plane2;
	double** Ey_save_plane1, ** Ey_save_plane2;
	double** Ez_save_plane1, ** Ez_save_plane2;
	double** AmpEx, ** PhaseEx;
	double** AmpEy, ** PhaseEy;
	double** AmpEz, ** PhaseEz;

	// double **Hx_save_plane1, **Hx_save_plane2;
	// double **Hy_save_plane1, **Hy_save_plane2;
	// double **Hz_save_plane1, **Hz_save_plane2;

}struct_Obs_Plane_datasave;

class TSaveSnapshot_expend : public TSaveSnapshot
{
	// 외부 선언시 사용할수 있는 함수.
public:
	// 전류분포, 필드분포 공통.
	TSaveSnapshot_expend(struct_Field_carrier* sFc, Material_properties *pMp);
	TSaveSnapshot_expend() { pi = acos(-1); };
	~TSaveSnapshot_expend() {};

public:
	// void set_Media(int*** Media_Ex, int*** Media_Ey, int*** Media_Ez);
	void set_Frequency_range(double f0, double sfre, double tfre, int sampling);
	void store_snapshot_Time(int Nststep_iter);
	void store_current_Time(int Nststep_iter);
	void set_steady_start(double steady_state_start);
	void set_T0(double T0);
	// 내부에서 돌아가는 부분.
protected:
	Material_properties* pMp;
	struct_Field_carrier* sFc;
	double  pi;
	double Ex_of, Ey_of, Ez_of, Hx_of, Hy_of, Hz_of;
	double* pSave_ExAmp, * pSave_ExPhase, * pSave_EyAmp, * pSave_EyPhase, * pSave_EzAmp, * pSave_EzPhase;
	int*** Media_Ex, *** Media_Ey, *** Media_Ez;
	double f0, sfre, tfre, df;
	int tt1, tt2;
	double T0;
	double Amp, Phase;
	int steady_state_start_discrete;
	int steady_state_step;
	int T0_discrete;
	int T0_stepSize;
	int j_offset;
	TFourierTransform* pFourier;

	void Amp_Phase_Tunning();
	void Voxel_center_average(int xx, int yy, int zz);

	//  저장 경로 폴더 생성
	void save_FieldDistribution_make_saveFolder(int kk);
	void save_FieldDistribution_assign_file_pointer(int ff);

	// 필드분포
	FILE* pFEx_Amp, * pFEx_Phase;
	FILE* pFEy_Amp, * pFEy_Phase;
	FILE* pFEz_Amp, * pFEz_Phase;
	char Fname_save_field_frequency_path[256];
	char Fname_save_field_frequency[256];

	char Fname_save_field_Path[256];
	char Fname_save_field_Pathcpy_amp[256];
	char Fname_save_field_Pathcpy_phase[256];
	char Fname_save_field_Name_temp[256];
	char Fname_save_field_Amp[256];
	char Fname_save_field_Phase[256];

	//switch 관련 분배
	void SaveFieldXY(int ff, int pp, int t1t2);
	void SaveFieldYZ(int ff, int pp, int t1t2);
	void SaveFieldXZ(int ff, int pp, int t1t2);
	struct_Obs_Plane_datasave** sopds;
	void MemoryAssign_for_FieldDistribution();
	void Calcuate_Amp_Phase();
	void Calcuate_Amp_PhaseXY(int ff, int pp);
	void Calcuate_Amp_PhaseYZ(int ff, int pp);
	void Calcuate_Amp_PhaseXZ(int ff, int pp);
	void save_FieldDistribution();

	// 전류분포 (곡면구조 금속이라 3차원 공간값을 따야함)
	// 현재는 xy 평면만 가능.

	struct_Current_density	**Jup, ** Jdown; // 카이랄 메타 물질의 윗면 아랫면에 대한 전류분포추출.
	// struct_Current_density** scd; // 주파수, 평면
	int*** pExtract_2D_plane;
	int** pExtract_2D_envelope;

	int ***pExtract_current_distribution_plane   ;
	int **pExtract_current_distribution_envelope;

	FILE* pFname_Current_density_Amp;
	FILE* pFname_Current_density_Phase;
	char Fname_Current_density_Path[256];
	char Fname_Current_density_Pathcpy[256];
	char Fname_Current_density_Name_temp[256];
	char Fname_Current_density_Amp[256];
	char Fname_Current_density_Phase[256];
	void Exctract_Metal_surface();
	// J1, J2 :  Amp, Phase 구하기 위한 지점
	void Single_Freq_Amp_Phase_J1(int xx, int yy);
	void Single_Freq_Amp_Phase_J2(int xx, int yy);
	void save_nXh();
	void save_Current_density();
};

