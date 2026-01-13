#pragma once
#include "TFunctions.h"
#include "TFourierTransform.h"
#include "TSaveFDTD_P_pros_Comm.h"
#include "TEnvelop_Extract.h"

struct struct_Current_density;
struct struct_Current_Position;
struct struct_Obs_Plane_datasave;

typedef struct struct_Current_Position
{
	double ff;
	double ii, jj, kk;
	struct_Current_density* sCd;
}struct_Current_Position;

typedef struct struct_Current_density
{
	double** Hx1, ** Hy1, ** Hz1;
	double** Hx2, ** Hy2, ** Hz2;

	double** Jx_save_plane1, ** Jx_save_plane2;
	double** Jy_save_plane1, ** Jy_save_plane2;
	double** Jz_save_plane1, ** Jz_save_plane2;

	double** AmpJx, ** PhaseJx;
	double** AmpJy, ** PhaseJy;
	double** AmpJz, ** PhaseJz;
}struct_Current_density; // 2Â÷¿ø °î¸é, 

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

class TSaveCurrentDistribution: public TSaveFDTD_P_pros_Comm
{
public:
	TSaveCurrentDistribution();
	~TSaveCurrentDistribution() {};
	void Save(int Nststep_iter) {};
	void Prepare2save() {};

private:
	void save_Current_density();
	void save_nXh();
	void Read_save_info();


	// reference func.
	FILE* pFJx_Amp, * pFJy_Amp, * pFJz_Amp;
	FILE* pFJx_Phase, *pFJy_Phase, *pFJz_Phase;
	double* pSave_ExAmp, * pSave_ExPhase, * pSave_EyAmp, * pSave_EyPhase, * pSave_EzAmp, * pSave_EzPhase;
	double* pSave_HxAmp, * pSave_HxPhase, * pSave_HyAmp, * pSave_HyPhase, * pSave_HzAmp, * pSave_HzPhase;
	struct_no_planes_envelop* psEnvelop_Material;
	int pos;
	bool flag;
	int ActualNum;
	char matName[16];
	char Fname_save_field_Path[256];
	char Fname_save_field_Pathcpy_amp[256];
	char Fname_save_field_Pathcpy_phase[256];
	char Fname_save_field_Name_temp[256];
	char Fname_save_field_Amp[512];
	char Fname_save_field_Phase[512];
	void Save_POI_frequency(int Nststep_iter);
	void MemoryAssign_for_FieldDistribution();
	void Calcuate_Amp_Phase();
	void Calcuate_Amp_PhaseXY(int pp);
	void Calcuate_Amp_PhaseYZ(int pp);
	void Calcuate_Amp_PhaseXZ(int pp);
	void store_snapshot_Time(int Nststep_iter);
	void SaveFieldXY(int pp, int t1t2);
	void SaveFieldYZ(int pp, int t1t2);
	void SaveFieldXZ(int pp, int t1t2);
	void save_FieldDistribution();
	void save_FieldDistribution_assign_file_pointer(int nn);

};

