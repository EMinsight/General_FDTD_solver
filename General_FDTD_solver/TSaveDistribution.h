#pragma once
#include "TFunctions.h"
// f모드, t모드 둘다 해야 함.
// #include "TFourierTransform.h"
#include "TSaveFDTD_P_pros_Comm.h"
#include "TEnvelop_Extract.h"

typedef struct struct_Obs_Plane
{
	int s_plane; // 1: xy, 2: yz, 3: zx
	double r_th_real;
	int r_th;	// n th 평면
	bool is_sopdf = false;
}  struct_Obs_Plane;

typedef struct struct_Obs_Plane_frequency
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
}struct_Obs_Plane_frequency;

class TSaveDistribution : public TSaveFDTD_P_pros_Comm
{
public:
	TSaveDistribution();
	~TSaveDistribution() {};
	void Save(int Nststep_iter);
	void Prepare2save();

private:
	FILE* pf_Plane;
	// Document Jdocument;
	// char JsonFileName[1024];
	char PulseName[512];
	char** Plane_save_name;
	char pFileName[1024];
	char pFileName_tmp[1024];

	int no_planes;

	int total_no_planes;
	int no_snap;
	int snap_iter;
	int snap_count;
	struct_Obs_Plane* sopd;

	// TFourierTransform* pFourier;
	struct_Obs_Plane_frequency* sopds;
	string* Obs_save_name;

	char Domain[16];
	double*** E_value;
	double** E_value_2D;
	double ExV, EyV, EzV, HxV, HyV, HzV;
	void Read_save_info();
	void Read_save_mod();
	void XY_plane_plot(int kk, char* Fname);
	void YZ_plane_plot(int ii, char* Fname);
	void XZ_plane_plot(int jj, char* Fname);
	void Save_POI(int Nststep_iter);
	void Field_XY_2D_Each(int kk, int ss, char* Fname);
	void Field_YZ_2D_Each(int ii, int ss, char* Fname);
	void Field_XZ_2D_Each(int jj, int ss, char* Fname);

	// frequency
	FILE* pFEx_Amp, * pFEy_Amp, * pFEz_Amp;
	FILE* pFEx_Phase, * pFEy_Phase, * pFEz_Phase;
	double* pSave_ExAmp, * pSave_ExPhase, * pSave_EyAmp, * pSave_EyPhase, * pSave_EzAmp, * pSave_EzPhase;
	double* pSave_HxAmp, * pSave_HxPhase, * pSave_HyAmp, * pSave_HyPhase, * pSave_HzAmp, * pSave_HzPhase;
	struct_no_planes_envelop* psEnvelop_Material;
	int findex;
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

	// 
	void assign_electric_flux_D();
};