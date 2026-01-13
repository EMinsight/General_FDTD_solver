#pragma once
#include "TSaveFDTD_P_pros_Comm.h"
#include "TFourierTransform.h"
#include "TEnvelop_Extract.h"

typedef struct struct_optical_chirality
{
	char Name[16];
	char Position[16]; // Transmission or Reflection
	double **Ex_save_array_t1, **Ex_save_array_t2;
	double **Ey_save_array_t1, **Ey_save_array_t2;
	double **Ez_save_array_t1, **Ez_save_array_t2;
	double **Hx_save_array_t1, **Hx_save_array_t2;
	double **Hy_save_array_t1, **Hy_save_array_t2;
	double **Hz_save_array_t1, **Hz_save_array_t2;
	double **pOptical_chirality;
} struct_optical_chirality;


class Toptical_chirality : public TSaveFDTD_P_pros_Comm
{
public:
	Toptical_chirality();
	~Toptical_chirality() {};
	void Save(int Nststep_iter);
	void Prepare2save();
	// optical chirality : 공용
private:
	double eps0;
	double mu0;
	double pi;
	double c0;

// optical chirality : Time
private:
	int i_PBC, j_PBC, aa1m1, aa1m2;
	int ActualNum;
	char matName[16];
	bool flag;
	double pivot_min;
	double pivot_max;

	FILE *pF_CC;
	char Fname_CC[256];
	void optical_chirality_caclulate();
	double curl_E(int ii, int jj, int kk, char comp);
	double curl_H(int ii, int jj, int kk, char comp);
	double E_dot_p(int ii, int jj, int kk, char comp);
	double H_dot_p(int ii, int jj, int kk, char comp);
	
	double T0;
	double Amp, Phase;
	int j_offset;
	void Amp_Phase_Tunning();
	struct_no_planes_envelop* psEnvelop_Material;
	int*** pExtract_Optical_Chirality_plane;
	double ***pOptical_chirality;
	double ***pOptical_chirality_frequency;
	char **Plane_save_name;

	int Layer;
	int sampling;
	double f0  ;
	double sfre;
	double tfre;
	double df;
	int tt1, tt2;
	int source_pos;

	TFourierTransform *pFourier;
	FILE *pFname_optical_chirality;
	char Fname_optical_chirality_path[256];
	char Fname_optical_chirality[256];

	// 2차원 공간, optical chirality를 추출할 평면 개수 : 3차원 공간 할당.
	struct_optical_chirality *pOptical_plane;
	void MemoryAssign_for_calculate_optical_chirality();
	void store_2D_field_frequency(int Nststep_iter);
	// void DFT_chiral_surface_wide(int xx, int yy, int Nststep_iter);
	void Single_Freq_Amp_Phase1(int xx, int yy);
	void Single_Freq_Amp_Phase2(int xx, int yy);

	double *pSave_Temp_re, *pSave_Temp_im;
	
	void save_optical_chirality_Frequency();
	void Initiate_optical_option_Frequency();
	// void Calculate_optical_chirality_Frequency();
	void Exctract_chiral_surface();
};

