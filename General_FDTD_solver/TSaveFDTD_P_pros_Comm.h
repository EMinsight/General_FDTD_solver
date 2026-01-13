#pragma once
#include "TFunctions.h"

typedef struct struct_no_planes_envelop
{
	int*** pEnvelop;
	char Name[16];
	char Position[16]; // Transmission or Reflection
	//	struct_no_planes_envelop* plist;
}struct_no_planes_envelop;

// postprocessing에서 공통으로 사용되는 함수나, 변수를 추가.
class TSaveFDTD_P_pros_Comm
{
public:
	TSaveFDTD_P_pros_Comm() {};
	~TSaveFDTD_P_pros_Comm() {};
public:
	virtual void Save() {};
	void set_postprocessing_info(struct_Field_carrier* pFc, Material_properties* pMp);
protected:
	struct_Field_carrier* pFc;
	Document Jdocument;
	char JsonFileName[1024];
	int txsize, tysize, tzsize;
	double pi;
	double eps_0;
	double f0, BW;
	double dt;
	double T0;
	int T0_discrete, T0_stepSize;

	double* dx, * dy, * dz;
	double*** EOx, *** EOy, *** EOz;
	double*** HOx, *** HOy, *** HOz;
	
	double ExV, EyV, EzV, HxV, HyV, HzV, JxV, JyV, JzV;
	Material_properties* pMp;
	int*** Media;
	int Nststep_Max, Nststep_save, ss_loc, * time_obs;

	// steady-state 
	int save_sta;
	int save_end;

	void Voxel_center_average(int xx, int yy, int zz, string who);
	void Voxel_face_average(int xx, int yy, int zz, string who);
	void set_TimeDiff_NstepMax(double dt, int Nststep_max);
	void set_steady_start(double percent);
	void set_steady_end(double percent);
	void set_T0(double T0);

	// Transform : Rotation
	int tt1, tt2, sampling;
	double sfre, tfre, df;
	double Amp, Phase;
	void makeRotation(double fhi, double theta, double psi);
	void Amp_Phase_Tunning();
	double** RotationZXZp;
};

