#pragma once

#include "TFunctions.h"
#include "TFourierTransform.h"
#include "TSaveFDTD_P_pros_Comm.h"
//using namespace std;

typedef struct struct_Obs_Composit
{
	int Oz;
	double Roz;
	double* SumEx;
	double* SumEy;
	double* SumEz;
	double* SumHx;
	double* SumHy;
	double* SumHz;
} struct_Obs_Composit;

class TSaveComposit:public TSaveFDTD_P_pros_Comm
{
public:
	TSaveComposit();
	~TSaveComposit(){};
	void Save(int Nststep_iter);
	void Prepare2save();
private:
	FILE* pf_Obs;
	int no_composit;
	struct_Obs_Composit* pSOC;
	double* time_comp;
	double sumEx, sumEy, sumEz, sumHx, sumHy, sumHz;
	void Read_save_info();
	char** Comp_save_name;
	void Assing_save_length();
	
};

