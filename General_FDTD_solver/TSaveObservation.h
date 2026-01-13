#pragma once
#include "TFunctions.h"
#include "TSaveFDTD_P_pros_Comm.h"
//using namespace std;

typedef struct struct_Obs_Point
{
	int Ox, Oy, Oz;
	double Rox, Roy, Roz;
	double* Ex_save_array;
	double* Ey_save_array;
	double* Ez_save_array;

	double* Hx_save_array;
	double* Hy_save_array;
	double* Hz_save_array;
}  struct_Obs_Point;

class TSaveObservation:public TSaveFDTD_P_pros_Comm
{
public:
	TSaveObservation();
	~TSaveObservation() {}
	void Save(int Nststep_iter);
	void Prepare2save();
private:
	FILE* pf_Obs;
	int no_points;
	struct_Obs_Point* sop;
	char** Obs_save_name;
	double* time_obs;
	void Read_save_info();
	void Assing_save_length();

};