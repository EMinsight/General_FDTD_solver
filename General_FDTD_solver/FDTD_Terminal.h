#pragma once

#define WIN32_LEAN_AND_MEAN   // GDI, RPC, Shell, WinSock, COM 등 불필요한 부분 제외
#define NOMINMAX  
#include <Windows.h>


#include "TFDTD_Generic_PBC.h"
#include "TFDTD_GenericAlgorithm.h"
#include "TFDTDparameters.h"
#include "TFunctions.h"
#include "TPostprocessing.h"
#include "TSource_Excitation.h"
#include "Toptical_chirality.h"
#include "FDTD_Boundary_Condition.h"

#include <memory.h>
#include <direct.h>
#include <time.h>

class Tmain
{
public:
	Tmain();
	~Tmain();
protected:
	double m_pi;
private: // Structure parameters
	double dt;
	double* dx, * dy, * dz;
	int PMLX, PMLY, PMLZ;
	int txsize, tysize, tzsize;

	double eps0;
	char JsonFileName[256];
private: // Field Components
	unsigned long ss, NSTEPS;
	int i, j, k;
	int total_pol;
	double*** Ex, *** Ey, *** Ez;
	double**** Jx, **** Jy, **** Jz;
	double*** Hx, *** Hy, *** Hz;
	int*** media;
	int*** media_Ex;
	int*** media_Ey;
	int*** media_Ez;
	double*** E_value;
	int index;
	int snap_index;
	int nosnap; // number of snapshot
	int snapshot_iter; // snapshot period

private: //Source parameters
	int opt_source;
	double f0;
	double BW;
	double Pulse;
	int number_of_source;
	int SourceType;
	int IsPlanewave;
	int selContents;
	double**** TestPmlField;

	double* Phase;
	double* PW_Theta;
	double* PW_Fhi;
	// 90 0
	double* PW_Polar;

public:		// User declarations
	void TFDTD_main_routine();

private:
	bool PBCmod;
	// step 1
	void Source_infomation_upload();
	// step 2
	void Generic_FDTD_mode();
	void Generic_PBC_FDTD_mode();
	// step 3
	void Initiallized_Fields();
	void Post_processing_info();
	void Time_marching();
	void Time_marching_PBC();

private: // input parameters from manager
	void setSourceShape(int os) { opt_source = os; };
	void setCenterFrequency(double f0) { this->f0 = f0; };
	void setNosnapshot(int nosnap) { this->nosnap = nosnap; };

private:
	struct_Source_type* sst;
	struct_Field_carrier sFc;
	Material_properties *pMp;
	TSource_Excitation* Source1;
	TPostprocessing* pPost;
	TFDTD_GenericAlgorithm* pGeneric;
	TFDTD_Generic_PBC* pGeneric_PBC;
	Toptical_chirality* pOptical;
	
	// PML
	FDTD_Boundary_Condition* pBC;
	int dum1, dum2, dum3, planewave_pos;

	// omp
	int active_thread_num;
};