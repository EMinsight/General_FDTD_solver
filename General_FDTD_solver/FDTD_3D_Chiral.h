#pragma once

#include "FDTD_Boundary_Condition.h"
#include "TFDTD_Generic_PBC.h"
#include "TFDTD_GenericAlgorithm.h"
#include "TFDTD_YeeAlgorithm.h"
#include "TFDTDparameters.h"
#include "TFunctions.h"
#include "TSaveSnapshot.h"
#include "TSaveSnapshot_expend.h"
#include "TSource_Excitation.h"
#include "Toptical_chirality.h"

#include <Windows.h>
#include <memory.h>
#include <direct.h>
#include <time.h>

class TChiral
{
public:
	TChiral();
	~TChiral();
protected:
	double m_pi;
private: // Structure parameters
	double dt;
	double *dx, *dy, *dz;
	int PMLX, PMLY, PMLZ;
	int txsize, tysize, tzsize;
	int PbufferX, PbufferY, PbufferZ;
	double eps0;
	char JsonFileName[256];
private: // Field Components
	unsigned long ss, NSTEPS;
	int i, j, k;
	int total_pol;
	double ***Ex, ***Ey, ***Ez;
	double ****Jx, ****Jy, ****Jz;
	double ***Hx, ***Hy, ***Hz;
	int*** media;
	int*** media_Ex;
	int*** media_Ey;
	int*** media_Ez;
	double ***E_value;
	int index;
	int snap_index;
	int nosnap; // number of snapshot
	int snapshot_iter; // snapshot period

private: //Source parameters
	int opt_source;
	double f0;
	double BW;
	double Peroid;
	double Pulse;
	int number_of_source;
	int SourceType;
	int IsPlanewave;
	int selContents;
	double ****TestPmlField;

	double* Phase;
	double *PW_Theta;
	double *PW_Fhi;
		// 90 0
	double *PW_Polar;


public:		// User declarations
	void TFDTD_main_routine();
private:
	int mod;
	int isCFS;
	void Source_infomation_upload();
	void Initiallized_Fields();
	void Field_data_assign_Generic();
	void Run_Initiate();
	void Run_Chiral();

private: // input parameters from manager
	void setSourceShape(int os) { opt_source = os; };
	void setCenterFrequency(double f0) { this->f0 = f0; };
	void setNosnapshot(int nosnap) { this->nosnap = nosnap; };

private:
	struct_Field_carrier sFc;
	struct_Source_type* sst;
	TFDTD_YeeAlgorithm *pFDTD1;
	TSource_Excitation *Source1;
	TSaveSnapshot_expend *pSnap;
	TFDTD_GenericAlgorithm *pGeneric;
	TFDTD_Generic_PBC *pChiral;
	Toptical_chirality *pOptical;
	FDTD_Boundary_Condition* pBC;

	int dum1, dum2, dum3, planewave_pos;
};

