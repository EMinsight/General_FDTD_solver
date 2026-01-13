#pragma once

#include "TFDTD_PML.h"
#include "TFDTD_2PML.h"
#include "TFDTD_2PML_PBC.h"
#include "TFDTD_CFS_PML.h"
#include "TFDTD_Chiral_Contents.h"
#include "TFDTD_GenericAlgorithm.h"
#include "TFDTD_PBC_Boundary_normal.h"
#include "TFDTD_YeeAlgorithm.h"
#include "TFDTDparameters.h"
#include "TFunctions.h"
#include "TSaveSnapshot.h"
#include "TSaveSnapshot_expend.h"
#include "TSource_Excitation.h"
#include "TFDTD_Boundary_TFSF.h"
#include "Toptical_chirality.h"
#include "TFDTD_Hybrid_mod_with_MOM.h"

#include <Windows.h>
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
	double *dx, *dy, *dz;
	int PMLX, PMLY, PMLZ;
	int txsize, tysize, tzsize;
	
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
<<<<<<<< HEAD:General_FDTD_solver/FDTD_Terminal.h
	// step 1
	void Assign_Generic_with_Get_mesh_parameters();
	// step 2
========
	void Source_infomation_upload_Hybrid();
>>>>>>>> Hybrid_2024.0:General_FDTD_solver/FDTD_3D_Hybrid.h
	void Initiallized_Fields();

	void Run_Initiate();
<<<<<<<< HEAD:General_FDTD_solver/FDTD_Terminal.h
	void Time_marching();
========
	// void Run_cavity_contents();
	void TFDTD_view_source();
	void Run_Hybrid();
>>>>>>>> Hybrid_2024.0:General_FDTD_solver/FDTD_3D_Hybrid.h

private: // input parameters from manager
	void setSourceShape(int os) { opt_source = os; };
	void setCenterFrequency(double f0) { this->f0 = f0; };
	void setNosnapshot(int nosnap) { this->nosnap = nosnap; };

private:
	struct_Source_type* sst;
	struct_Hybrid_source_info* pHsi;
	TFDTD_YeeAlgorithm *pFDTD1;
	TSource_Excitation *Source1;
	TSaveSnapshot_expend *pSnap;
	TFDTD_PBC_Boundary_normal *pPBC1;
	TFDTD_GenericAlgorithm *pGeneric;
	TFDTD_Chiral_Contents *pChiral;
	Toptical_chirality *pOptical;
	TFDTD_Hybrid_mod_with_MOM* pHybrid;
	// PML
	TFDTD_2PML* PML1;
	TFDTD_CFS_PML* pCFSPML;
	TFDTD_2PML_PBC* pPMLPBC;

	int dum1, dum2, dum3, planewave_pos;

};

