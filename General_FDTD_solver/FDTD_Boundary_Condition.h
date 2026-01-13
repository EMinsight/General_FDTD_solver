#pragma once

#include "TFDTD_2PML.h"
#include "TFDTD_CFS_PML.h"
#include "TFDTD_PBC_Boundary_normal.h"
#include "TFDTD_2PML_PBC.h"
#include "TFDTD_PBC_Screwed.h"
class FDTD_Boundary_Condition
{
	public:
		FDTD_Boundary_Condition() {};
		FDTD_Boundary_Condition(struct_Field_carrier* pFc);
		~FDTD_Boundary_Condition();

	public:
		// void Set_EHJ_Fields(struct_Field_carrier* pFc);
		void Initicate();
		bool Get_BC_mode() { return PBC_flag; }
		void Boundary_Condition_PMLF();
		void Boundary_Condition_PBC();
		void Boundary_Condition_PMLG();
		int get_PMLX() { return PMLX; };
		int get_PMLY() { return PMLY; };
		int get_PMLZ() { return PMLZ; };

	private:
		void Read_Json();
		void Set_Boundary_Conditions_parameters();

		int*** Media_Ex, *** Media_Ey, *** Media_Ez;
		struct_Field_carrier *pFc;
		double*** Ex, *** Ey, *** Ez, *** Hx, *** Hy, *** Hz;
		int PMLX, PMLY, PMLZ;
		bool PBCX, PBCY, PBCZ;
		
		bool PML_flag;
		bool PBC_flag;
		bool CFS_flag;
		char JsonFileName[512];
		char nx[6], px[6], ny[6], py[6], nz[6], pz[6];

	private:
		TFDTD_2PML* pPML;
		TFDTD_CFS_PML* pCFSPML;
		TFDTD_2PML_PBC* pPMLPBC;
		TFDTD_PBC_Screwed* pPBCScrewed;
};