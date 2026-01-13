#include "FDTD_Boundary_Condition.h"

FDTD_Boundary_Condition::FDTD_Boundary_Condition(struct_Field_carrier* pFc)
{
	this->pFc = pFc;
	strcpy(JsonFileName, "./IEM_solver_input_format.json");
	PBC_flag = false;
	PML_flag = false;
	CFS_flag = false;
	Read_Json();
	Set_Boundary_Conditions_parameters();

}

FDTD_Boundary_Condition::~FDTD_Boundary_Condition() {}

void FDTD_Boundary_Condition::Initicate()
{
	if (PML_flag == true && PBC_flag == false)
	{
		pPML  = new TFDTD_2PML(pFc);
	}
	else if (CFS_flag == true && PBC_flag == false)
	{
		pCFSPML = new TFDTD_CFS_PML(pFc);
		pCFSPML->Set_Kappa_coefficient(0, 1);
	}
	else if (PML_flag == true && PBC_flag == true)
	{
		pPBCScrewed = new TFDTD_PBC_Screwed(pFc);
		pPBCScrewed->Set_Kappa_coefficient(0, 1);
	}
}

void FDTD_Boundary_Condition::Set_Boundary_Conditions_parameters() 
{
	if (strcmp(nx, "PML") == 0 && strcmp(px, "PML") == 0)
	{
		PML_flag = true;
		PMLX = 10;
		printf("X Boundary: PML\n");
	}
	else if (strcmp(nx, "CFSPML") == 0 && strcmp(px, "CFSPML") == 0)
	{
		CFS_flag = true;
		PMLX = 10;
		printf("X Boundary: CFSPML\n");
	}
	else if (strcmp(nx, "PBC") == 0 && strcmp(px, "PBC") == 0)
	{
		PBCX = true;
		PBC_flag = true;
		PMLX = 0;
		printf("X Boundary: PBC\n");
	}
	else
		; // 예외처리 필요.

	if (strcmp(ny, "PML") == 0 && strcmp(py, "PML") == 0)
	{
		PML_flag = true;
		PMLY = 10;
		printf("Y Boundary: PML\n");
	}
	else if (strcmp(ny, "CFSPML") == 0 && strcmp(py, "CFSPML") == 0)
	{
		CFS_flag = true;
		PMLY = 10;		
		printf("Y Boundary: CFSPML\n");
	}
	else if (strcmp(ny, "PBC") == 0 && strcmp(py, "PBC") == 0)
	{
		PBCY = true;
		PBC_flag = true;
		PMLY = 0;
		printf("Y Boundary: PBC\n");
	}
	else
		; // 예외처리 필요.

	if (strcmp(nz, "PML") == 0 && strcmp(pz, "PML") == 0)
	{
		PML_flag = true;
		PMLZ = 10;
		printf("Z Boundary: PML\n");
	}
	else if (strcmp(nz, "CFSPML") == 0 && strcmp(pz, "CFSPML") == 0)
	{
		CFS_flag = true;
		PMLZ = 10;
		printf("Z Boundary: CFSPML\n");
	}
	else if (strcmp(nz, "PBC") == 0 && strcmp(pz, "PBC") == 0)
	{
		PBCZ = true;
		PBC_flag = true;
		PMLZ = 0;
		printf("Z Boundary: PBC\n");
	}
	else
		; // 예외처리 필요.

	pFc->PML[0] = PMLX;
	pFc->PML[1] = PMLY;
	pFc->PML[2] = PMLZ;
}

void FDTD_Boundary_Condition::Read_Json()
{
	Document Jdocument;
	ifstream ifs(JsonFileName);
	IStreamWrapper isw(ifs);
	Jdocument.ParseStream(isw);
	strcpy(nx, Jdocument["BoundaryCondition"]["nx"].GetString());
	strcpy(px, Jdocument["BoundaryCondition"]["px"].GetString());
	strcpy(ny, Jdocument["BoundaryCondition"]["ny"].GetString());
	strcpy(py, Jdocument["BoundaryCondition"]["py"].GetString());
	strcpy(nz, Jdocument["BoundaryCondition"]["nz"].GetString());
	strcpy(pz, Jdocument["BoundaryCondition"]["pz"].GetString());
	ifs.close();
}

void FDTD_Boundary_Condition::Boundary_Condition_PMLF()
{
	if (PML_flag == true && PBC_flag==false)
	{
		pPML->Fx_update();
		pPML->Fy_update();
		pPML->Fz_update();
	}
	else if (CFS_flag == true && PBC_flag == false)
	{
		pCFSPML->Fx_update();
		pCFSPML->Fy_update();
		pCFSPML->Fz_update();
	}
	else if (PML_flag == true && PBC_flag == true)
	{
		pPBCScrewed->PBC_Z_F_PML_minus();
		pPBCScrewed->PBC_Z_F_PML_plus();
	}
}
void FDTD_Boundary_Condition::Boundary_Condition_PMLG() 
{
	if (PML_flag == true && PBC_flag == false)
	{
		pPML->Gx_update();
		pPML->Gy_update();
		pPML->Gz_update();
	}
	else if (CFS_flag == true && PBC_flag == false)
	{
		pCFSPML->Gx_update();
		pCFSPML->Gy_update();
		pCFSPML->Gz_update();
	}
	else if (PML_flag == true && PBC_flag == true)
	{
		pPBCScrewed->PBC_Z_G_PML_plus();
		pPBCScrewed->PBC_Z_G_PML_minus();
	}
}

void FDTD_Boundary_Condition::Boundary_Condition_PBC() 
{
	pPBCScrewed->FDTD_PBC_Chiral_BC1_omp();
	pPBCScrewed->FDTD_PBC_Chiral_BC2_omp();
}
