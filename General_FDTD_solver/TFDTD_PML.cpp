//---------------------------------------------------------------------------
#pragma hdrstop
#include "TFDTD_PML.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)

__fastcall TFDTD_PML::TFDTD_PML()
{
    isFDTD=0;
    this->eps_0 = 1e-9 / 36 / m_pi;
    this->mu_0 = 4 * m_pi * 1e-7;

	
}

TFDTD_PML::~TFDTD_PML()
{
    memfree(&fxy1, txsize, PMLY+1);
	memfree(&fxz1, txsize, tysize);
	memfree(&fxy2, txsize, PMLY+1);
	memfree(&fxz2, txsize, tysize);

	memfree(&fyx1, PMLX+1, tysize);
	memfree(&fyz1, txsize, tysize);
	memfree(&fyx2, PMLX+1, tysize);
	memfree(&fyz2, txsize, tysize);

	memfree(&fzx1, PMLX+1, tysize);
	memfree(&fzy1, txsize, PMLY+1);
	memfree(&fzx2, PMLX+1, tysize);
	memfree(&fzy2, txsize, PMLY+1);

	memfree(&gxy1, txsize+1, PMLY+1);
	memfree(&gxz1, txsize+1, tysize);
	memfree(&gxy2, txsize+1, PMLY+1);
	memfree(&gxz2, txsize+1, tysize);

	memfree(&gyz1, txsize, tysize+1);
	memfree(&gyx1, PMLX, tysize+1);
	memfree(&gyz2, txsize, tysize+1);
	memfree(&gyx2, PMLX, tysize+1);

	memfree(&gzx1, PMLX, tysize);
	memfree(&gzy1, txsize, PMLY);
	memfree(&gzx2, PMLX, tysize);
	memfree(&gzy2, txsize, PMLY);

	memfree(&Cxa1);  memfree(&Cxa2);
	memfree(&Cxb1);  memfree(&Cxb2);
	memfree(&Cya1);  memfree(&Cya2);
	memfree(&Cyb1);  memfree(&Cyb2);
	memfree(&Cza1);  memfree(&Cza2);
	memfree(&Czb1);  memfree(&Czb2);

	memfree(&Dxa1);   memfree(&Dxa2);
	memfree(&Dxb1);   memfree(&Dxb2);
	memfree(&Dya1);   memfree(&Dya2);
	memfree(&Dyb1);   memfree(&Dyb2);
	memfree(&Dza1);   memfree(&Dza2);
	memfree(&Dzb1);   memfree(&Dzb2);
}


void TFDTD_PML::Set_Kappa_coefficient(double alphaT, double kappaT)
{
	this->alpha_max = (double)alphaT;
	this->kappa_max = (double)kappaT;
}

//void TFDTD_PML::Set_Structure_and_Electric_property()
//{
//	Read_Material_Data();
//	Read_Media_Data();
//	dt = 0.99 / c_0 / sqrt((1 / min1D(dx)) * (1 / min1D(dx)) + (1 / min1D(dy)) * (1 / min1D(dy)) + (1 / min1D(dz)) * (1 / min1D(dz)));
//
//	// Cd, Ch ����
//	double temp = 0;
//	C_d = Tptr1<double>(num_material);
//	for (int ii = 0; ii < num_material ; ii++)
//	{
//		// ��л�
//		if ( pMp[ii].TypeOfMaterials == 'N')
//		{
//			C_d[ii] = dt / ((2 * eps_0 * pMp[ii].Epsr) + (pMp[ii].Sigma * dt));
//		}
//		//�л�
//		if (pMp[ii].TypeOfMaterials == 'D')
//		{
//			for (int pp = 0; pp < total_pol; pp++)
//				temp = pMp[ii].Cc[pp];
//			C_d[ii] = 1 / (pMp[ii].Epsr*eps_0 / dt + temp * 0.5);	
//		}	
//	}
//	//if (child == 0)
//	//{
//	//	this->C_d = dt /( (2*eps_0*eps_r)+(sigma_e*dt) );
//	//	this->C_h = dt / 2.0 / mu_0;
//	//}
//	//else if (child == 1)
//	//{
//	//	this->C_d = dt / 2.0; // dt /( (2*eps_0*eps_r)+(sigma_e*dt) );
//	//	this->C_h = dt / mu_0;
//	//}
//	
//	// PML_Initial_Setting();
//
//}

void TFDTD_PML::PML_field_initiate()
{
	// PML - Coef
	Cxa1 = Tptr1<double>(PMLX + 1);  Cxa2 = Tptr1<double>(PMLX + 1);
	Cxb1 = Tptr1<double>(PMLX + 1);  Cxb2 = Tptr1<double>(PMLX + 1);
	Cya1 = Tptr1<double>(PMLY + 1);  Cya2 = Tptr1<double>(PMLY + 1);
	Cyb1 = Tptr1<double>(PMLY + 1);  Cyb2 = Tptr1<double>(PMLY + 1);
	Cza1 = Tptr1<double>(PMLZ + 1);  Cza2 = Tptr1<double>(PMLZ + 1);
	Czb1 = Tptr1<double>(PMLZ + 1);  Czb2 = Tptr1<double>(PMLZ + 1);

	Dxa1 = Tptr1<double>(PMLX + 1);   Dxa2 = Tptr1<double>(PMLX + 1);
	Dxb1 = Tptr1<double>(PMLX + 1);   Dxb2 = Tptr1<double>(PMLX + 1);
	Dya1 = Tptr1<double>(PMLY + 1);   Dya2 = Tptr1<double>(PMLY + 1);
	Dyb1 = Tptr1<double>(PMLY + 1);   Dyb2 = Tptr1<double>(PMLY + 1);
	Dza1 = Tptr1<double>(PMLZ + 1);   Dza2 = Tptr1<double>(PMLZ + 1);
	Dzb1 = Tptr1<double>(PMLZ + 1);   Dzb2 = Tptr1<double>(PMLZ + 1);

	// PML - Fields
	fxy1 = Tptr3<double>(txsize, PMLY + 1, tzsize);
	fxz1 = Tptr3<double>(txsize, tysize, PMLZ + 1);
	fxy2 = Tptr3<double>(txsize, PMLY + 1, tzsize);
	fxz2 = Tptr3<double>(txsize, tysize, PMLZ + 1);

	fyx1 = Tptr3<double>(PMLX + 1, tysize, tzsize);
	fyz1 = Tptr3<double>(txsize, tysize, PMLZ + 1);
	fyx2 = Tptr3<double>(PMLX + 1, tysize, tzsize);
	fyz2 = Tptr3<double>(txsize, tysize, PMLZ + 1);

	fzx1 = Tptr3<double>(PMLX + 1, tysize, tzsize);
	fzy1 = Tptr3<double>(txsize, PMLY + 1, tzsize);
	fzx2 = Tptr3<double>(PMLX + 1, tysize, tzsize);
	fzy2 = Tptr3<double>(txsize, PMLY + 1, tzsize);

	gxy1 = Tptr3<double>(txsize + 1, PMLY + 1, tzsize);
	gxz1 = Tptr3<double>(txsize + 1, tysize, PMLZ + 1);
	gxy2 = Tptr3<double>(txsize + 1, PMLY + 1, tzsize);
	gxz2 = Tptr3<double>(txsize + 1, tysize, PMLZ + 1);

	gyz1 = Tptr3<double>(txsize, tysize + 1, PMLZ + 1);
	gyx1 = Tptr3<double>(PMLX + 1, tysize + 1, tzsize);
	gyz2 = Tptr3<double>(txsize, tysize + 1, PMLZ + 1);
	gyx2 = Tptr3<double>(PMLX + 1, tysize + 1, tzsize);

	gzx1 = Tptr3<double>(PMLX + 1, tysize, tzsize + 1);
	gzy1 = Tptr3<double>(txsize, PMLY + 1, tzsize + 1);
	gzx2 = Tptr3<double>(PMLX + 1, tysize, tzsize + 1);
	gzy2 = Tptr3<double>(txsize, PMLY + 1, tzsize + 1);
}

void TFDTD_PML::PML_Initial_Setting()
{
	this->Ex = pFc->Ex_address;
	this->Ey = pFc->Ey_address;
	this->Ez = pFc->Ez_address;

	this->Hx = pFc->Hx_address;
	this->Hy = pFc->Hy_address;
	this->Hz = pFc->Hz_address;

	this->dx = pFc->dx;
	this->dy = pFc->dy;
	this->dz = pFc->dz;
	this->E_dx = pFc->E_dx;
	this->E_dy = pFc->E_dy;
	this->E_dz = pFc->E_dz;
	this->H_dx = pFc->H_dx;
	this->H_dy = pFc->H_dy;
	this->H_dz = pFc->H_dz;

	this->Media_Ex = pFc->Media_Ex;
	this->Media_Ey = pFc->Media_Ey;
	this->Media_Ez = pFc->Media_Ez;

	this->dt = pFc->dt;

	this->txsize = pFc->txsize;
	this->tysize = pFc->tysize;
	this->tzsize = pFc->tzsize;
	this->PMLX = pFc->PML[0];
	this->PMLY = pFc->PML[1];
	this->PMLZ = pFc->PML[2];

	Set_Structure_and_Electric_property();
	PML_field_initiate();
	PML_parameter_initiate();

}

void TFDTD_PML::PML_parameter_initiate()
{
	if (PMLX != 0)
	{
		PML_parameter_initiateX();
	}
	if (PMLY != 0)
	{
		PML_parameter_initiateY();
	}
	if (PMLZ != 0)
	{
		PML_parameter_initiateZ();
	}
}

void TFDTD_PML::Fx_update()
{
	X_NFpml();
	X_PFpml();
} 
void TFDTD_PML::Gx_update()
{
	X_NGpml();
	X_PGpml();
}
void TFDTD_PML::Fy_update()
{
	Y_NFpml();
	Y_PFpml();
}

void TFDTD_PML::Gy_update()
{
	Y_NGpml();
	Y_PGpml();
}
void TFDTD_PML::Fz_update()
{
	Z_NFpml();
	Z_PFpml();
}
void TFDTD_PML::Gz_update()
{
	Z_NGpml();
	Z_PGpml();
}
