//---------------------------------------------------------------------------
#pragma hdrstop
#include "TFDTD_CFS_PML.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)

TFDTD_CFS_PML::TFDTD_CFS_PML(struct_Field_carrier* pFc)
{
	this->pFc = pFc;
	this->dt = pFc->dt;
	Read_Material_Data();
	PML_Initial_Setting();
}

TFDTD_CFS_PML::~TFDTD_CFS_PML(){}

void TFDTD_CFS_PML::Set_Structure_and_Electric_property()
{
// Cd, Ch 정의
	double temp = 0;
	C_d = Tptr1<double>(num_material);
	for (int ii = 0; ii < num_material; ii++)
	{
		// 비분산
		if (pMp[ii].TypeOfMaterials == 'N')
		{
			C_d[ii] = 2*dt / ((2 * eps_0 * pMp[ii].Epsr) + (pMp[ii].Sigma * dt));
			C_h = dt / mu_0 ;
		}
		//분산
		if (pMp[ii].TypeOfMaterials == 'D')
		{
			for (int pp = 0; pp < total_pol; pp++)
				temp = pMp[ii].Cc[pp];
			C_d[ii] = 1 / (pMp[ii].Epsr * eps_0 / dt + temp * 0.5);
			C_h = dt / mu_0;
		}
	}
	//if (child == 0)
	//{
	//	this->C_d = dt /( (2*eps_0*eps_r)+(sigma_e*dt) );
	//	this->C_h = dt / 2.0 / mu_0;
	//}
	//else if (child == 1)
	//{
	//	this->C_d = dt / 2.0; // dt /( (2*eps_0*eps_r)+(sigma_e*dt) );
	//	this->C_h = dt / mu_0;
	//}

	// PML_Initial_Setting();

}

void TFDTD_CFS_PML::PML_parameter_initiateX()
{
    m_sigma=4; // Degree of polynomial
	m_kappa=4;
    m = 4.0;
    eta_0 = sqrt(mu_0 / eps_0);
	eps_r4pml = 1.0;
	mu_r = 1;
	
	kappa_e_xprof = Tptr1<double>(txsize); // profile in x-direction
	kappa_h_xprof = Tptr1<double>(txsize);
	
	Kappa_parameter_initiateX();
	
	double sigma_max1 = 0.8 * (m + 1.0) / (eta_0 * dx[0] * sqrt(mu_r * eps_r4pml)); // maximum conductivity
	double sigma_max2 = 0.8 * (m + 1.0) / (eta_0 * dx[txsize-1] * sqrt(mu_r * eps_r4pml)); // maximum conductivity
	// profile in x-direction
	double* sigma_e_xprof1 = Tptr1<double>(PMLX + 1); 
	double* sigma_e_xprof2 = Tptr1<double>(PMLX + 1); 
	double* sigma_h_xprof1 = Tptr1<double>(PMLX + 1);
	double* sigma_h_xprof2 = Tptr1<double>(PMLX + 1);
	double* alpha_e_xprof = Tptr1<double>(PMLX + 1);
	double* alpha_h_xprof = Tptr1<double>(PMLX + 1);

	for(i=0;i<PMLX+1;i++) 
	{
		// setup conductivity profile in x-direction
		sigma_e_xprof1[i] = sigma_max1 * pow((i + 0.0) / PMLX, m_sigma);
		sigma_e_xprof2[i] = sigma_max2 * pow((i + 0.0) / PMLX, m_sigma);
		sigma_h_xprof1[i] = sigma_max1 * pow((i + 1.0 / 2.0) / PMLX, m_sigma);
		sigma_h_xprof2[i] = sigma_max2 * pow((i + 1.0 / 2.0) / PMLX, m_sigma);

		alpha_e_xprof[i] = alpha_max * (1. - (i + 0.0) / PMLX); 
		alpha_h_xprof[i] = alpha_max * (1. - (i + 1.0 / 2.0) / PMLX);
	}

	// X-direction
	for(i=0;i<PMLX+1;i++) 
    {
        ai=PMLX-i; //recursive index
		temp= exp(-(sigma_e_xprof1[ai]/kappa_e_xprof[ai] + alpha_e_xprof[ai])*dt/(eps_0));
		Cxa1[i] = temp;
		Cxb1[i] = (temp-1) * sigma_e_xprof1[ai]/(dx[0]* kappa_e_xprof[ai]*(sigma_e_xprof1[ai] + kappa_e_xprof[ai]*alpha_e_xprof[ai]));
	}

	for(i=0;i<PMLX+1;i++) 
	{
		temp= exp(-(sigma_e_xprof2[i]/kappa_e_xprof[i] + alpha_e_xprof[i])*dt/(eps_0));
		Cxa2[i] = temp;
		Cxb2[i] = (temp-1) * sigma_e_xprof2[i]/(dx[txsize-1]* kappa_e_xprof[i]*(sigma_e_xprof2[i] + kappa_e_xprof[i]*alpha_e_xprof[i]));
	}

	// X-direction
	for (i = 0; i < PMLX; i++) 
	{
		bi = PMLX - 1 - i; //recursive index
		temp = exp(-(sigma_h_xprof1[bi] / kappa_h_xprof[bi] + alpha_h_xprof[bi]) * dt / (eps_0));
		Dxa1[i] = temp;
		Dxb1[i] = (temp - 1)
			* sigma_h_xprof1[bi] / (dx[0] * kappa_h_xprof[bi] * (sigma_h_xprof1[bi] + kappa_h_xprof[bi] * alpha_h_xprof[bi]));

		temp = exp(-(sigma_h_xprof2[i] / kappa_h_xprof[i] + alpha_h_xprof[i]) * dt / (eps_0));
		Dxa2[i] = temp;
		Dxb2[i] = (temp - 1)
			* sigma_h_xprof2[i] / (dx[txsize-1] * kappa_h_xprof[i] * (sigma_h_xprof2[i] + kappa_h_xprof[i] * alpha_h_xprof[i]));
	}

     memfree(&sigma_e_xprof1);
     memfree(&sigma_h_xprof1);
	 memfree(&sigma_e_xprof2);
	 memfree(&sigma_h_xprof2);
     memfree(&alpha_e_xprof);
     memfree(&alpha_h_xprof);

}

void TFDTD_CFS_PML::PML_parameter_initiateY()
{
    m_sigma=4; // Degree of polynomial
	m_kappa=4;
    m = 4.0;
    eta_0 = sqrt(mu_0 / eps_0);
	eps_r4pml = 1.0;
	mu_r = 1;

	kappa_e_yprof = Tptr1<double>(tysize); // profile in x-direction
	kappa_h_yprof = Tptr1<double>(tysize);

	Kappa_parameter_initiateY();

	double sigma_max1 = 0.8 * (m + 1.0) / (eta_0 * dy[0] * sqrt(mu_r * eps_r4pml)); // maximum conductivity
	double sigma_max2 = 0.8 * (m + 1.0) / (eta_0 * dy[tysize - 1] * sqrt(mu_r * eps_r4pml)); // maximum conductivity
	// profile in x-direction
	double* sigma_e_yprof1 = Tptr1<double>(PMLY + 1);
	double* sigma_e_yprof2 = Tptr1<double>(PMLY + 1);
	double* sigma_h_yprof1 = Tptr1<double>(PMLY + 1);
	double* sigma_h_yprof2 = Tptr1<double>(PMLY + 1);
	double* alpha_e_yprof  = Tptr1<double>(PMLY + 1);
	double* alpha_h_yprof  = Tptr1<double>(PMLY + 1);

	for (i = 0; i < PMLY + 1; i++)
	{
		// setup conductivity profile in x-direction
		sigma_e_yprof1[i] = sigma_max1 * pow((i + 0.0) / PMLY, m_sigma);
		sigma_e_yprof2[i] = sigma_max2 * pow((i + 0.0) / PMLY, m_sigma);
		sigma_h_yprof1[i] = sigma_max1 * pow((i + 1.0 / 2.0) / PMLY, m_sigma);
		sigma_h_yprof2[i] = sigma_max2 * pow((i + 1.0 / 2.0) / PMLY, m_sigma);
		alpha_e_yprof[i] = alpha_max * (1. - (i + 0.0) / PMLY);
		alpha_h_yprof[i] = alpha_max * (1. - (i + 1.0 / 2.0) / PMLY);
	}

   // Y-direction
	for(j=1;j<PMLY+1;j++) {
		aj=PMLY-j; //recursive index
			temp= exp(-(sigma_e_yprof1[aj]/ kappa_e_yprof[aj] + alpha_e_yprof[aj])*dt/(eps_0));
			Cya1[j] = temp;
			Cyb1[j] = (temp-1)
				       * sigma_e_yprof1[aj]/(dy[0]* kappa_e_yprof[aj]*(sigma_e_yprof1[aj]+kappa_e_yprof[aj]*alpha_e_yprof[aj]));

	}

	for(j=0;j<PMLY+1;j++) {
			temp= exp(-(sigma_e_yprof2[j]/ kappa_e_yprof[j] + alpha_e_yprof[j])*dt/(eps_0));
			Cya2[j] = temp;
			Cyb2[j] = (temp-1)
				         * sigma_e_yprof2[j]/(dy[tysize-1]* kappa_e_yprof[j]*(sigma_e_yprof2[j]+ kappa_e_yprof[j]*alpha_e_yprof[j]));
	}

	for (j = 0; j < PMLY; j++) 
	{
		bj = PMLY - 1 - j; //recursive index
		temp = exp(-(sigma_h_yprof1[bj] / kappa_h_yprof[bj] + alpha_h_yprof[bj]) * dt / (eps_0));
		Dya1[j] = temp;
		Dyb1[j] = (temp - 1)
			* sigma_h_yprof1[bj] / (dy[0] * kappa_h_yprof[bj] * (sigma_h_yprof1[bj] + kappa_h_yprof[bj] * alpha_h_yprof[bj]));

		temp = exp(-(sigma_h_yprof2[j] / kappa_h_yprof[j] + alpha_h_yprof[j]) * dt / (eps_0));
		Dya2[j] = temp;
		Dyb2[j] = (temp - 1)
			* sigma_h_yprof2[j] / (dy[tysize-1] * kappa_h_yprof[j] * (sigma_h_yprof2[j] + kappa_h_yprof[j] * alpha_h_yprof[j]));
	}
     memfree(&sigma_e_yprof1);
     memfree(&sigma_h_yprof1);
	 memfree(&sigma_e_yprof2);
	 memfree(&sigma_h_yprof2);
     memfree(&kappa_e_yprof);
     memfree(&kappa_h_yprof);
     memfree(&alpha_e_yprof);
     memfree(&alpha_h_yprof);
}

void TFDTD_CFS_PML::PML_parameter_initiateZ()
{
	m_sigma = 4; // Degree of polynomial
	m_kappa = 4;
	m = 4.0;
	eta_0 = sqrt(mu_0 / eps_0);
	eps_r4pml = 1.0;
	mu_r = 1;

	kappa_e_zprof = Tptr1<double>(tzsize); // profile in x-direction
	kappa_h_zprof = Tptr1<double>(tzsize);

	Kappa_parameter_initiateY();

	double sigma_max1 = 0.8 * (m + 1.0) / (eta_0 * dz[0]          * sqrt(mu_r * eps_r4pml)); // maximum conductivity
	double sigma_max2 = 0.8 * (m + 1.0) / (eta_0 * dz[tzsize - 1] * sqrt(mu_r * eps_r4pml)); // maximum conductivity
	// profile in x-direction
	double* sigma_e_zprof1 = Tptr1<double>(PMLZ + 1);
	double* sigma_e_zprof2 = Tptr1<double>(PMLZ + 1);
	double* sigma_h_zprof1 = Tptr1<double>(PMLZ + 1);
	double* sigma_h_zprof2 = Tptr1<double>(PMLZ + 1);
	double* alpha_e_zprof  = Tptr1<double>(PMLZ + 1);
	double* alpha_h_zprof  = Tptr1<double>(PMLZ + 1);

	for (i = 0; i < PMLZ + 1; i++)
	{
		// setup conductivity profile in x-direction
		sigma_e_zprof1[i] = sigma_max1 * pow((i + 0.0) / PMLZ, m_sigma);
		sigma_e_zprof2[i] = sigma_max2 * pow((i + 0.0) / PMLZ, m_sigma);
		sigma_h_zprof1[i] = sigma_max1 * pow((i + 1.0 / 2.0) / PMLZ, m_sigma);
		sigma_h_zprof2[i] = sigma_max2 * pow((i + 1.0 / 2.0) / PMLZ, m_sigma);
		alpha_e_zprof[i] = alpha_max * (1. - (i + 0.0) / PMLZ);
		alpha_h_zprof[i] = alpha_max * (1. - (i + 1.0 / 2.0) / PMLZ);
	}

    // Z-direction
	for(k=0;k<PMLZ+1;k++) 
	{
		ak=PMLZ-k; //recursive index
		temp= exp(-(sigma_e_zprof1[ak]/kappa_e_zprof[ak] + alpha_e_zprof[ak])*dt/(eps_0));
		Cza1[k] = temp;
		Czb1[k] = (temp-1)
				        * sigma_e_zprof1[ak]/(dz[0]* kappa_e_zprof[ak]*(sigma_e_zprof1[ak]+ kappa_e_zprof[ak]*alpha_e_zprof[ak]));

	}
	for(k=0;k<PMLZ+1;k++) {
			temp= exp(-(sigma_e_zprof2[k]/ kappa_e_zprof[k] + alpha_e_zprof[k])*dt/(eps_0));
			Cza2[k] = temp;
			Czb2[k] = (temp-1)
				         * sigma_e_zprof2[k]/(dz[tzsize-1]*kappa_e_zprof[k]*(sigma_e_zprof2[k]+ kappa_e_zprof[k]*alpha_e_zprof[k]));
	}

    for(k=0;k<PMLZ;k++) {

	    bk=PMLZ-1-k; //recursive index
		temp= 	exp(-(sigma_h_zprof1[bk]/kappa_h_zprof[bk] + alpha_h_zprof[bk])*dt/(eps_0));
		Dza1[k] = temp;
		Dzb1[k] = (temp-1)
			         * sigma_h_zprof1[bk]/(dz[0]* kappa_h_zprof[bk]*(sigma_h_zprof1[bk]+ kappa_h_zprof[bk]*alpha_h_zprof[bk]));

		//printf("%e\t%e\t%\n",Dza1[k],Dzb1[k]);

		temp=	exp(-(sigma_h_zprof2[k]/ kappa_h_zprof[k] + alpha_h_zprof[k])*dt/(eps_0));
		Dza2[k] = temp;
		Dzb2[k] = (temp-1)
			         * sigma_h_zprof2[k]/(dz[tzsize-1]* kappa_h_zprof[k]*(sigma_h_zprof2[k]+ kappa_h_zprof[k]*alpha_h_zprof[k]));
	}

	memfree(&sigma_e_zprof1);
	memfree(&sigma_h_zprof1);
	memfree(&sigma_e_zprof2);
	memfree(&sigma_h_zprof2);
	memfree(&kappa_e_zprof);
	memfree(&kappa_h_zprof);
	memfree(&alpha_e_zprof);
	memfree(&alpha_h_zprof);
}

void TFDTD_CFS_PML::X_NFpml()
{
	int i = 0, j = 0, k = 0;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i = 1; i < PMLX + 1; i++)
		for (j = 0; j < tysize; j++)
		for (k = 0; k < tzsize; k++)
		{
			my = Media_Ey[i][j][k];
			mz = Media_Ez[i][j][k];;
			fyx1[i][j][k] = Cxa1[i] * fyx1[i][j][k] + Cxb1[i] * (Hz[i][j][k] - Hz[i - 1][j][k]);
			fzx1[i][j][k] = Cxa1[i] * fzx1[i][j][k] + Cxb1[i] * (Hy[i][j][k] - Hy[i - 1][j][k]);
			Ey[i][j][k] = Ey[i][j][k] - C_d[my] * fyx1[i][j][k];
			Ez[i][j][k] = Ez[i][j][k] + C_d[mz] * fzx1[i][j][k];
		}
	}
}

void TFDTD_CFS_PML::X_NGpml()
{
	int i = 0, j = 0, k = 0;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i = 0; i < PMLX; i++)
		for (j = 0; j < tysize; j++)
		for (k = 0; k < tzsize; k++)
		{
			gyx1[i][j][k] = Dxa1[i] * gyx1[i][j][k] + Dxb1[i] * (Ez[i + 1][j][k] - Ez[i][j][k]);
			gzx1[i][j][k] = Dxa1[i] * gzx1[i][j][k] + Dxb1[i] * (Ey[i + 1][j][k] - Ey[i][j][k]);
			Hy[i][j][k] = Hy[i][j][k] + C_h * gyx1[i][j][k];
			Hz[i][j][k] = Hz[i][j][k] - C_h * gzx1[i][j][k];
		}
	}
}

void TFDTD_CFS_PML::X_PFpml()
{
	int ai = 0, i = 0, j = 0, k = 0;
#pragma omp parallel private(i,j,k, ai)
	{
#pragma omp for schedule(dynamic)
		for (ai = 0; ai < PMLX; ai++)
		{
			i = txsize - PMLX + ai;
			for (j = 0; j < tysize; j++)
			for (k = 0; k < tzsize; k++)
			{
				my = Media_Ey[i][j][k];
				mz = Media_Ez[i][j][k];
				fyx2[ai][j][k] = Cxa2[ai] * fyx2[ai][j][k] + Cxb2[ai] * (Hz[i][j][k] - Hz[i - 1][j][k]);
				fzx2[ai][j][k] = Cxa2[ai] * fzx2[ai][j][k] + Cxb2[ai] * (Hy[i][j][k] - Hy[i - 1][j][k]);
				Ey[i][j][k] = Ey[i][j][k] - C_d[my] * fyx2[ai][j][k];
				Ez[i][j][k] = Ez[i][j][k] + C_d[mz] * fzx2[ai][j][k];
			}
		}
	}
}

void TFDTD_CFS_PML::X_PGpml()
{
	int bi = 0, i = 0, j = 0, k = 0;
#pragma omp parallel private(i,j,k, bi)
	{
#pragma omp for schedule(dynamic)
		for (bi = 0; bi < PMLX; bi++) {
			i = txsize - PMLX + bi;
			for (j = 0; j < tysize; j++)
				for (k = 0; k < tzsize; k++)
				{
					gyx2[bi][j][k] = Dxa2[bi] * gyx2[bi][j][k] + Dxb2[bi] * (Ez[i + 1][j][k] - Ez[i + 1][j][k]);
					gzx2[bi][j][k] = Dxa2[bi] * gzx2[bi][j][k] + Dxb2[bi] * (Ey[i + 1][j][k] - Ey[i + 1][j][k]);
					Hy[i][j][k] = Hy[i][j][k] + C_h * gyx2[bi][j][k];
					Hz[i][j][k] = Hz[i][j][k] - C_h * gzx2[bi][j][k];
				}

		}
	}
}

void TFDTD_CFS_PML::Y_NFpml()
{
	int i = 0, j = 0, k = 0;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i = 0; i < txsize; i++)
			for (j = 0; j < PMLY; j++)
				for (k = 0; k < tzsize; k++)
				{
					mx = Media_Ex[i][j][k];
					mz = Media_Ez[i][j][k];
					fxy1[i][j][k] = Cya1[j] * fxy1[i][j][k] + Cyb1[j] * (Hz[i][j][k] - Hz[i][j - 1][k]);
					fzy1[i][j][k] = Cya1[j] * fzy1[i][j][k] + Cyb1[j] * (Hx[i][j][k] - Hx[i][j - 1][k]);
					Ex[i][j][k] = Ex[i][j][k] + C_d[mx] * fxy1[i][j][k];
					Ez[i][j][k] = Ez[i][j][k] - C_d[mz] * fzy1[i][j][k];
				}
	}
}

void TFDTD_CFS_PML::Y_NGpml()
{
	int i = 0, j = 0, k = 0;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i = 1; i < txsize; i++)
			for (j = 0; j < PMLY; j++)
				for (k = 0; k < tzsize; k++)
				{
					gxy1[i][j][k] = Dya1[j] * gxy1[i][j][k] + Dyb1[j] * (Ez[i][j + 1][k] - Ez[i][j][k]);
					gzy1[i][j][k] = Dya1[j] * gzy1[i][j][k] + Dyb1[j] * (Ex[i][j + 1][k] - Ex[i][j][k]);
					Hz[i][j][k] = Hz[i][j][k] + C_h * gzy1[i][j][k];
					Hx[i][j][k] = Hx[i][j][k] - C_h * gxy1[i][j][k];
				}
	}
}

void TFDTD_CFS_PML::Y_PFpml()
{
	int aj = 0, i = 0, j = 0, k = 0;
#pragma omp parallel private(i,j,k, aj)
	{
#pragma omp for schedule(dynamic)
		for (i = 0; i < txsize; i++)
			for (aj = 0; aj < PMLY; aj++)
			{
				j = tysize - PMLY + aj;
				for (k = 0; k < tzsize; k++)
				{
					mx = Media_Ex[i][j][k];
					mz = Media_Ez[i][j][k];
					fxy2[i][aj][k] = Cya2[aj] * fxy2[i][aj][k] + Cyb2[aj] * (Hz[i][j][k] - Hz[i][j - 1][k]);
					fzy2[i][aj][k] = Cya2[aj] * fzy2[i][aj][k] + Cyb2[aj] * (Hx[i][j][k] - Hx[i][j - 1][k]);
					Ex[i][j][k] = Ex[i][j][k] + C_d[mx] * fxy2[i][aj][k];
					Ez[i][j][k] = Ez[i][j][k] - C_d[mz] * fzy2[i][aj][k];
				}
			}
	}
}

void TFDTD_CFS_PML::Y_PGpml()
{
	int bj = 0, i = 0, j = 0, k = 0;
#pragma omp parallel private(i,j,k, bj)
	{
#pragma omp for schedule(dynamic)
		for (i = 0; i < txsize; i++)
			for (bj = 0; bj < PMLY; bj++)
			{
				j = tysize - PMLY + bj;
				for (k = 0; k < tzsize; k++)
				{
					gxy2[i][bj][k] = Dya2[bj] * gxy2[i][bj][k] + Dyb2[bj] * (Ez[i][j + 1][k] - Ez[i][j][k]);
					gzy2[i][bj][k] = Dya2[bj] * gzy2[i][bj][k] + Dyb2[bj] * (Ex[i][j + 1][k] - Ex[i][j][k]);
					Hx[i][j][k] = Hx[i][j][k] - C_h * gxy2[i][bj][k];
					Hz[i][j][k] = Hz[i][j][k] + C_h * gzy2[i][bj][k];
				}
			}
	}

}

void TFDTD_CFS_PML::Z_NFpml()
{
	int i = 0, j = 0, k = 0;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i = 0; i < txsize; i++)
			for (j = 0; j < tysize; j++)
				for (k = 1; k < PMLZ + 1; k++)
				{
					mx = Media_Ex[i][j][k];
					my = Media_Ey[i][j][k];
					fxz1[i][j][k] = Cza1[k] * fxz1[i][j][k] + Czb1[k] * (Hy[i][j][k] - Hy[i][j][k - 1]);
					fyz1[i][j][k] = Cza1[k] * fyz1[i][j][k] + Czb1[k] * (Hx[i][j][k] - Hx[i][j][k - 1]);
					Ex[i][j][k] = Ex[i][j][k] - C_d[mx] * fxz1[i][j][k];
					Ey[i][j][k] = Ey[i][j][k] + C_d[my] * fyz1[i][j][k];
				}
	}
}

void TFDTD_CFS_PML::Z_NGpml()
{
	int i = 0, j = 0, k = 0;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i = 0; i < txsize; i++)
			for (j = 0; j < tysize; j++)
				for (k = 0; k < PMLZ; k++)
				{
					gxz1[i][j][k] = Dza1[k] * gxz1[i][j][k] + Dzb1[k] * (Ey[i][j][k + 1] - Ey[i][j][k]);
					gyz1[i][j][k] = Dza1[k] * gyz1[i][j][k] + Dzb1[k] * (Ex[i][j][k + 1] - Ex[i][j][k]);
					Hx[i][j][k] = Hx[i][j][k] + C_h * gxz1[i][j][k];
					Hy[i][j][k] = Hy[i][j][k] - C_h * gyz1[i][j][k];
				}
	}
}

void TFDTD_CFS_PML::Z_PFpml()
{
	int ak = 0, i = 0, j = 0, k = 0;
#pragma omp parallel private(i,j,k, ak)
	{
#pragma omp for schedule(dynamic)
		for (i = 0; i < txsize; i++)
			for (j = 0; j < tysize; j++)
				for (ak = 0; ak < PMLZ; ak++)
				{
					k = tzsize - PMLZ + ak;
					mx = Media_Ex[i][j][k];
					my = Media_Ey[i][j][k];
					fxz2[i][j][ak] = Cza2[ak] * fxz2[i][j][ak] + Czb2[ak] * (Hy[i][j][k] - Hy[i][j][k - 1]);
					fyz2[i][j][ak] = Cza2[ak] * fyz2[i][j][ak] + Czb2[ak] * (Hx[i][j][k] - Hx[i][j][k - 1]);
					Ex[i][j][k] = Ex[i][j][k] - C_d[mx] * fxz2[i][j][ak];
					Ey[i][j][k] = Ey[i][j][k] + C_d[my] * fyz2[i][j][ak];
				}
	}
}

void TFDTD_CFS_PML::Z_PGpml()
{
	int bk = 0, i = 0, j = 0, k = 0;
#pragma omp parallel private(i,j,k, bk)
	{
#pragma omp for schedule(dynamic)
		for (i = 0; i < txsize; i++)
			for (j = 0; j < tysize; j++)
				for (bk = 0; bk < PMLZ; bk++)
				{
					k = tzsize - PMLZ + bk;
					gxz2[i][j][bk] = Dza2[bk] * gxz2[i][j][bk] + Dzb2[bk] * (Ey[i][j][k + 1] - Ey[i][j][k]);
					gyz2[i][j][bk] = Dza2[bk] * gyz2[i][j][bk] + Dzb2[bk] * (Ex[i][j][k + 1] - Ex[i][j][k]);
					Hx[i][j][k] = Hx[i][j][k] + C_h * gxz2[i][j][bk];
					Hy[i][j][k] = Hy[i][j][k] - C_h * gyz2[i][j][bk];
				}
	}
}
