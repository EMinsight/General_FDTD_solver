//---------------------------------------------------------------------------
#pragma hdrstop
#include "TFDTD_2PML.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)

TFDTD_2PML::TFDTD_2PML(struct_Field_carrier* pFc)
{
	this->pFc = pFc;
	this->dt = pFc->dt;
	Read_Material_Data();
	PML_Initial_Setting();
}

TFDTD_2PML::TFDTD_2PML(){}
TFDTD_2PML::~TFDTD_2PML(){}
void TFDTD_2PML::Set_Structure_and_Electric_property()
{
	// Cd, Ch ����
	double temp = 0;
	C_d = Tptr1<double>(num_material);
	for (int ii = 0; ii < num_material; ii++)
	{
		// ��л�
		if (pMp[ii].TypeOfMaterials == 'N')
		{
			C_d[ii] = dt / ((2 * eps_0 * pMp[ii].Epsr) + (pMp[ii].Sigma * dt));
			C_h = dt / mu_0 / 2;
		}
		//�л�
		if (pMp[ii].TypeOfMaterials == 'D')
		{
			for (int pp = 0; pp < total_pol; pp++)
				temp = pMp[ii].Cc[pp];
			C_d[ii] = 0.5 / (pMp[ii].Epsr * eps_0 / dt + temp * 0.5);
			C_h = dt / mu_0 / 2;
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

void TFDTD_2PML::PML_parameter_initiateX()
{
	double eps_r4pml = 1.0;
	double mu_r = 1;
	double m = 4.0;
	double eta_0 = sqrt(mu_0 / eps_0);

	double sigma_max1 = 0.8 * (m + 1.0) / (eta_0 * dx[0] * sqrt(mu_r * eps_r4pml)); // maximum conductivity
	double sigma_max2 = 0.8 * (m + 1.0) / (eta_0 * dx[txsize - 1] * sqrt(mu_r * eps_r4pml)); // maximum conductivity

	double* sigma_e_x1 = Tptr1<double>(PMLX + 1); // profile in x-direction
	double* sigma_e_x2 = Tptr1<double>(PMLX + 1); // profile in x-direction
	double* sigma_h_x1 = Tptr1<double>(PMLX + 1);
	double* sigma_h_x2 = Tptr1<double>(PMLX + 1);

	// Degree of polynomial
	for (i = 0; i<PMLX + 1; i++)
	{
		sigma_e_x1[i] = sigma_max1 * pow((i + 0.0) / PMLX, m);
		sigma_e_x2[i] = sigma_max2 * pow((i + 0.0) / PMLX, m);
		sigma_h_x1[i] = sigma_max1 * pow((i + 1.0 / 2.0) / PMLX, m);
		sigma_h_x2[i] = sigma_max2 * pow((i + 1.0 / 2.0) / PMLX, m);
	}

	// setup kappa_e, Ca,Cb for updating E fields in PML
	// X-direction
	for (i = 0; i<PMLX + 1; i++)
	{
		ai = PMLX - i; //recursive index
		Cxa1[i] = (1 - (sigma_e_x1[ai] * dt) / (2 * eps_0)) / (1 + (sigma_e_x1[ai] * dt) / (2 * eps_0));
		Cxb1[i] = ((sigma_e_x1[ai] * dt) / (eps_0*dx[0])) / (1 + (sigma_e_x1[ai] * dt) / (2 * eps_0)); //X-

		Cxa2[i] = (1 - (sigma_e_x2[i] * dt) / (2 * eps_0)) / (1 + (sigma_e_x2[i] * dt) / (2 * eps_0));
		Cxb2[i] = ((sigma_e_x2[i] * dt) / (eps_0*dx[txsize-1])) / (1 + (sigma_e_x2[i] * dt) / (2 * eps_0));   //X+
	}

	// setup kappa_h, Da,Db for updating H fields in PML (i = i+1/2)
	// X-direction
	for (i = 0; i<PMLX; i++)
	{
		bi = PMLX - 1 - i; //recursive index
		Dxa1[i] = (1 - (sigma_h_x1[bi] * dt) / (2 * eps_0)) / (1 + (sigma_h_x1[bi] * dt) / (2 * eps_0));
		Dxb1[i] = ((sigma_h_x1[bi] * dt) / (eps_0*dx[0])) / (1 + (sigma_h_x1[bi] * dt) / (2 * eps_0));  //X-

		Dxa2[i] = (1 - (sigma_h_x2[i] * dt) / (2 * eps_0)) / (1 + (sigma_h_x2[i] * dt) / (2 * eps_0));
		Dxb2[i] = ((sigma_h_x2[i] * dt) / (eps_0*dx[txsize-1])) / (1 + (sigma_h_x2[i] * dt) / (2 * eps_0));    //X+
	}
	memfree(&sigma_e_x1);
	memfree(&sigma_e_x2);
	memfree(&sigma_h_x1);
	memfree(&sigma_h_x2);
}

void TFDTD_2PML::PML_parameter_initiateY()
{
	double eps_r4pml = 1.0;
	double mu_r = 1;
	double m = 4.0;
	double eta_0 = sqrt(mu_0 / eps_0);

	double sigma_max1 = 0.8 * (m + 1.0) / (eta_0 * dy[0] * sqrt(mu_r * eps_r4pml)); // maximum conductivity
	double sigma_max2 = 0.8 * (m + 1.0) / (eta_0 * dy[tysize - 1] * sqrt(mu_r * eps_r4pml)); // maximum conductivity
	double* sigma_e_y1 = Tptr1<double>(PMLY + 1);
	double* sigma_e_y2 = Tptr1<double>(PMLY + 1);
	double* sigma_h_y1 = Tptr1<double>(PMLY + 1);
	double* sigma_h_y2 = Tptr1<double>(PMLY + 1);
	// Degree of polynomial

	for (j = 0; j<PMLY + 1; j++)
	{
		sigma_e_y1[i] = sigma_max1 * pow((i + 0.0) / PMLY, m);
		sigma_e_y2[i] = sigma_max2 * pow((i + 0.0) / PMLY, m);
		sigma_h_y1[i] = sigma_max1 * pow((i + 1.0 / 2.0) / PMLY, m);
		sigma_h_y2[i] = sigma_max2 * pow((i + 1.0 / 2.0) / PMLY, m);
	}

	// setup kappa_e, Ca,Cb for updating E fields in PML
	// Y-direction
	for (j = 0; j<PMLY + 1; j++)
	{
		aj = PMLY - j; //recursive index
		Cya1[j] = (1 - (sigma_e_y1[aj] * dt) / (2 * eps_0)) / (1 + (sigma_e_y1[aj] * dt) / (2 * eps_0));
		Cyb1[j] = ((sigma_e_y1[aj] * dt) / (eps_0*dy[0])) / (1 + (sigma_e_y1[aj] * dt) / (2 * eps_0));  //Y-

		Cya2[j] = (1 - (sigma_e_y2[j] * dt) / (2 * eps_0)) / (1 + (sigma_e_y2[j] * dt) / (2 * eps_0));
		Cyb2[j] = ((sigma_e_y2[j] * dt) / (eps_0*dy[tysize-1])) / (1 + (sigma_e_y2[j] * dt) / (2 * eps_0));    //Y+
	}
	// setup kappa_h, Da,Db for updating H fields in PML (i = i+1/2)
	// Y-direction
	for (j = 0; j<PMLY; j++)
	{
		bj = PMLY - 1 - j; //recursive index
		Dya1[j] = (1 - (sigma_h_y1[bj] * dt) / (2 * eps_0)) / (1 + (sigma_h_y1[bj] * dt) / (2 * eps_0));
		Dyb1[j] = ((sigma_h_y1[bj] * dt) / (eps_0*dy[0])) / (1 + (sigma_h_y1[bj] * dt) / (2 * eps_0));  //Y-

		Dya2[j] = (1 - (sigma_h_y2[j] * dt) / (2 * eps_0)) / (1 + (sigma_h_y2[j] * dt) / (2 * eps_0));
		Dyb2[j] = ((sigma_h_y2[j] * dt) / (eps_0*dy[tysize-1])) / (1 + (sigma_h_y2[j] * dt) / (2 * eps_0));    //Y+
	}
	memfree(&sigma_e_y1);
	memfree(&sigma_h_y1);
	memfree(&sigma_e_y2);
	memfree(&sigma_h_y2);
}

void TFDTD_2PML::PML_parameter_initiateZ()
{
	double eps_r4pml = 1.0;
	double mu_r = 1;
	double m = 4.0;
	double eta_0 = sqrt(mu_0 / eps_0);

	double sigma_max1 = 0.8 * (m + 1.0) / (eta_0 * dz[0] * sqrt(mu_r * eps_r4pml)); // maximum conductivity
	double sigma_max2 = 0.8 * (m + 1.0) / (eta_0 * dz[tzsize - 1] * sqrt(mu_r * eps_r4pml)); // maximum conductivity
	double* sigma_e_z1 = Tptr1<double>(PMLZ + 1); // profile in z-direction
	double* sigma_e_z2 = Tptr1<double>(PMLZ + 1); // profile in z-direction
	double* sigma_h_z1 = Tptr1<double>(PMLZ + 1);
	double* sigma_h_z2 = Tptr1<double>(PMLZ + 1);

	for (k = 0; k<PMLZ + 1; k++)
	{
		sigma_e_z1[k] = sigma_max1 * pow((k + 0.0) / PMLZ, m);
		sigma_h_z1[k] = sigma_max1 * pow((k + 1.0 / 2.0) / PMLZ, m);
		sigma_e_z2[k] = sigma_max2 * pow((k + 0.0) / PMLZ, m);
		sigma_h_z2[k] = sigma_max2 * pow((k + 1.0 / 2.0) / PMLZ, m);
	}

	// setup kappa_e, Ca,Cb for updating E fields in PML
	// Z-direction
	for (k = 0; k<PMLZ + 1; k++)
	{
		ak = PMLZ - k; //recursive index
		Cza1[k] = (1 - (sigma_e_z1[ak] * dt) / (2 * eps_0)) / (1 + (sigma_e_z1[ak] * dt) / (2 * eps_0));
		Czb1[k] = ((sigma_e_z1[ak] * dt) / (eps_0*dz[0])) / (1 + (sigma_e_z1[ak] * dt) / (2 * eps_0));  //Z-

		Cza2[k] = (1 - (sigma_e_z2[k] * dt) / (2 * eps_0)) / (1 + (sigma_e_z2[k] * dt) / (2 * eps_0));
		Czb2[k] = ((sigma_e_z2[k] * dt) / (eps_0*dz[tzsize-1])) / (1 + (sigma_e_z2[k] * dt) / (2 * eps_0));    //Z+
	}

	// setup kappa_h, Da,Db for updating H fields in PML (i = i+1/2)
	// Z-direction
	for (k = 0; k<PMLZ; k++)
	{
		bk = PMLZ - 1 - k; //recursive index
		Dza1[k] = (1 - (sigma_h_z1[bk] * dt) / (2 * eps_0)) / (1 + (sigma_h_z1[bk] * dt) / (2 * eps_0));
		Dzb1[k] = ((sigma_h_z1[bk] * dt) / (eps_0*dz[0])) / (1 + (sigma_h_z1[bk] * dt) / (2 * eps_0));  //Z-

		Dza2[k] = (1 - (sigma_h_z2[k] * dt) / (2 * eps_0)) / (1 + (sigma_h_z2[k] * dt) / (2 * eps_0));
		Dzb2[k] = ((sigma_h_z2[k] * dt) / (eps_0*dz[tzsize-1])) / (1 + (sigma_h_z2[k] * dt) / (2 * eps_0));    //Z+
	}
	// meme_Free
	memfree(&sigma_e_z1);
	memfree(&sigma_h_z1);
	memfree(&sigma_e_z2);
	memfree(&sigma_h_z2);
}

void TFDTD_2PML::X_NFpml()
{
	// Fy, Fz
	//  (X-)
	int i, j, k;
	int mx, my, mz;
	double temp;
#pragma omp parallel private(i, j, k, temp, mx, my, mz)
	{
#pragma omp for schedule(static) collapse(3)
		for (i = 1; i < PMLX + 1; i++)
		for (j = 0; j < tysize; j++)
		for (k = 1; k < tzsize; k++)
		{
			my = Media_Ey[i][j][k];
			temp = fyx1[i][j][k];
			fyx1[i][j][k] = Cxa1[i] * fyx1[i][j][k] + Cxb1[i] * (Hz[i][j][k] - Hz[i - 1][j][k]);
			Ey[i][j][k] = Ey[i][j][k] + C_d[my] * (temp + fyx1[i][j][k]);
		}

		//  (X-)
#pragma omp for schedule(static) collapse(3)
		for (i = 1; i < PMLX + 1; i++)
		for (j = 1; j < tysize; j++)
		for (k = 0; k < tzsize; k++)
		{
			mz = Media_Ez[i][j][k];
			temp = fzx1[i][j][k];
			fzx1[i][j][k] = Cxa1[i] * fzx1[i][j][k] + Cxb1[i] * (Hy[i][j][k] - Hy[i - 1][j][k]);
			Ez[i][j][k] = Ez[i][j][k] - C_d[mz] * (temp + fzx1[i][j][k]);
		}
	}
}

void TFDTD_2PML::X_NGpml()
{
	// Gy, Gz
	//   (X-)
	int i, j, k;
	double temp;
#pragma omp parallel private(i,j,k, temp)
	{
#pragma omp for schedule(static) collapse(3)
		for (i = 0; i < PMLX; i++)
		for (j = 1; j < tysize; j++)
		for (k = 0; k < tzsize; k++)
		{
			temp = gyx1[i][j][k];
			gyx1[i][j][k] = Dxa1[i] * gyx1[i][j][k] + Dxb1[i] * (Ez[i + 1][j][k] - Ez[i][j][k]);
			Hy[i][j][k] = Hy[i][j][k] - C_h * (temp + gyx1[i][j][k]);
		}
		//   (X-)
#pragma omp for schedule(static) collapse(3)
		for (i = 0; i < PMLX; i++)
		for (j = 0; j < tysize; j++)
		for (k = 1; k < tzsize; k++)
		{
			temp = gzx1[i][j][k];
			gzx1[i][j][k] = Dxa1[i] * gzx1[i][j][k] + Dxb1[i] * (Ey[i + 1][j][k] - Ey[i][j][k]);
			Hz[i][j][k] = Hz[i][j][k] + C_h * (temp + gzx1[i][j][k]);
		}
	}
}

void TFDTD_2PML::X_PFpml()
{
	//   (X+)
	int ai, i, j, k;
	int mx, my, mz;
	double temp;
#pragma omp parallel private(ai, i,j,k, temp, mx, my, mz)
	{
#pragma omp for schedule(static) 
		for (ai = 0; ai < PMLX; ai++)
		{
			i = txsize - PMLX + ai;
			for (j = 0; j < tysize; j++)
			for (k = 1; k < tzsize; k++)
			{
				my = Media_Ey[i][j][k];
				temp = fyx2[ai][j][k];
				fyx2[ai][j][k] = Cxa2[ai] * fyx2[ai][j][k] + Cxb2[ai] * (Hz[i][j][k] - Hz[i - 1][j][k]);
				Ey[i][j][k] = Ey[i][j][k] + C_d[my] * (temp + fyx2[ai][j][k]);
			}
		}

		//   (X+)
#pragma omp for schedule(static)
		for (ai = 0; ai < PMLX; ai++)
		{
			i = txsize - PMLX + ai;
			for (j = 1; j < tysize; j++)
			for (k = 0; k < tzsize; k++)
			{
				mz = Media_Ez[i][j][k];
				temp = fzx2[ai][j][k];
				fzx2[ai][j][k] = Cxa2[ai] * fzx2[ai][j][k] + Cxb2[ai] * (Hy[i][j][k] - Hy[i - 1][j][k]);
				Ez[i][j][k] = Ez[i][j][k] - C_d[mz] * (temp + fzx2[ai][j][k]);
			}
		}
	}
}
void TFDTD_2PML::X_PGpml()
{
	int bi, i, j, k;
	double temp;
#pragma omp parallel private(bi, i,j,k, temp)
	{
		//   (X+)
#pragma omp for schedule(static)
		for (bi = 0; bi < PMLX; bi++)
		{
			i = txsize - PMLX + bi;
			for (j = 1; j < tysize; j++)
			for (k = 0; k < tzsize; k++)
			{
				temp = gyx2[bi][j][k];
				gyx2[bi][j][k] = Dxa2[bi] * gyx2[bi][j][k] + Dxb2[bi] * (Ez[i + 1][j][k] - Ez[i][j][k]);
				Hy[i][j][k] = Hy[i][j][k] - C_h * (temp + gyx2[bi][j][k]);
			}
		}

		//   (X+)
#pragma omp for schedule(static)
		for (bi = 0; bi < PMLX; bi++)
		{
			i = txsize - PMLX + bi;
			for (j = 0; j < tysize; j++)
			for (k = 1; k < tzsize; k++)
			{
				temp = gzx2[bi][j][k];
				gzx2[bi][j][k] = Dxa2[bi] * gzx2[bi][j][k] + Dxb2[bi] * (Ey[i + 1][j][k] - Ey[i][j][k]);
				Hz[i][j][k] = Hz[i][j][k] + C_h * (temp + gzx2[bi][j][k]);

			}
		}
	}

}

void TFDTD_2PML::Y_NFpml()
{
	int i, j, k;
	int mx, my, mz;
	double temp;
#pragma omp parallel private(i,j,k,temp,mx,my,mz)
	{
		//   (Y-)
#pragma omp for schedule(static) collapse(3)
		for (i = 0; i < txsize; i++)
		for (j = 1; j < PMLY + 1; j++)
		for (k = 1; k < tzsize; k++)
		{
			mx = Media_Ex[i][j][k];
			temp = fxy1[i][j][k];
			fxy1[i][j][k] = Cya1[j] * fxy1[i][j][k] + Cyb1[j] * (Hz[i][j][k] - Hz[i][j - 1][k]);
			Ex[i][j][k] = Ex[i][j][k] - C_d[mx] * (temp + fxy1[i][j][k]);
		}
		//  (Y-)
#pragma omp for schedule(static) collapse(3)
		for (i = 1; i < txsize; i++)
		for (j = 1; j < PMLY + 1; j++)
		for (k = 0; k < tzsize; k++)
		{
			mz = Media_Ez[i][j][k];
			temp = fzy1[i][j][k];
			fzy1[i][j][k] = Cya1[j] * fzy1[i][j][k] + Cyb1[j] * (Hx[i][j][k] - Hx[i][j - 1][k]);
			Ez[i][j][k] = Ez[i][j][k] + C_d[mz] * (temp + fzy1[i][j][k]);
		}
	}
}
void TFDTD_2PML::Y_NGpml()
{
	int i, j, k;
	double temp;
#pragma omp parallel private(i,j,k, temp)
	{
		//   (Y-)
#pragma omp for schedule(static) collapse(3)
		for (i = 1; i < txsize; i++)
		for (j = 0; j < PMLY; j++)
		for (k = 0; k < tzsize; k++)
		{
			temp = gxy1[i][j][k];
			gxy1[i][j][k] = Dya1[j] * gxy1[i][j][k] + Dyb1[j] * (Ez[i][j + 1][k] - Ez[i][j][k]);
			Hx[i][j][k] = Hx[i][j][k] + C_h * (temp + gxy1[i][j][k]);
		}

		//   (Y-)
#pragma omp for schedule(static) collapse(3)
		for (i = 0; i < txsize; i++)
		for (j = 0; j < PMLY; j++)
		for (k = 1; k < tzsize; k++)
		{
			temp = gzy1[i][j][k];
			gzy1[i][j][k] = Dya1[j] * gzy1[i][j][k] + Dyb1[j] * (Ex[i][j + 1][k] - Ex[i][j][k]);
			Hz[i][j][k] = Hz[i][j][k] - C_h * (temp + gzy1[i][j][k]);
		}
	}
}

void TFDTD_2PML::Y_PFpml()
{
	int i, aj, j, k;
	int mx, my, mz;
	double temp;
#pragma omp parallel private(i,aj,j,k,temp,mx,my,mz)
	{
		//   (Y+)
#pragma omp for schedule(static) collapse(2)
		for (i = 0; i < txsize; i++)
		for (aj = 0; aj < PMLY; aj++)
		{
			j = tysize - (PMLY)+aj;
			for (k = 1; k < tzsize; k++)
			{
				mx = Media_Ex[i][j][k];
				temp = fxy2[i][aj][k];
				fxy2[i][aj][k] = Cya2[aj] * fxy2[i][aj][k] + Cyb2[aj] * (Hz[i][j][k] - Hz[i][j - 1][k]);
				Ex[i][j][k] = Ex[i][j][k] - C_d[mx] * (temp + fxy2[i][aj][k]);
			}
		}


		// (Y+)
#pragma omp for schedule(static) collapse(2)
		for (i = 1; i < txsize; i++)
		for (aj = 0; aj < PMLY; aj++)
		{
			j = tysize - (PMLY)+aj;
			for (k = 0; k < tzsize; k++)
			{
				mz = Media_Ez[i][j][k];
				temp = fzy2[i][aj][k];
				fzy2[i][aj][k] = Cya2[aj] * fzy2[i][aj][k] + Cyb2[aj] * (Hx[i][j][k] - Hx[i][j - 1][k]);
				Ez[i][j][k] = Ez[i][j][k] + C_d[mz] * (temp + fzy2[i][aj][k]);
			}
		}
	}
}

void TFDTD_2PML::Y_PGpml()
{
	int i, bj, j, k;
	double temp;
#pragma omp parallel private(i,bj,j,k, temp)
	{
		//   (Y+)
#pragma omp for schedule(static) collapse(2)
		for (i = 1; i < txsize; i++)
		for (bj = 0; bj < PMLY; bj++)
		{
			j = tysize - PMLY + bj;
			for (k = 0; k < tzsize; k++)
			{
				temp = gxy2[i][bj][k];
				gxy2[i][bj][k] = Dya2[bj] * gxy2[i][bj][k] + Dyb2[bj] * (Ez[i][j + 1][k] - Ez[i][j][k]);
				Hx[i][j][k] = Hx[i][j][k] + C_h * (temp + gxy2[i][bj][k]);
			}
		}

		//   (Y+)
#pragma omp for schedule(static) collapse(2)
		for (i = 0; i < txsize; i++)
		for (bj = 0; bj < PMLY; bj++)
		{
			j = tysize - PMLY + bj;
			for (k = 1; k < tzsize; k++)
			{
				temp = gzy2[i][bj][k];
				gzy2[i][bj][k] = Dya2[bj] * gzy2[i][bj][k] + Dyb2[bj] * (Ex[i][j + 1][k] - Ex[i][j][k]);
				Hz[i][j][k] = Hz[i][j][k] - C_h * (temp + gzy2[i][bj][k]);
			}
		}
	}
}

void TFDTD_2PML::Z_NFpml()
{
	int i, j, k;
	int mx, my, mz;
	double temp;
#pragma omp parallel private(i,j,k,temp,mx,my,mz)
	{
		//  (Z-)
#pragma omp for schedule(static) collapse(3)
		for (i = 0; i < txsize; i++)
		for (j = 1; j < tysize; j++)
		for (k = 1; k < PMLZ + 1; k++)
		{
			mx = Media_Ex[i][j][k];
			temp = fxz1[i][j][k];
			fxz1[i][j][k] = Cza1[k] * fxz1[i][j][k] + Czb1[k] * (Hy[i][j][k] - Hy[i][j][k - 1]);
			Ex[i][j][k] = Ex[i][j][k] + C_d[mx] * (temp + fxz1[i][j][k]);
		}

		//   (Z-)
#pragma omp for schedule(static) collapse(3)
		for (i = 1; i < txsize; i++)
		for (j = 0; j < tysize; j++)
		for (k = 1; k < PMLZ + 1; k++)
		{
			my = Media_Ey[i][j][k];
			temp = fyz1[i][j][k];
			fyz1[i][j][k] = Cza1[k] * fyz1[i][j][k] + Czb1[k] * (Hx[i][j][k] - Hx[i][j][k - 1]);
			Ey[i][j][k] = Ey[i][j][k] - C_d[my] * (temp + fyz1[i][j][k]);
		}
	}
}
void TFDTD_2PML::Z_NGpml()
{
	int i, j, k;
	double temp;
#pragma omp parallel private(i,j,k,temp)
	{
	//   (Z-)
#pragma omp for schedule(static) collapse(3)
		for (i = 1; i<txsize; i++) 
		for (j = 0; j<tysize; j++) 
		for (k = 0; k<PMLZ; k++) 
		{
			temp = gxz1[i][j][k];
			gxz1[i][j][k] = Dza1[k] * gxz1[i][j][k] + Dzb1[k] * (Ey[i][j][k + 1] - Ey[i][j][k]);
			Hx[i][j][k] = Hx[i][j][k] - C_h*(temp + gxz1[i][j][k]);
		}

	//   (Z-)
#pragma omp for schedule(static) collapse(3)
		for (i = 0; i<txsize; i++) 
		for (j = 1; j<tysize; j++) 
		for (k = 0; k<PMLZ; k++) 
		{
			temp = gyz1[i][j][k];
			gyz1[i][j][k] = Dza1[k] * gyz1[i][j][k] + Dzb1[k] * (Ex[i][j][k + 1] - Ex[i][j][k]);
			Hy[i][j][k] = Hy[i][j][k] + C_h*(temp + gyz1[i][j][k]);
		}
	}
}

void TFDTD_2PML::Z_PFpml()
{
	int i, j, ak, k;
	int mx, my, mz;
	double temp;
#pragma omp parallel private(i,j,ak,k, temp, mx, my,mz)
	{
		//   (Z+)
#pragma omp for schedule(static) collapse(3)
		for (i = 0; i < txsize; i++)
		for (j = 1; j < tysize; j++)
		for (ak = 0; ak < PMLZ; ak++)
		{
			k = tzsize - (PMLZ)+ak;
			mx = Media_Ex[i][j][k];
			temp = fxz2[i][j][ak];
			fxz2[i][j][ak] = Cza2[ak] * fxz2[i][j][ak] + Czb2[ak] * (Hy[i][j][k] - Hy[i][j][k - 1]);
			Ex[i][j][k] = Ex[i][j][k] + C_d[mx] * (temp + fxz2[i][j][ak]);
		}

		//   (Z+)
#pragma omp for schedule(static) collapse(3)
		for (i = 1; i < txsize; i++)
		for (j = 0; j < tysize; j++)
		for (ak = 0; ak < PMLZ; ak++)
		{
			k = tzsize - (PMLZ)+ak;
			my = Media_Ey[i][j][k];
			temp = fyz2[i][j][ak];
			fyz2[i][j][ak] = Cza2[ak] * fyz2[i][j][ak] + Czb2[ak] * (Hx[i][j][k] - Hx[i][j][k - 1]);
			Ey[i][j][k] = Ey[i][j][k] - C_d[my] * (temp + fyz2[i][j][ak]);
		}
	}
}

void TFDTD_2PML::Z_PGpml()
{
	int i, j, bk, k;
	double temp;
#pragma omp parallel private(i,j,bk,k,temp)
	{
		//   (Z+)
#pragma omp for schedule(static) collapse(3)
		for (i = 1; i < txsize; i++)
		for (j = 0; j < tysize; j++)
		for (bk = 0; bk < PMLZ; bk++)
		{
			k = tzsize - PMLZ + bk;
			temp = gxz2[i][j][bk];
			gxz2[i][j][bk] = Dza2[bk] * gxz2[i][j][bk] + Dzb2[bk] * (Ey[i][j][k + 1] - Ey[i][j][k]);
			Hx[i][j][k] = Hx[i][j][k] - C_h * (temp + gxz2[i][j][bk]);
		}

		//   (Z+)
#pragma omp for schedule(static) collapse(3)
		for (i = 0; i < txsize; i++)
		for (j = 1; j < tysize; j++)
		for (bk = 0; bk < PMLZ; bk++)
		{
			k = tzsize - PMLZ + bk;
			temp = gyz2[i][j][bk];
			gyz2[i][j][bk] = Dza2[bk] * gyz2[i][j][bk] + Dzb2[bk] * (Ex[i][j][k + 1] - Ex[i][j][k]);
			Hy[i][j][k] = Hy[i][j][k] + C_h * (temp + gyz2[i][j][bk]);
		}
	}
}


