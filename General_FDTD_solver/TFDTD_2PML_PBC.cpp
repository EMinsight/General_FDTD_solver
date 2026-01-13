//---------------------------------------------------------------------------
#pragma hdrstop
#include "TFDTD_2PML_PBC.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
TFDTD_2PML_PBC::~TFDTD_2PML_PBC(){}

void TFDTD_2PML_PBC::PBC_Z_F_PML_minus()
{
	int i, j, k, mx, my;
	int ai, aj, ak;
	double temp;
	//  (Z-)
#pragma omp parallel private(i,j,k, temp, ak, mx, my)
	{
#pragma omp for schedule(dynamic)
		for (i = 0; i < txsize; i++) 
		for (j = 0; j < tysize; j++) 
		for (k = 1; k < PMLZ + 1; k++)
		{
			temp = fxz1[i][j][k];
			mx=Media_Ex[i][j][k];
			fxz1[i][j][k] = Cza1[k] * fxz1[i][j][k] + Czb1[k] * (Hy[i][j][k] - Hy[i][j][k - 1]);
			Ex[i][j][k]= Ex[i][j][k] + C_d[mx]*(temp + fxz1[i][j][k]);
		}

#pragma omp for schedule(dynamic)
		for (i = 0; i < txsize; i++) 
		for (j = 0; j < tysize; j++) 
		for (k = 1; k < PMLZ + 1; k++)
		{
			temp = fyz1[i][j][k];
			my=Media_Ey[i][j][k];
			fyz1[i][j][k] = Cza1[k] * fyz1[i][j][k] + Czb1[k] * (Hx[i][j][k] - Hx[i][j][k - 1]);
			Ey[i][j][k] = Ey[i][j][k] - C_d[my]*(temp + fyz1[i][j][k]);
		}
	}
}

void TFDTD_2PML_PBC::PBC_Z_G_PML_minus()
{
	int i, j, k;
	int bi, bj, bk;
	double temp;
	//   (Z-)
#pragma omp parallel private(i,j,k, temp, bk)
	{
#pragma omp for schedule(dynamic)
		for (i = 0; i < txsize; i++) 
		for (j = 0; j < tysize; j++) 
		for (k = 0; k < PMLZ; k++) 
		{
			temp = gxz1[i][j][k];
			gxz1[i][j][k] = Dza1[k] * gxz1[i][j][k] + Dzb1[k] * (Ey[i][j][k + 1] - Ey[i][j][k]);
			Hx[i][j][k] = Hx[i][j][k] - C_h*(temp + gxz1[i][j][k]);
		}
		//   (Z-)
#pragma omp for schedule(dynamic)
		for (i = 0; i < txsize; i++) 
		for (j = 0; j < tysize; j++) 
		for (k = 0; k < PMLZ; k++) 
		{
			temp = gyz1[i][j][k];
			gyz1[i][j][k] = Dza1[k] * gyz1[i][j][k] + Dzb1[k] * (Ex[i][j][k + 1] - Ex[i][j][k]);
			Hy[i][j][k] = Hy[i][j][k] + C_h*(temp + gyz1[i][j][k]);
		}
	}
}

void TFDTD_2PML_PBC::PBC_Z_F_PML_plus()
{
	//   (Z+)
	int i, j, k, mx, my;
	int ai, aj, ak;
	double temp;
	//  (Z-)
#pragma omp parallel private(i,j,k, temp, ak, mx, my)
	{
#pragma omp for schedule(dynamic)
		for (i = 0;   i < txsize; i++)
		for (j = 0;   j < tysize; j++)
		for (ak = 0; ak < PMLZ; ak++)
		{
			k = tzsize - PMLZ + ak;
			temp = fxz2[i][j][ak];
			mx=Media_Ex[i][j][k];
			fxz2[i][j][ak] = Cza2[ak] * fxz2[i][j][ak] + Czb2[ak] * (Hy[i][j][k] - Hy[i][j][k - 1]);
			Ex[i][j][k] = Ex[i][j][k] + C_d[mx]*(temp + fxz2[i][j][ak]);
		}

		//   (Z+)
#pragma omp for schedule(dynamic)
		for (i = 0; i < txsize; i++)
		for (j = 0; j < tysize; j++)
		for (ak = 0; ak<PMLZ; ak++)
		{
			k = tzsize - PMLZ + ak;
			temp = fyz2[i][j][ak];
			my=Media_Ey[i][j][k];
			fyz2[i][j][ak] = Cza2[ak] * fyz2[i][j][ak] + Czb2[ak] * (Hx[i][j][k] - Hx[i][j][k - 1]);
			Ey[i][j][k] = Ey[i][j][k] - C_d[my]*(temp + fyz2[i][j][ak]);
		}
	}
}

void TFDTD_2PML_PBC::PBC_Z_G_PML_plus()
{
	int i, j, k;
	int bi, bj, bk;
	double temp;
	//   (Z-)
#pragma omp parallel private(i,j,k, temp, bk)
	{
		//   (Z+)
#pragma omp for schedule(dynamic)
		for (i = 0; i < txsize; i++)
		for (j = 0; j < tysize; j++)
		for (bk = 0; bk < PMLZ; bk++)
		{
			k = tzsize - PMLZ + bk;
			temp = gxz2[i][j][bk];
			gxz2[i][j][bk] = Dza2[bk] * gxz2[i][j][bk] + Dzb2[bk] * (Ey[i][j][k + 1] - Ey[i][j][k]);
			Hx[i][j][k] = Hx[i][j][k] - C_h*(temp + gxz2[i][j][bk]);
		}
		//   (Z+)
#pragma omp for schedule(dynamic)
		for (i = 0; i < txsize; i++)
		for (j = 0; j < tysize; j++)
		for (bk = 0; bk < PMLZ; bk++)
		{
			k = tzsize - PMLZ + bk;
			temp = gyz2[i][j][bk];
			gyz2[i][j][bk] = Dza2[bk] * gyz2[i][j][bk] + Dzb2[bk] * (Ex[i][j][k + 1] - Ex[i][j][k]);
			Hy[i][j][k] = Hy[i][j][k] + C_h*(temp + gyz2[i][j][bk]);
		}
	}
}
