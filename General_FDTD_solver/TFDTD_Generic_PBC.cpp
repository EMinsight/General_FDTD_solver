//---------------------------------------------------------------------------
#pragma hdrstop
#include "TFDTD_Generic_PBC.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)

TFDTD_Generic_PBC::TFDTD_Generic_PBC()
{
	eps_0 = 1e-9 / 36 / m_pi;
	mu_0 = 4 * m_pi * 1e-7;
	c_0 = 1 / sqrt(eps_0 * mu_0);
}

TFDTD_Generic_PBC::~TFDTD_Generic_PBC() {}

void TFDTD_Generic_PBC::set_PBC_Screwed_Parameters()
{
	//Read_Media_Data();
	//dt = 0.99 / c_0 / sqrt((1 / min1D(dx)) * (1 / min1D(dx)) + (1 / min1D(dy)) * (1 / min1D(dy)) + (1 / min1D(dz)) * (1 / min1D(dz)));
	//Initiate_MediaSpace();
	//Basic_parameter_Initiate();

	FILE* pf1;
	pf1 = fopen("CHIRALPBC.dat", "r");
	fscanf(pf1, "%d,%d\n", &aa1m1, &aa1m2);
	fclose(pf1);
	printf("@FieldUpdate\naa1m1=%d\naa1m2=%d\n", aa1m1, aa1m2);
}

void TFDTD_Generic_PBC::PBC_Save_cal()
{
	int i = 0, j = 0, k = 0, pol;
	int mx, my, mz;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i = 0; i < txsize; i++)
			for (j = 0; j < tysize + 1; j++)
				//for (k = 0; k < tzsize; k++)
			{
				memcpy(Ex_pp[i][j], Ex_p[i][j], sizeof(double) * (tzsize+1));
				memcpy(Ex_p[i][j], Ex[i][j], sizeof(double) * (tzsize + 1));
				// Ex_pp[i][j][k] = Ex_p[i][j][k];
				// Ex_p[i][j][k] = Ex[i][j][k];
			}
#pragma omp for schedule(dynamic) 
		for (i = 0; i < txsize + 1; i++)
			for (j = 0; j < tysize; j++)
				//for (k = 0; k < tzsize; k++)
			{
				memcpy(Ey_pp[i][j], Ey_p[i][j], sizeof(double) * (tzsize + 1));
				memcpy(Ey_p[i][j], Ey[i][j], sizeof(double) * (tzsize + 1));
				//Ey_pp[i][j][k] = Ey_p[i][j][k];
				//Ey_p[i][j][k] = Ey[i][j][k];
			}
#pragma omp for schedule(dynamic)
		for (i = 0; i < txsize + 1; i++)
			for (j = 0; j < tysize + 1; j++)
				//for (k = 0; k < tzsize; k++)
			{
				memcpy(Ez_pp[i][j], Ez_p[i][j], sizeof(double) * tzsize);
				memcpy(Ez_p[i][j], Ez[i][j], sizeof(double) * tzsize);
				//Ez_pp[i][j][k] = Ez_p[i][j][k];
				//Ez_p[i][j][k] = Ez[i][j][k];
			}
	}
}

void TFDTD_Generic_PBC::E_field_update()
{
	//n_Thread = omp_get_num_procs(); // omp_get_max_thread();
	//omp_set_num_threads(n_Thread);
	int i = 0, j = 0, k = 0, pol = 0;
	int mx, my, mz;
	double temp1 = 0, temp2 = 0, temp3 = 0, temp4 = 0, temp5 = 0, temp6 = 0;
#pragma omp parallel private(i,j,k, pol, temp1, temp2, temp3, temp4, temp5, temp6, mx, my, mz)
	{
#pragma omp for schedule(dynamic)
		for (i = 0; i < txsize - 1; i++)
			for (j = 1; j < tysize - 1; j++)
				for (k = 1; k < tzsize; k++)
				{
					mx = Media_Ex[i][j][k];
					temp1 = 0.; temp2 = 0.;
					if (pMp[mx].TypeOfMaterials == 'D')
					{
						for (pol = 0; pol < total_pol; pol++)
						{
							temp1 += Ca_temp[mx][pol] * Jx[i][j][k][pol];
							temp2 += Cb_temp[mx][pol] * Jx_p[i][j][k][pol];
						}
					}
					Ex[i][j][k] = C1[mx] * Ex_p[i][j][k] - C2[mx] * Ex_pp[i][j][k]
						+ Cb[mx] / E_dy[j] * (Hz[i][j][k] - Hz[i][j - 1][k])
						- Cb[mx] / E_dz[k] * (Hy[i][j][k] - Hy[i][j][k - 1]) - temp1 - temp2;
				}

#pragma omp for schedule(dynamic)
		for (i = 1; i < txsize - 1; i++)
			for (j = 0; j < tysize - 1; j++)
				for (k = 1; k < tzsize; k++)
				{
					my = Media_Ey[i][j][k];
					temp3 = 0.; temp4 = 0.;
					if (pMp[my].TypeOfMaterials == 'D')
					{
						for (pol = 0; pol < total_pol; pol++)
						{
							temp3 += Ca_temp[my][pol] * Jy[i][j][k][pol];
							temp4 += Cb_temp[my][pol] * Jy_p[i][j][k][pol];
						}
					}

					Ey[i][j][k] = C1[my] * Ey_p[i][j][k] - C2[my] * Ey_pp[i][j][k]
						+ Cb[my] / E_dz[k] * (Hx[i][j][k] - Hx[i][j][k - 1])
						- Cb[my] / E_dx[i] * (Hz[i][j][k] - Hz[i - 1][j][k]) - temp3 - temp4;

				}
#pragma omp for schedule(dynamic)
		for (i = 1; i < txsize - 1; i++)
			for (j = 1; j < tysize - 1; j++)
				for (k = 0; k < tzsize; k++)
				{
					mz = Media_Ez[i][j][k];
					temp5 = 0.; temp6 = 0.;
					if (pMp[mz].TypeOfMaterials == 'D')
					{
						for (pol = 0; pol < total_pol; pol++)
						{
							temp5 += Ca_temp[mz][pol] * Jz[i][j][k][pol];
							temp6 += Cb_temp[mz][pol] * Jz_p[i][j][k][pol];
						}
					}
					Ez[i][j][k] = C1[mz] * Ez_p[i][j][k] - C2[mz] * Ez_pp[i][j][k]
						+ Cb[mz] / E_dx[i] * (Hy[i][j][k] - Hy[i - 1][j][k])
						- Cb[mz] / E_dy[j] * (Hx[i][j][k] - Hx[i][j - 1][k]) - temp5 - temp6;
				}
	}
}

void TFDTD_Generic_PBC::J_field_update()
{
	int i, j, k, pol;
	int mx, my, mz;
	double tempx, tempy, tempz;

#pragma omp parallel private(i,j,k,pol,mx, tempx)
	{
		//#pragma omp for schedule(dynamic)
#pragma omp for schedule(dynamic)
		for (i = 0; i < txsize; i++)
			for (j = 0; j < tysize + 1; j++)
				for (k = 0; k < tzsize; k++)
				{
					mx = Media_Ex[i][j][k];
					if (pMp[mx].TypeOfMaterials == 'D')
					{
						for (pol = 0; pol < total_pol; pol++)
						{
							tempx = Jx[i][j][k][pol];
							Jx[i][j][k][pol] = pMp[mx].Ca[pol] * Jx[i][j][k][pol] + pMp[mx].Cc[pol] * Ex[i][j][k] +
								pMp[mx].Cb[pol] * Jx_p[i][j][k][pol] + pMp[mx].Cd[pol] * Ex_p[i][j][k]
								+ pMp[mx].Ce[pol] * Ex_pp[i][j][k];
							Jx_p[i][j][k][pol] = tempx;
						}
					}
				}
	}

#pragma omp parallel private(i,j,k,pol,my, tempy)
	{
#pragma omp for schedule(dynamic)
		for (i = 0; i < txsize + 1; i++)
			for (j = 0; j < tysize; j++)
				for (k = 0; k < tzsize; k++)
				{
					my = Media_Ey[i][j][k];
					if (pMp[my].TypeOfMaterials == 'D')
					{
						for (pol = 0; pol < total_pol; pol++)
						{
							tempy = Jy[i][j][k][pol];
							Jy[i][j][k][pol] = pMp[my].Ca[pol] * Jy[i][j][k][pol] + pMp[my].Cc[pol] * Ey[i][j][k] +
								pMp[my].Cb[pol] * Jy_p[i][j][k][pol] + pMp[my].Cd[pol] * Ey_p[i][j][k]
								+ pMp[my].Ce[pol] * Ey_pp[i][j][k];
							Jy_p[i][j][k][pol] = tempy;
						}
					}
				}
	}

#pragma omp parallel private(i,j,k, pol,mz, tempz)
	{
#pragma omp for schedule(dynamic)
		for (i = 0; i < txsize + 1; i++)
			for (j = 0; j < tysize + 1; j++)
				for (k = 0; k < tzsize; k++)
				{
					mz = Media_Ez[i][j][k];
					if (pMp[mz].TypeOfMaterials == 'D')
					{
						for (pol = 0; pol < total_pol; pol++)
						{
							tempz = Jz[i][j][k][pol];
							Jz[i][j][k][pol] = pMp[mz].Ca[pol] * Jz[i][j][k][pol] + pMp[mz].Cc[pol] * Ez[i][j][k] +
								pMp[mz].Cb[pol] * Jz_p[i][j][k][pol] + pMp[mz].Cd[pol] * Ez_p[i][j][k]
								+ pMp[mz].Ce[pol] * Ez_pp[i][j][k];
							Jz_p[i][j][k][pol] = tempz;
						}
					}
				}
	}
}

void TFDTD_Generic_PBC::H_field_update()
{
	//n_Thread = omp_get_max_threads(); // omp_get_max_thread();
	//omp_set_num_threads(n_Thread);
	int i, j, k;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i = 0; i < txsize; i++)
			for (j = 0; j < tysize - 1; j++)
				for (k = 0; k < tzsize; k++)
				{
					Hx[i][j][k] = Hx[i][j][k]
						- Db / H_dy[j] * (Ez[i][j + 1][k] - Ez[i][j][k])
						+ Db / H_dz[k] * (Ey[i][j][k + 1] - Ey[i][j][k]);
				}

#pragma omp for schedule(dynamic)
		for (i = 0; i < txsize - 1; i++)
			for (j = 0; j < tysize; j++)
				for (k = 0; k < tzsize; k++)
				{
					Hy[i][j][k] = Hy[i][j][k]
						- Db / H_dz[k] * (Ex[i][j][k + 1] - Ex[i][j][k])
						+ Db / H_dx[i] * (Ez[i + 1][j][k] - Ez[i][j][k]);
				}

#pragma omp for schedule(dynamic)
		for (i = 0; i < txsize - 1; i++)
			for (j = 0; j < tysize - 1; j++)
				for (k = 0; k < tzsize; k++)
				{
					Hz[i][j][k] = Hz[i][j][k]
						- Db / H_dx[i] * (Ey[i + 1][j][k] - Ey[i][j][k])
						+ Db / H_dy[j] * (Ex[i][j + 1][k] - Ex[i][j][k]);
				}
	}
}

void TFDTD_Generic_PBC::Modeling_YeeCell()
{
	// int locX = log10(txsize)
	// int locY = log10(tysize)
	// int aaLoc1 = log10(aa1m1)
	// int aaLoc2 = log10(aa1m2)
	// for (int ii = 0 ; ii <= locX ; ii++)
	//  for (int jj = 0 ; jj <= locY ; hh++)
	//		if (ii == 0 || ii == locX-1)
	//				printf("*")
	// 


	for (int num_ind = 0; num_ind < num_material; num_ind++)
	{
		for (i = 0; i < txsize; i++)
			for (j = 0; j < tysize; j++)
				for (k = 0; k < tzsize; k++)
				{
					if (Media[i][j][k] == num_ind)
					{
						Media_Ex[i][j][k] = Media[i][j][k];
						Media_Ex[i][j + 1][k] = Media[i][j][k];
						Media_Ex[i][j][k + 1] = Media[i][j][k];
						Media_Ex[i][j + 1][k + 1] = Media[i][j][k];

						Media_Ey[i][j][k] = Media[i][j][k];
						Media_Ey[i + 1][j][k] = Media[i][j][k];
						Media_Ey[i][j][k + 1] = Media[i][j][k];
						Media_Ey[i + 1][j][k + 1] = Media[i][j][k];

						Media_Ez[i][j][k] = Media[i][j][k];
						Media_Ez[i + 1][j][k] = Media[i][j][k];
						Media_Ez[i][j + 1][k] = Media[i][j][k];
						Media_Ez[i + 1][j + 1][k] = Media[i][j][k];

						Media_Hx[i][j][k] = Media[i][j][k];
						Media_Hx[i + 1][j][k] = Media[i][j][k];

						Media_Hy[i][j][k] = Media[i][j][k];
						Media_Hy[i][j + 1][k] = Media[i][j][k];

						Media_Hz[i][j][k] = Media[i][j][k];
						Media_Hz[i][j][k + 1] = Media[i][j][k];
					}
				}
	}

	int aa1m1, aa1m2;
	FILE* pf1;
	pf1 = fopen("CHIRALPBC.dat", "r");
	fscanf(pf1, "%d,%d\n", &aa1m1, &aa1m2);
	fclose(pf1);
	printf("@Media\naa1m1=%d\naa1m2=%d\n", aa1m1, aa1m2);

	i = txsize - 1;
	int j_offset = 0;
	for (j = 0; j < aa1m1; j++)
	{
		j_offset = j + aa1m2;
		for (k = 0; k < tzsize; k++)
		{
			Media_Ey[0][j][k] = Media_Ey[i][j_offset][k];
			Media_Ez[0][j][k] = Media_Ez[i][j_offset][k];
		}
	}
	for (j = aa1m1; j < tysize; j++)
	{
		j_offset = j - aa1m1;
		for (k = 0; k < tzsize; k++)
		{
			Media_Ey[0][j][k] = Media_Ey[i][j_offset][k];
			Media_Ez[0][j][k] = Media_Ez[i][j_offset][k];
		}
	}

	j = tysize - 1;
	for (i = 0; i < txsize; i++)
	{
		for (k = 0; k < tzsize; k++)
		{
			Media_Ex[i][0][k] = Media_Ex[i][j][k];
			Media_Ez[i][0][k] = Media_Ez[i][j][k];
		}
	}

	//	printf("%d\n", Media);
	// memfree(&Media, txsize, tysize);

}