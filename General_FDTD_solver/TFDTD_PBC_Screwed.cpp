#include "TFDTD_PBC_Screwed.h"

TFDTD_PBC_Screwed::TFDTD_PBC_Screwed(struct_Field_carrier* pFc)
{
	eps_0 = 1e-9 / 36 / m_pi;
	mu_0 = 4 * m_pi * 1e-7;
	c_0 = 1 / sqrt(eps_0 * mu_0);
	this->pFc = pFc;

	this->Ex_p = pFc->Ex_p_address;
	this->Ey_p = pFc->Ey_p_address;
	this->Ez_p = pFc->Ez_p_address;
	
	this->Ex_pp = pFc->Ex_pp_address;
	this->Ey_pp = pFc->Ey_pp_address;
	this->Ez_pp = pFc->Ez_pp_address;

	this->Jx = pFc->Jx_address;
	this->Jy = pFc->Jy_address;
	this->Jz = pFc->Jz_address;

	this->Jx_p = pFc->Jx_p_address;
	this->Jy_p = pFc->Jy_p_address;
	this->Jz_p = pFc->Jz_p_address;
	this->dt = pFc->dt;

	Read_Material_Data();
	set_PBC_Screwed_Parameters();
	PML_Initial_Setting();
}

TFDTD_PBC_Screwed::TFDTD_PBC_Screwed()
{
	eps_0 = 1e-9 / 36 / m_pi;
	mu_0 = 4 * m_pi * 1e-7;
	c_0 = 1 / sqrt(eps_0 * mu_0);
}

void TFDTD_PBC_Screwed::set_PBC_Screwed_Parameters()
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

void TFDTD_PBC_Screwed::FDTD_PBC_Chiral_BC1_omp()
{
	//n_Thread = omp_get_max_threads(); // omp_get_max_thread();
	//omp_set_num_threads(n_Thread);
	int i = 0, j = 0, k = 0, pol = 0;
	double temp1 = 0, temp2 = 0, temp3 = 0, temp4 = 0, temp5 = 0, temp6 = 0;
	int mx, my, mz;
#pragma omp parallel private(i,k,pol, temp1, temp2, temp5, temp6, mx, my, mz)
	{
		j = 0;
#pragma omp for schedule(dynamic)
		for (i = 0; i < txsize - 1; i++)
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
					+ Cb[mx] / E_dy[j] * (Hz[i][j][k] - Hz[i][tysize - 2][k])
					- Cb[mx] / E_dz[k] * (Hy[i][j][k] - Hy[i][j][k - 1]) - temp1 - temp2;

			}

		j = 0;
#pragma omp for schedule(dynamic)
		for (i = 1; i < txsize - 1; i++)
			for (k = 1; k < tzsize; k++)
			{
				mz = Media_Ez[i][j][k];
				temp5 = 0.; temp6 = 0.;
				if (pMp[mz].TypeOfMaterials == 'D')
				{
					for (pol = 0; pol < total_pol; pol++)
					{
						temp5 += Ca_temp[mz][pol] * Jz[i][j][k][pol];;
						temp6 += Cb_temp[mz][pol] * Jz_p[i][j][k][pol];;
					}
				}
				Ez[i][j][k] = C1[mz] * Ez_p[i][j][k] - C2[mz] * Ez_pp[i][j][k]
					+ Cb[mz] / E_dx[i] * (Hy[i][j][k] - Hy[i - 1][j][k])
					- Cb[mz] / E_dy[j] * (Hx[i][j][k] - Hx[i][tysize - 2][k]) - temp5 - temp6;

			}

		j = tysize - 1;
#pragma omp for schedule(dynamic)
		for (i = 0; i < txsize - 1; i++)
			for (k = 0; k < tzsize; k++)
			{
				Ex[i][j][k] = Ex[i][0][k];
			}


		j = tysize - 1;
#pragma omp for schedule(dynamic)
		for (i = 0; i < txsize; i++)
			for (k = 0; k < tzsize; k++)
			{
				Ez[i][j][k] = Ez[i][0][k];
			}
	}
}

void TFDTD_PBC_Screwed::FDTD_PBC_Chiral_BC2_omp()
{
	int j_offset = 0;
	int media = 0;
	int comp = 0;
	int i = 0, j = 0, k = 0, pol = 0;
	double temp1 = 0, temp2 = 0, temp3 = 0, temp4 = 0, temp5 = 0, temp6 = 0;
	int mx, my, mz;
	
// 	cout << 1 / E_dx[0] << " " << 1 / E_dy[0] << endl;
// 	cout << "Hdiff = " << (Hy[0][10][10] - Hy[txsize - 2][10][10]) << endl;

#pragma omp parallel private(i,j,k, pol, j_offset, temp1, temp2, temp3, temp4, temp5, temp6, mx, my, mz)
	{
		i = 0;
#pragma omp for schedule(dynamic)
		for (j = 0; j < aa1m1; j++)
		{
			j_offset = j + aa1m2;
			for (k = 1; k < tzsize; k++)
			{
				my = Media_Ey[i][j][k];
				temp3 = 0.; temp4 = 0.;
				if (pMp[my].TypeOfMaterials == 'D')
				{
					for (int pol = 0; pol < total_pol; pol++)
					{
						temp3 += Ca_temp[my][pol] * Jy[i][j][k][pol];
						temp4 += Cb_temp[my][pol] * Jy_p[i][j][k][pol];
					}
				}
				Ey[i][j][k] = C1[my] * Ey_p[i][j][k] - C2[my] * Ey_pp[i][j][k]
					+ Cb[my] / E_dz[k] * (Hx[i][j][k] - Hx[i][j][k - 1])
					- Cb[my] / E_dx[i] * (Hz[i][j][k] - Hz[txsize - 2][j_offset][k]) - temp3 - temp4;
			}
		}

		i = 0;
#pragma omp for schedule(dynamic)
		for (j = aa1m1; j < tysize - 1; j++)
		{
			j_offset = j - aa1m1;
			for (k = 1; k < tzsize; k++)
			{
				my = Media_Ey[i][j][k];
				temp3 = 0.; temp4 = 0.;
				if (pMp[my].TypeOfMaterials == 'D')
				{
					for (int pol = 0; pol < total_pol; pol++)
					{
						temp3 += Ca_temp[my][pol] * Jy[i][j][k][pol];
						temp4 += Cb_temp[my][pol] * Jy_p[i][j][k][pol];
					}
				}
				Ey[i][j][k] = C1[my] * Ey_p[i][j][k] - C2[my] * Ey_pp[i][j][k]
					+ Cb[my] / E_dz[k] * (Hx[i][j][k] - Hx[i][j][k - 1])
					- Cb[my] / E_dx[i] * (Hz[i][j][k] - Hz[txsize - 2][j_offset][k]) - temp3 - temp4;

			}
		}

		i = 0;
#pragma omp for schedule(dynamic)
		for (j = 1; j < aa1m1; j++)
		{
			j_offset = j + aa1m2;
			for (k = 1; k < tzsize; k++)
			{
				mz = Media_Ez[i][j][k];
				temp5 = 0.; temp6 = 0.;
				if (pMp[mz].TypeOfMaterials == 'D')
				{
					for (int pol = 0; pol < total_pol; pol++)
					{
						temp5 += Ca_temp[mz][pol] * Jz[i][j][k][pol];;
						temp6 += Cb_temp[mz][pol] * Jz_p[i][j][k][pol];;
					}
				}
				Ez[i][j][k] = C1[mz] * Ez_p[i][j][k] - C2[mz] * Ez_pp[i][j][k]
					+ Cb[mz] / E_dx[i] * (Hy[i][j][k] - Hy[txsize - 2][j_offset][k])
					- Cb[mz] / E_dy[j] * (Hx[i][j][k] - Hx[i][j - 1][k]) - temp5 - temp6;
			}
		}

		i = 0;
		int j_offset2 = 0;
#pragma omp for schedule(dynamic)
		for (j = aa1m1; j < tysize - 1; j++)
		{
			j_offset = j - aa1m1;
			if (j == 0)
				j_offset2 = tysize - 2;
			else
				j_offset2 = j - 1;
			for (k = 1; k < tzsize; k++)
			{
				mz = Media_Ez[i][j][k];
				temp5 = 0.; temp6 = 0.;
				if (pMp[mz].TypeOfMaterials == 'D')
				{
					for (int pol = 0; pol < total_pol; pol++)
					{
						temp5 += Ca_temp[mz][pol] * Jz[i][j][k][pol];
						temp6 += Cb_temp[mz][pol] * Jz_p[i][j][k][pol];
					}
				}
				Ez[i][j][k] = C1[mz] * Ez_p[i][j][k] - C2[mz] * Ez_pp[i][j][k]
					+ Cb[mz] / E_dx[i] * (Hy[i][j][k] - Hy[txsize - 2][j_offset][k])
					- Cb[mz] / E_dy[j] * (Hx[i][j][k] - Hx[i][j_offset2][k]) - temp5 - temp6;
			}
		}

		i = 0;
		j = 0;
#pragma omp for schedule(dynamic)
		for (k = 1; k < tzsize; k++)
		{
			mz = Media_Ez[i][j][k];
			temp5 = 0.; temp6 = 0.;
			if (pMp[mz].TypeOfMaterials == 'D')
			{
				for (int pol = 0; pol < total_pol; pol++)
				{
					temp5 += Ca_temp[mz][pol] * Jz[i][j][k][pol];;
					temp6 += Cb_temp[mz][pol] * Jz_p[i][j][k][pol];;
				}
			}
			Ez[i][j][k] = C1[mz] * Ez_p[i][j][k] - C2[mz] * Ez_pp[i][j][k]
				+ Cb[mz] / E_dx[i] * (Hy[i][j][k] - Hy[txsize - 2][aa1m2][k])
				- Cb[mz] / E_dy[j] * (Hx[i][j][k] - Hx[txsize - 1][aa1m2 - 1][k]) - temp5 - temp6;

		}


		//i==txsize PBC
		i = txsize - 1;
#pragma omp for schedule(dynamic)
		for (j = 0; j < aa1m2; j++) {
			j_offset = j + aa1m1;
			for (k = 0; k < tzsize; k++) {
				//            if( *(*(*(*Medi_Ey+i)+j)+k) - *(*(*(*Medi_Ey+0)+j_offset)+k) != 0)
				//            ShowMessage("Not Matching");
				Ey[i][j][k] = Ey[0][j_offset][k];
			}
		}

		i = txsize - 1;
#pragma omp for schedule(dynamic)
		for (j = aa1m2; j < tysize - 1; j++) {
			j_offset = j - aa1m2;
			for (k = 0; k < tzsize; k++) {
				//            if( *(*(*(*Medi_Ey+i)+j)+k) - *(*(*(*Medi_Ey+0)+j_offset)+k) != 0)
				//            ShowMessage("Not Matching");
				Ey[i][j][k] = Ey[0][j_offset][k];
			}
		}
		i = txsize - 1;
#pragma omp for schedule(dynamic)
		for (j = 0; j < aa1m2; j++) {
			j_offset = j + aa1m1;
			for (k = 1; k < tzsize; k++) {
				//            if( *(*(*(*Medi_Ez+i)+j)+k) - *(*(*(*Medi_Ez+0)+j_offset)+k) != 0)
				//            ShowMessage("Not Matching");
				Ez[i][j][k] = Ez[0][j_offset][k];
			}
		}

		i = txsize - 1;
#pragma omp for schedule(dynamic)
		for (j = aa1m2; j < tysize - 1; j++) {
			j_offset = j - aa1m2;
			for (k = 1; k < tzsize; k++) {
				//            if( *(*(*(*Medi_Ez+i)+j)+k) - *(*(*(*Medi_Ez+0)+j_offset)+k) != 0)
				//            ShowMessage("Not Matching");
				Ez[i][j][k] = Ez[0][j_offset][k];
			}
		}
	}
}

//void TFDTD_Generic_PBC::FDTD_PBC_Chiral_BC1_omp()
//{
//	//n_Thread = omp_get_max_threads(); // omp_get_max_thread();
//	//omp_set_num_threads(n_Thread);
//	int i = 0, j = 0, k = 0, pol = 0;
//	double temp1 = 0, temp2 = 0, temp3 = 0, temp4 = 0, temp5 = 0, temp6 = 0;
//#pragma omp parallel private(i,k,pol, temp1, temp2, temp5, temp6)
//	{
//		j = 0;
//#pragma omp for schedule(dynamic)
//		for (i = 0; i < txsize - 1; i++)
//			for (k = 1; k < tzsize; k++)
//			{
//				mx = Media_Ex[i][j][k];
//				temp1 = 0.; temp2 = 0.;
//				if (pMp[mx].TypeOfMaterials == 'D')
//				{
//					for (pol = 0; pol < total_pol; pol++)
//					{
//						temp1 += Ca_temp[mx][pol] * Jx[i][j][k][pol];
//						temp2 += Cb_temp[mx][pol] * Jx_p[i][j][k][pol];
//					}
//				}
//				Ex[i][j][k] = C1[mx] * Ex_p[i][j][k] - C2[mx] * Ex_pp[i][j][k]
//					+ Cb[mx] / E_dy[j] * (Hz[i][j][k] - Hz[i][tysize - 2][k])
//					- Cb[mx] / E_dz[k] * (Hy[i][j][k] - Hy[i][j][k - 1]) - temp1 - temp2;
//
//			}
//
//		j = 0;
//#pragma omp for schedule(dynamic)
//		for (i = 1; i < txsize - 1; i++)
//			for (k = 1; k < tzsize; k++)
//			{
//				mz = Media_Ez[i][j][k];
//				temp5 = 0.; temp6 = 0.;
//				if (pMp[mz].TypeOfMaterials == 'D')
//				{
//					for (pol = 0; pol < total_pol; pol++)
//					{
//						temp5 += Ca_temp[mz][pol] * Jz[i][j][k][pol];;
//						temp6 += Cb_temp[mz][pol] * Jz_p[i][j][k][pol];;
//					}
//				}
//				Ez[i][j][k] = C1[mz] * Ez_p[i][j][k] - C2[mz] * Ez_pp[i][j][k]
//					+ Cb[mz] / E_dx[i] * (Hy[i][j][k] - Hy[i - 1][j][k])
//					- Cb[mz] / E_dy[j] * (Hx[i][j][k] - Hx[i][tysize - 2][k]) - temp5 - temp6;
//
//			}
//
//		j = tysize - 1;
//#pragma omp for schedule(dynamic)
//		for (i = 0; i < txsize - 1; i++)
//			for (k = 0; k < tzsize; k++)
//			{
//				Ex[i][j][k] = Ex[i][0][k];
//			}
//
//
//		j = tysize - 1;
//#pragma omp for schedule(dynamic)
//		for (i = 0; i < txsize; i++)
//			for (k = 0; k < tzsize; k++)
//			{
//				Ez[i][j][k] = Ez[i][0][k];
//			}
//	}
//}
//
//void TFDTD_Generic_PBC::FDTD_PBC_Chiral_BC2_omp()
//{
//	int j_offset = 0;
//	int media = 0;
//	int comp = 0;
//	int i = 0, j = 0, k = 0, pol = 0;
//	double temp1 = 0, temp2 = 0, temp3 = 0, temp4 = 0, temp5 = 0, temp6 = 0;
//#pragma omp parallel private(i,j,k, pol, j_offset, temp1, temp2, temp3, temp4, temp5, temp6)
//	{
//		i = 0;
//#pragma omp for schedule(dynamic)
//		for (j = 0; j < aa1m1; j++)
//		{
//			j_offset = j + aa1m2;
//			for (k = 1; k < tzsize; k++)
//			{
//				my = Media_Ey[i][j][k];
//				temp3 = 0.; temp4 = 0.;
//				if (pMp[my].TypeOfMaterials == 'D')
//				{
//					for (int pol = 0; pol < total_pol; pol++)
//					{
//						temp3 += Ca_temp[my][pol] * Jy[i][j][k][pol];
//						temp4 += Cb_temp[my][pol] * Jy_p[i][j][k][pol];
//					}
//				}
//				Ey[i][j][k] = C1[my] * Ey_p[i][j][k] - C2[my] * Ey_pp[i][j][k]
//					+ Cb[my] / E_dz[k] * (Hx[i][j][k] - Hx[i][j][k - 1])
//					- Cb[my] / E_dx[i] * (Hz[i][j][k] - Hz[txsize - 2][j_offset][k]) - temp3 - temp4;
//			}
//		}
//
//		i = 0;
//#pragma omp for schedule(dynamic)
//		for (j = aa1m1; j < tysize - 1; j++)
//		{
//			j_offset = j - aa1m1;
//			for (k = 1; k < tzsize; k++)
//			{
//				my = Media_Ey[i][j][k];
//				temp3 = 0.; temp4 = 0.;
//				if (pMp[my].TypeOfMaterials == 'D')
//				{
//					for (int pol = 0; pol < total_pol; pol++)
//					{
//						temp3 += Ca_temp[my][pol] * Jy[i][j][k][pol];
//						temp4 += Cb_temp[my][pol] * Jy_p[i][j][k][pol];
//					}
//				}
//				Ey[i][j][k] = C1[my] * Ey_p[i][j][k] - C2[my] * Ey_pp[i][j][k]
//					+ Cb[my] / E_dz[k] * (Hx[i][j][k] - Hx[i][j][k - 1])
//					- Cb[my] / E_dx[i] * (Hz[i][j][k] - Hz[txsize - 2][j_offset][k]) - temp3 - temp4;
//
//			}
//		}
//
//		i = 0;
//#pragma omp for schedule(dynamic)
//		for (j = 1; j < aa1m1; j++)
//		{
//			j_offset = j + aa1m2;
//			for (k = 1; k < tzsize; k++)
//			{
//				mz = Media_Ez[i][j][k];
//				temp5 = 0.; temp6 = 0.;
//				if (pMp[mz].TypeOfMaterials == 'D')
//				{
//					for (int pol = 0; pol < total_pol; pol++)
//					{
//						temp5 += Ca_temp[mz][pol] * Jz[i][j][k][pol];;
//						temp6 += Cb_temp[mz][pol] * Jz_p[i][j][k][pol];;
//					}
//				}
//				Ez[i][j][k] = C1[mz] * Ez_p[i][j][k] - C2[mz] * Ez_pp[i][j][k]
//					+ Cb[mz] / E_dx[i] * (Hy[i][j][k] - Hy[txsize - 2][j_offset][k])
//					- Cb[mz] / E_dy[j] * (Hx[i][j][k] - Hx[i][j - 1][k]) - temp5 - temp6;
//
//			}
//		}
//
//		i = 0;
//		int j_offset2 = 0;
//#pragma omp for schedule(dynamic)
//		for (j = aa1m1; j < tysize - 1; j++)
//		{
//			j_offset = j - aa1m1;
//			if (j == 0)
//				j_offset2 = tysize - 2;
//			else
//				j_offset2 = j - 1;
//			for (k = 1; k < tzsize; k++)
//			{
//				mz = Media_Ez[i][j][k];
//				temp5 = 0.; temp6 = 0.;
//				if (pMp[mz].TypeOfMaterials == 'D')
//				{
//					for (int pol = 0; pol < total_pol; pol++)
//					{
//						temp5 += Ca_temp[mz][pol] * Jz[i][j][k][pol];
//						temp6 += Cb_temp[mz][pol] * Jz_p[i][j][k][pol];
//					}
//				}
//				Ez[i][j][k] = C1[mz] * Ez_p[i][j][k] - C2[mz] * Ez_pp[i][j][k]
//					+ Cb[mz] / E_dx[i] * (Hy[i][j][k] - Hy[txsize - 2][j_offset][k])
//					- Cb[mz] / E_dy[j] * (Hx[i][j][k] - Hx[i][j_offset2][k]) - temp5 - temp6;
//			}
//		}
//
//		i = 0;
//		j = 0;
//#pragma omp for schedule(dynamic)
//		for (k = 1; k < tzsize; k++)
//		{
//			mz = Media_Ez[i][j][k];
//			temp5 = 0.; temp6 = 0.;
//			if (pMp[mz].TypeOfMaterials == 'D')
//			{
//				for (int pol = 0; pol < total_pol; pol++)
//				{
//					temp5 += Ca_temp[mz][pol] * Jz[i][j][k][pol];;
//					temp6 += Cb_temp[mz][pol] * Jz_p[i][j][k][pol];;
//				}
//			}
//			Ez[i][j][k] = C1[mz] * Ez_p[i][j][k] - C2[mz] * Ez_pp[i][j][k]
//				+ Cb[mz] / E_dx[i] * (Hy[i][j][k] - Hy[txsize - 2][aa1m2][k])
//				- Cb[mz] / E_dy[j] * (Hx[i][j][k] - Hx[txsize - 1][aa1m2 - 1][k]) - temp5 - temp6;
//
//		}
//
//
//		//i==txsize PBC
//		i = txsize - 1;
//#pragma omp for schedule(dynamic)
//		for (j = 0; j < aa1m2; j++) {
//			j_offset = j + aa1m1;
//			for (k = 0; k < tzsize; k++) {
//				//            if( *(*(*(*Medi_Ey+i)+j)+k) - *(*(*(*Medi_Ey+0)+j_offset)+k) != 0)
//				//            ShowMessage("Not Matching");
//				Ey[i][j][k] = Ey[0][j_offset][k];
//			}
//		}
//
//		i = txsize - 1;
//#pragma omp for schedule(dynamic)
//		for (j = aa1m2; j < tysize - 1; j++) {
//			j_offset = j - aa1m2;
//			for (k = 0; k < tzsize; k++) {
//				//            if( *(*(*(*Medi_Ey+i)+j)+k) - *(*(*(*Medi_Ey+0)+j_offset)+k) != 0)
//				//            ShowMessage("Not Matching");
//				Ey[i][j][k] = Ey[0][j_offset][k];
//			}
//		}
//		i = txsize - 1;
//#pragma omp for schedule(dynamic)
//		for (j = 0; j < aa1m2; j++) {
//			j_offset = j + aa1m1;
//			for (k = 1; k < tzsize; k++) {
//				//            if( *(*(*(*Medi_Ez+i)+j)+k) - *(*(*(*Medi_Ez+0)+j_offset)+k) != 0)
//				//            ShowMessage("Not Matching");
//				Ez[i][j][k] = Ez[0][j_offset][k];
//			}
//		}
//
//		i = txsize - 1;
//#pragma omp for schedule(dynamic)
//		for (j = aa1m2; j < tysize - 1; j++) {
//			j_offset = j - aa1m2;
//			for (k = 1; k < tzsize; k++) {
//				//            if( *(*(*(*Medi_Ez+i)+j)+k) - *(*(*(*Medi_Ez+0)+j_offset)+k) != 0)
//				//            ShowMessage("Not Matching");
//				Ez[i][j][k] = Ez[0][j_offset][k];
//			}
//		}
//	}
//}

//void TFDTD_Generic_PBC::FDTD_PBC_Chiral_BC1()
//{
//	//n_Thread = omp_get_max_threads(); // omp_get_max_thread();
//	//omp_set_num_threads(n_Thread);
//	j = 0;
//	for (i = 0; i < txsize - 1; i++) {
//		for (k = 0; k < tzsize; k++) {
//			Dx[i][j][k] = Dx[i][j][k] + Cb_dy * (Hz[i][j][k] - Hz[i][tysize - 2][k])
//				- g_cb_dz[k] * (Hy[i][j][k] - Hy[i][j][k - 1]);
//		}
//	}
//
//	for (i = 1; i < txsize - 1; i++) {
//		for (k = 0; k < tzsize; k++) {
//			Dz[i][j][k] = Dz[i][j][k] + Cb_dx*(Hy[i][j][k] - Hy[i - 1][j][k])
//				- Cb_dy*(Hx[i][j][k] - Hx[i][tysize - 2][k]);
//		}
//	}
//
//	j = tysize - 1;
//	for (i = 0; i < txsize - 1; i++) {
//		for (k = 0; k < tzsize; k++) {
//			//        if( *(*(*(*Medi_Ex+i)+j)+k) - *(*(*(*Medi_Ex+i)+0)+k) != 0)
//			//            ShowMessage("Not Matching");
//			Dx[i][j][k] = Dx[i][0][k];
//		}
//	}
//
//	for (i = 0; i < txsize; i++) {
//		for (k = 0; k < tzsize; k++) {
//			//            if( *(*(*(*Medi_Ez+i)+j)+k) - *(*(*(*Medi_Ez+i)+0)+k) != 0)
//			//                ShowMessage("Not Matching");
//			Dz[i][j][k] = Dz[i][0][k];
//		}
//	}
//}
//
//void TFDTD_Generic_PBC::FDTD_PBC_Chiral_BC2()
//{
//	int j_offset = 0;
//	int media = 0;
//	int comp = 0;
//	//n_Thread = omp_get_max_threads(); // omp_get_max_thread();
//	//omp_set_num_threads(n_Thread);
//
//	i = 0;
//	for (j = 0; j < aa1m1; j++) 
//	{
//		j_offset = j + aa1m2;
//		for (k = 0; k < tzsize; k++) 
//			Dy[i][j][k] = Dy[i][j][k] + g_cb_dz[k] * (Hx[i][j][k] - Hx[i][j][k - 1])
//			- Cb_dx * (Hz[i][j][k] - Hz[txsize - 2][j_offset][k]);
//	}
//
//
//		for (j = aa1m1; j < tysize - 1; j++) {
//			j_offset = j - aa1m1;
//			for (k = 0; k < tzsize; k++) {
//				Dy[i][j][k] = Dy[i][j][k] + g_cb_dz[k] * (Hx[i][j][k] - Hx[i][j][k - 1])
//					- Cb_dx * (Hz[i][j][k] - Hz[txsize - 1][j_offset][k]);
//			}
//		}
//
//		for (j = 1; j < aa1m1; j++) {
//			j_offset = j + aa1m2;
//			for (k = 0; k < tzsize; k++) {
//				Dz[i][j][k] = Dz[i][j][k] + Cb_dx*(Hy[i][j][k] - Hy[txsize - 2][j_offset][k])
//					- Cb_dy*(Hx[i][j][k] - Hx[i][j - 1][k]);
//			}
//		}
//		int j_offset2 = 0;
//		for (j = aa1m1; j < tysize - 1; j++) {
//			j_offset = j - aa1m1;
//			if (j == 0)
//				j_offset2 = tysize - 2;
//			else
//				j_offset2 = j - 1;
//			for (k = 0; k < tzsize; k++) 
//			{
//				Dz[i][j][k] = Dz[i][j][k] + Cb_dx*(Hy[i][j][k] - Hy[txsize - 2][j_offset][k])
//										  - Cb_dy*(Hx[i][j][k] - Hx[i][j_offset2][k]);
//			}
//		}
//	
//
//
//	i = 0;
//	j = 0;
//	
//		for (k = 0; k < tzsize; k++)
//		{
//			Dz[i][j][k] = Dz[i][j][k] 
//				+ Cb_dx*(Hy[i][j][k] - Hy[txsize - 2][aa1m2][k])
//				- Cb_dy*(Hx[i][j][k] - Hx[txsize - 1][aa1m2 - 1][k]);
//		}
//	
//	//i==txsize PBC
//	i = txsize - 1;
//	
//		for (j = 0; j < aa1m2; j++) {
//			j_offset = j + aa1m1;
//			for (k = 0; k < tzsize; k++) {
//				//            if( *(*(*(*Medi_Ey+i)+j)+k) - *(*(*(*Medi_Ey+0)+j_offset)+k) != 0)
//				//            ShowMessage("Not Matching");
//				Dy[i][j][k] = Dy[0][j_offset][k];
//			}
//		}
//
//		for (j = aa1m2; j < tysize - 1; j++) {
//			j_offset = j - aa1m2;
//			for (k = 0; k < tzsize; k++) {
//				//            if( *(*(*(*Medi_Ey+i)+j)+k) - *(*(*(*Medi_Ey+0)+j_offset)+k) != 0)
//				//            ShowMessage("Not Matching");
//				Dy[i][j][k] = Dy[0][j_offset][k];
//			}
//		}
//
//			
//		for (j = 0; j < aa1m2; j++) {
//			j_offset = j + aa1m1;
//			for (k = 1; k < tzsize; k++) {
//				//            if( *(*(*(*Medi_Ez+i)+j)+k) - *(*(*(*Medi_Ez+0)+j_offset)+k) != 0)
//				//            ShowMessage("Not Matching");
//				Dz[i][j][k] = Dz[0][j_offset][k];
//			}
//		}
//	
//		for (j = aa1m2; j < tysize - 1; j++) {
//			j_offset = j - aa1m2;
//			for (k = 1; k < tzsize; k++) {
//				//            if( *(*(*(*Medi_Ez+i)+j)+k) - *(*(*(*Medi_Ez+0)+j_offset)+k) != 0)
//				//            ShowMessage("Not Matching");
//				Dz[i][j][k] = Dz[0][j_offset][k];
//			}
//		}
//	
//}

