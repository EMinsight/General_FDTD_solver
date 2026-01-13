//---------------------------------------------------------------------------
#pragma hdrstop
#include "TFDTD_GenericAlgorithm.h"

//---------------------------------------------------------------------------
#pragma package(smart_init)

TFDTD_GenericAlgorithm::TFDTD_GenericAlgorithm()
{
	eps_0 = 1e-9 / 36 / m_pi;
	mu_0 = 4 * m_pi * 1e-7;
	c_0 = 1 / sqrt(eps_0*mu_0);
    isFDTD=1;
}

TFDTD_GenericAlgorithm::~TFDTD_GenericAlgorithm(){}

struct_Field_carrier TFDTD_GenericAlgorithm::get_Field_address(struct_Field_carrier* Fc)
{
	Fc->Ex_p_address = Ex_p;
	Fc->Ey_p_address = Ey_p;
	Fc->Ez_p_address = Ez_p;

	Fc->Ex_pp_address = Ex_pp;
	Fc->Ey_pp_address = Ey_pp;
	Fc->Ez_pp_address = Ez_pp;

	Fc->Jx_address = Jx;
	Fc->Jy_address = Jy;
	Fc->Jz_address = Jz;

	Fc->Jx_p_address = Jx_p;
	Fc->Jy_p_address = Jy_p;
	Fc->Jz_p_address = Jz_p;

	Fc->media = Media;
	Fc->Media_Ex = Media_Ex;
	Fc->Media_Ey = Media_Ey;
	Fc->Media_Ez = Media_Ez;

	Fc->dt = dt;
	Fc->dx = dx;
	Fc->dy = dy;
	Fc->dz = dz;
//	Fc->total_pol = total_pol;
	return *Fc;
}

void TFDTD_GenericAlgorithm::Modeling_YeeCell()
{
	for (int num_ind = 0; num_ind < num_material; num_ind++)
	{
		for (i = 0; i<txsize; i++)
        for (j = 0; j<tysize; j++)
        for (k = 0; k<tzsize; k++)
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
	//	printf("%d\n", Media);
	// memfree(&Media, txsize, tysize);
}

void TFDTD_GenericAlgorithm::Modeling_Meta_Cell()
{
    for (int num_ind = 0; num_ind < num_material; num_ind++)
	{
		for (i = 0; i<txsize; i++)
        for (j = 0; j<tysize; j++)
        for (k = 0; k<tzsize; k++)
        {
            if (Media[i][j][k] ==  num_material-1)
            {
                Media_Ex[i][j][k] = Media[i][j][k];
                Media_Ex[i][j + 1][k] = Media[i][j][k];

                Media_Ey[i][j][k] = Media[i][j][k];
                Media_Ey[i + 1][j][k] = Media[i][j][k];


                Media_Ez[i][j][k] = Media[i][j][k+1];
                Media_Ez[i + 1][j][k] = Media[i][j][k+1];
                Media_Ez[i][j + 1][k] = Media[i][j][k+1];
                Media_Ez[i + 1][j + 1][k] = Media[i][j][k+1];


                Media_Hx[i][j][k] = Media[i][j][k];
                Media_Hx[i + 1][j][k] = Media[i][j][k];

                Media_Hy[i][j][k] = Media[i][j][k];
                Media_Hy[i][j + 1][k] = Media[i][j][k];

                Media_Hz[i][j][k] = Media[i][j][k];
                Media_Hz[i][j][k + 1] = Media[i][j][k];
            }
            else if (Media[i][j][k] == num_ind)
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
}
void TFDTD_GenericAlgorithm::Field_Initiate(struct_Field_carrier* sFc)
{
	this->Ex = sFc->Ex_address;
	this->Ey = sFc->Ey_address;
	this->Ez = sFc->Ez_address;
	
	this->Jx = sFc->Jx_address;
	this->Jy = sFc->Jy_address;
	this->Jz = sFc->Jz_address;
	
	this->Hx = sFc->Hx_address;
	this->Hy = sFc->Hy_address;
	this->Hz = sFc->Hz_address;

	Jx_p = Tptr4<double>(txsize, tysize + 1, tzsize + 1, total_pol);
	Jy_p = Tptr4<double>(txsize + 1, tysize, tzsize + 1, total_pol);
	Jz_p = Tptr4<double>(txsize + 1, tysize + 1, tzsize, total_pol);

	Ex_p = Tptr3<double>(txsize, tysize + 1, tzsize + 1);
	Ey_p = Tptr3<double>(txsize + 1, tysize, tzsize + 1);
	Ez_p = Tptr3<double>(txsize + 1, tysize + 1, tzsize);

	Ex_pp = Tptr3<double>(txsize, tysize + 1, tzsize + 1);
	Ey_pp = Tptr3<double>(txsize + 1, tysize, tzsize + 1);
	Ez_pp = Tptr3<double>(txsize + 1, tysize + 1, tzsize);

	sFc->Jx_p_address = Jx_p;
	sFc->Jy_p_address = Jy_p;
	sFc->Jz_p_address = Jz_p;
	sFc->Ex_p_address = Ex_p;
	sFc->Ey_p_address = Ey_p;
	sFc->Ez_p_address = Ez_p;
	sFc->Ex_pp_address = Ex_pp;
	sFc->Ey_pp_address = Ey_pp;
	sFc->Ez_pp_address = Ez_pp;

	Sfc = sFc;
	num_thread = Sfc->num_thread;
	num_thread_pef_Fnc = (num_thread-3)/3.0;
	printf("num_thread : %d per update Fnc\n", num_thread_pef_Fnc);
	omp_set_max_active_levels(2);
}

void TFDTD_GenericAlgorithm::Field_Initiate(double ***Ex, double ***Ey, double ***Ez, double ****Jx, double ****Jy, double ****Jz, double ***Hx, double ***Hy, double ***Hz)
{
	this->Ex = Ex;
	this->Ey = Ey;
	this->Ez = Ez;

	this->Jx = Jx;
    this->Jy = Jy;
    this->Jz = Jz;

	this->Hx = Hx;
	this->Hy = Hy;
	this->Hz = Hz;

	Jx_p = Tptr4<double>(txsize, tysize + 1, tzsize + 1, total_pol);
	Jy_p = Tptr4<double>(txsize + 1, tysize, tzsize + 1, total_pol);
	Jz_p = Tptr4<double>(txsize + 1, tysize + 1, tzsize, total_pol);

	Ex_p  = Tptr3<double>(txsize, tysize + 1, tzsize + 1);
	Ey_p  = Tptr3<double>(txsize + 1, tysize, tzsize + 1);
	Ez_p  = Tptr3<double>(txsize + 1, tysize + 1, tzsize);

	Ex_pp = Tptr3<double>(txsize, tysize + 1, tzsize + 1);
	Ey_pp = Tptr3<double>(txsize + 1, tysize, tzsize + 1);
	Ez_pp = Tptr3<double>(txsize + 1, tysize + 1, tzsize);
}


void TFDTD_GenericAlgorithm::Ex_previous()
{
	int i=0, j=0;
#pragma omp parallel for collapse(2) schedule(static) private(i,j) //num_threads(num_thread_pef_Fnc)
		for (i = 0; i < txsize; i++)
		for (j = 0; j < tysize + 1; j++)
		{
			memcpy(Ex_pp[i][j], Ex_p[i][j], sizeof(double) * tzsize);
			memcpy(Ex_p[i][j], Ex[i][j], sizeof(double) * tzsize);
		}
}

void TFDTD_GenericAlgorithm::Ey_previous()
{
	int i=0, j=0;
#pragma omp parallel for collapse(2) schedule(static) private(i,j) //num_threads(num_thread_pef_Fnc)
		for (i = 0; i < txsize + 1; i++)
		for (j = 0; j < tysize; j++)
		{
			memcpy(Ey_pp[i][j], Ey_p[i][j], sizeof(double) * tzsize);
			memcpy(Ey_p[i][j], Ey[i][j], sizeof(double) * tzsize);
		}
}
void TFDTD_GenericAlgorithm::Ez_previous()
{
	int i=0, j=0;
#pragma omp parallel for collapse(2) schedule(static) private(i,j) //num_threads(num_thread_pef_Fnc)
		for (i = 0; i < txsize + 1; i++)
		for (j = 0; j < tysize + 1; j++)
		{
			memcpy(Ez_pp[i][j], Ez_p[i][j], sizeof(double) * tzsize);
			memcpy(Ez_p[i][j], Ez[i][j], sizeof(double) * tzsize);
		}

}

void TFDTD_GenericAlgorithm::Ex_field_update() 
{
	int i=0, j=0, k=0, pol = 0;
	int mx=0;
	double temp1 = 0, temp2 = 0;
#pragma omp parallel for collapse(3) schedule(static) private(i,j,k, pol, mx, temp1, temp2) //num_threads(num_thread_pef_Fnc)
		for (i = 0; i < txsize; i++)
		for (j = 1; j < tysize; j++)
		for (k = 1; k < tzsize; k++)
		{
			temp1 = 0;
			temp2 = 0;
			mx = Media_Ex[i][j][k];
			// if (total_pol > 0)
			// {
			// 	temp1 = cblas_ddot(total_pol, Ca_temp[mx], 1, Jx[i][j][k], 1);
			// 	temp2 = cblas_ddot(total_pol, Cb_temp[mx], 1, Jx_p[i][j][k], 1);
			// }

				temp1 = 0.;
				temp2 = 0.;
				for (pol = 0; pol < total_pol; pol++)
				{
					temp1 += Ca_temp[mx][pol] * Jx[i][j][k][pol];
					temp2 += Cb_temp[mx][pol] * Jx_p[i][j][k][pol];
				}

			Ex[i][j][k] = C1[mx] * Ex_p[i][j][k] - C2[mx] * Ex_pp[i][j][k] 				
				+ Cb[mx] / E_dy[j] * (Hz[i][j][k] - Hz[i][j - 1][k])				
				- Cb[mx] / E_dz[k] * (Hy[i][j][k] - Hy[i][j][k - 1]) - temp1 - temp2;
		}
}

void TFDTD_GenericAlgorithm::Ey_field_update() 
{
	int i=0, j=0, k=0, pol = 0;
	int my=0;
	double temp3 = 0, temp4 = 0;
#pragma omp parallel for collapse(3)  schedule(static) private(i,j,k, pol, my, temp3, temp4) //num_threads(num_thread_pef_Fnc)
		for (i = 1; i < txsize; i++)
		for (j = 0; j < tysize; j++)
		for (k = 1; k < tzsize; k++)
		{
			temp3 = 0;
			temp4 = 0;
			my = Media_Ey[i][j][k];
			// if (total_pol > 0)
			// {
			// 	temp3 = cblas_ddot(total_pol, Ca_temp[my], 1, Jy[i][j][k], 1);
			// 	temp4 = cblas_ddot(total_pol, Cb_temp[my], 1, Jy_p[i][j][k], 1);
			// }
			temp3 = 0.;
			temp4 = 0.;
			for (pol = 0; pol < total_pol; pol++)
			{
				temp3 += Ca_temp[my][pol] * Jy[i][j][k][pol];
				temp4 += Cb_temp[my][pol] * Jy_p[i][j][k][pol];
			}
			Ey[i][j][k] = C1[my] * Ey_p[i][j][k] - C2[my] * Ey_pp[i][j][k] 
						+ Cb[my] / E_dz[k] * (Hx[i][j][k] - Hx[i][j][k - 1])
						- Cb[my] / E_dx[i] * (Hz[i][j][k] - Hz[i - 1][j][k]) - temp3 - temp4;
		}
}
void TFDTD_GenericAlgorithm::Ez_field_update()
{
	int i=0, j=0, k=0, pol = 0;
	int mz=0;
	double temp5 = 0, temp6 = 0;
#pragma omp parallel for  collapse(3) schedule(static) private(i, j, k, pol, mz, temp5, temp6) //num_threads(num_thread_pef_Fnc)
		for (i = 1; i < txsize; i++)
		for (j = 1; j < tysize; j++)
		for (k = 0; k < tzsize; k++)
		{
			temp5 = 0;
			temp6 = 0;
			mz = Media_Ey[i][j][k];
			// if (total_pol > 0)
			// {
			// 	temp5 = cblas_ddot(total_pol, Ca_temp[mz], 1, Jz[i][j][k], 1);
			// 	temp6 = cblas_ddot(total_pol, Cb_temp[mz], 1, Jz_p[i][j][k], 1);
			// }
			temp5 = 0.;
			temp6 = 0.;
			mz = Media_Ez[i][j][k];
			for (pol = 0; pol < total_pol; pol++)
			{
				temp5 += Ca_temp[mz][pol] * Jz[i][j][k][pol];;
				temp6 += Cb_temp[mz][pol] * Jz_p[i][j][k][pol];;
			}

			Ez[i][j][k] = C1[mz] * Ez_p[i][j][k] - C2[mz] * Ez_pp[i][j][k]
				+ Cb[mz] / E_dx[i] * (Hy[i][j][k] - Hy[i - 1][j][k])
				- Cb[mz] / E_dy[j] * (Hx[i][j][k] - Hx[i][j - 1][k]) - temp5 - temp6;

		}
}

void TFDTD_GenericAlgorithm::E_field_update()
{
	//Ex_previous();
	//Ey_previous();
	//Ez_previous();
	//Ex_field_update();
	//Ey_field_update();
	//Ez_field_update();

	//cout << omp_get_max_active_levels() << endl;

#pragma omp parallel num_threads(3) 
	{
		#pragma omp single nowait
			#pragma omp task 
			Ex_previous();
			// Ex_field_update();
			#pragma omp task 
			Ey_previous();
			// Ey_field_update();
			#pragma omp task 
			Ez_previous();
			// Ez_field_update();
			#pragma omp taskwait
	}

#pragma omp parallel num_threads(3) 
	{
		#pragma omp single nowait
		{
			#pragma omp task 
			Ex_field_update();
			#pragma omp task 
			Ey_field_update();
			#pragma omp task 
			Ez_field_update();
		}
		#pragma omp taskwait
	}

}

// void TFDTD_GenericAlgorithm::E_field_update()
// {
// 	double temp_save;
// 	int i, j, k, pol;
// 	int mx, my, mz;
// 
// 	double temp1 = 0, temp2 = 0, temp3 = 0, temp4 = 0, temp5 = 0, temp6 = 0;
// #pragma omp parallel private(i,j,k, pol, temp1, temp2, temp3, temp4, temp5, temp6)
// 	{
// 
// #pragma omp for schedule(dynamic)
// 		for (i = 0; i < txsize; i++)
// 			for (j = 0; j < tysize + 1; j++)
// 			{
// 				memcpy(Ex_pp[i][j], Ex_p[i][j], sizeof(double) * tzsize);
// 				memcpy(Ex_p[i][j], Ex[i][j], sizeof(double) * tzsize);
// 			}
// 
// #pragma omp for schedule(dynamic) 
// 		for (i = 0; i < txsize + 1; i++)
// 			for (j = 0; j < tysize; j++)
// 			{
// 				memcpy(Ey_pp[i][j], Ey_p[i][j], sizeof(double) * tzsize);
// 				memcpy(Ey_p[i][j], Ey[i][j], sizeof(double) * tzsize);
// 			}
// 
// #pragma omp for schedule(dynamic)
// 		for (i = 0; i < txsize + 1; i++)
// 			for (j = 0; j < tysize + 1; j++)
// 			{
// 				memcpy(Ez_pp[i][j], Ez_p[i][j], sizeof(double) * tzsize);
// 				memcpy(Ez_p[i][j], Ez[i][j], sizeof(double) * tzsize);
// 			}
// 
// #pragma omp for schedule(dynamic)
// 		for (i = 0; i < txsize; i++)
// 			for (j = 1; j < tysize; j++)
// 				for (k = 1; k < tzsize; k++)
// 				{
// 					mx = Media_Ex[i][j][k];
// 					temp1 = 0.;
// 					temp2 = 0.;
// 					for (int pol = 0; pol < total_pol; pol++)
// 					{
// 						temp1 += Ca_temp[mx][pol] * Jx[i][j][k][pol];
// 						temp2 += Cb_temp[mx][pol] * Jx_p[i][j][k][pol];
// 					}
// 					Ex[i][j][k] = C1[mx] * Ex_p[i][j][k] - C2[mx] * Ex_pp[i][j][k] 
// 						+ Cb[mx] / E_dy[j] * (Hz[i][j][k] - Hz[i][j - 1][k])
// 						- Cb[mx] / E_dz[k] * (Hy[i][j][k] - Hy[i][j][k - 1]) - temp1 - temp2;
// 				}
// 
// 
// #pragma omp for schedule(dynamic)
// 		for (i = 1; i < txsize; i++)
// 			for (j = 0; j < tysize; j++)
// 				for (k = 1; k < tzsize; k++)
// 				{
// 					my = Media_Ey[i][j][k];
// 
// 					temp3 = 0.;
// 					temp4 = 0.;
// 
// 					for (int pol = 0; pol < total_pol; pol++)
// 					{
// 						temp3 += Ca_temp[my][pol] * Jy[i][j][k][pol];
// 						temp4 += Cb_temp[my][pol] * Jy_p[i][j][k][pol];
// 					}
// 					Ey[i][j][k] = C1[my] * Ey_p[i][j][k] - C2[my] * Ey_pp[i][j][k] 
// 						+ Cb[my] / E_dz[k] * (Hx[i][j][k] - Hx[i][j][k - 1])
// 						- Cb[my] / E_dx[i] * (Hz[i][j][k] - Hz[i - 1][j][k]) - temp3 - temp4;
// 					//Ey[i][j][k] = Ey[i][j][k] + Cb_dz[0] / kappa_e_zprof[k] * (Hx[i][j][k] - Hx[i][j][k - 1]) - Cb_dx[0] / kappa_e_xprof[i] * (Hz[i][j][k] - Hz[i - 1][j][k]);
// 
// 				}
// #pragma omp for schedule(dynamic)
// 		for (i = 1; i < txsize; i++)
// 			for (j = 1; j < tysize; j++)
// 				for (k = 0; k < tzsize; k++)
// 				{
// 					mz = Media_Ez[i][j][k];
// 
// 					temp5 = 0.;
// 					temp6 = 0.;
// 
// 					for (int pol = 0; pol < total_pol; pol++)
// 					{
// 						temp5 += Ca_temp[mz][pol] * Jz[i][j][k][pol];;
// 						temp6 += Cb_temp[mz][pol] * Jz_p[i][j][k][pol];;
// 					}
// 
// 					Ez[i][j][k] = C1[mz] * Ez_p[i][j][k] - C2[mz] * Ez_pp[i][j][k] 
// 						+ Cb[mz] / E_dx[i] * (Hy[i][j][k] - Hy[i - 1][j][k])
// 						- Cb[mz] / E_dy[j]  * (Hx[i][j][k] - Hx[i][j - 1][k]) - temp5 - temp6;
// 					//Ez[i][j][k] = Ez[i][j][k] + Cb_dx[0] / kappa_e_xprof[i] * (Hy[i][j][k] - Hy[i - 1][j][k]) - Cb_dy[0] / kappa_e_yprof[j] * (Hx[i][j][k] - Hx[i][j - 1][k]);
// 				}
// 	}
// }

void TFDTD_GenericAlgorithm::Jx_field_update()
{
	int i = 0, j = 0, k = 0, pol = 0;
	int mx=0;
	double tempx=0;
#pragma omp parallel for collapse(3) schedule(static) private(i,j,k,mx, pol,tempx) //num_threads(num_thread_pef_Fnc)
		for (i = 0; i < txsize; i++)
		for (j = 1; j < tysize; j++)
		for (k = 1; k < tzsize; k++)
		{
			mx = Media_Ex[i][j][k];
			if (pMp[mx].TypeOfMaterials == 'D')
			{
				for (int pol = 0; pol < total_pol; pol++)
				{

					tempx = Jx[i][j][k][pol];
					Jx[i][j][k][pol] = pMp[mx].Ca[pol] * Jx[i][j][k][pol] + pMp[mx].Cb[pol] * Jx_p[i][j][k][pol] + pMp[mx].Cc[pol] * Ex[i][j][k] + pMp[mx].Cd[pol] * Ex_p[i][j][k] + pMp[mx].Ce[pol] * Ex_pp[i][j][k];
					Jx_p[i][j][k][pol] = tempx;
				}
			}
		}
}
void TFDTD_GenericAlgorithm::Jy_field_update()
{
	int i = 0, j = 0, k = 0, pol = 0;
	int my=0;
	double tempy=0;
#pragma omp parallel for collapse(3) schedule(static) private(i,j,k,my, pol,tempy) //num_threads(num_thread_pef_Fnc)
		for (i = 1; i < txsize; i++)
		for (j = 0; j < tysize; j++)
		for (k = 1; k < tzsize; k++)
		{
			my = Media_Ey[i][j][k];
			if (pMp[my].TypeOfMaterials == 'D')
			{
				for (int pol = 0; pol < total_pol; pol++)
				{
					tempy = Jy[i][j][k][pol];
					Jy[i][j][k][pol] = pMp[my].Ca[pol] * Jy[i][j][k][pol] + pMp[my].Cb[pol] * Jy_p[i][j][k][pol] + pMp[my].Cc[pol] * Ey[i][j][k] + pMp[my].Cd[pol] * Ey_p[i][j][k] + pMp[my].Ce[pol] * Ey_pp[i][j][k];
					Jy_p[i][j][k][pol] = tempy;
				}
			}
		}
}
void TFDTD_GenericAlgorithm::Jz_field_update()
{
	int i = 0, j = 0, k = 0, pol = 0;
	int mz=0;
	double tempz=0;
#pragma omp parallel for collapse(3) schedule(static) private(i,j,k,mz, pol,tempz) //num_threads(num_thread_pef_Fnc)
		for (i = 1; i < txsize; i++)
		for (j = 1; j < tysize; j++)
		for (k = 0; k < tzsize; k++)
		{
			mz = Media_Ez[i][j][k];
			if (pMp[mz].TypeOfMaterials == 'D')
			{
				for (int pol = 0; pol < total_pol; pol++)
				{
					tempz = Jz[i][j][k][pol];
					Jz[i][j][k][pol] = pMp[mz].Ca[pol] * Jz[i][j][k][pol] + pMp[mz].Cb[pol] * Jz_p[i][j][k][pol] + pMp[mz].Cc[pol] * Ez[i][j][k] + pMp[mz].Cd[pol] * Ez_p[i][j][k] + pMp[mz].Ce[pol] * Ez_pp[i][j][k];
					Jz_p[i][j][k][pol] = tempz;
				}
			}
		}
}

void TFDTD_GenericAlgorithm::J_field_update()
{
	// Jx_field_update();
	// Jy_field_update();
	// Jz_field_update();
#pragma omp parallel num_threads(3) 
	{
		#pragma omp single nowait
		{
			#pragma omp task 
			Jx_field_update();
			#pragma omp task 
			Jy_field_update();
			#pragma omp task 
			Jz_field_update();
		}
		#pragma omp taskwait
	}
}
// void TFDTD_GenericAlgorithm::J_field_update()
// {
// 	int i = 0, j = 0, k = 0, pol = 0;
// 	int mx, my, mz;
// 	double tempx, tempy, tempz;
// #pragma omp parallel private(i,j,k,mx, my, mz, pol,tempx,tempy,tempz)
// 	{
// #pragma omp for schedule(static)
// 		for (i = 0; i < txsize; i++) 
// 		for (j = 1; j < tysize; j++) 
// 		for (k = 1; k < tzsize; k++) 
// 		for (int pol = 0; pol < total_pol; pol++)
// 		{
// 			mx = Media_Ex[i][j][k];
// 			if (pMp[mx].TypeOfMaterials == 'D')
// 			{
// 				tempx = Jx[i][j][k][pol];
// 				Jx[i][j][k][pol] = pMp[mx].Ca[pol] * Jx[i][j][k][pol] + pMp[mx].Cb[pol] * Jx_p[i][j][k][pol] + pMp[mx].Cc[pol] * Ex[i][j][k] + pMp[mx].Cd[pol] * Ex_p[i][j][k] + pMp[mx].Ce[pol] * Ex_pp[i][j][k];
// 				Jx_p[i][j][k][pol] = tempx;
// 			}
// 		}
// 
// 
// #pragma omp for schedule(static)
// 		for (i = 1; i < txsize; i++) 
// 		for (j = 0; j < tysize; j++) 
// 		for (k = 1; k < tzsize; k++) 
// 		for (int pol = 0; pol < total_pol; pol++)
// 		{
// 			my = Media_Ey[i][j][k];
// 			if (pMp[my].TypeOfMaterials == 'D')
// 			{
// 				tempy = Jy[i][j][k][pol];
// 				Jy[i][j][k][pol] = pMp[my].Ca[pol] * Jy[i][j][k][pol] + pMp[my].Cb[pol] * Jy_p[i][j][k][pol] + pMp[my].Cc[pol] * Ey[i][j][k] + pMp[my].Cd[pol] * Ey_p[i][j][k] + pMp[my].Ce[pol] * Ey_pp[i][j][k];
// 				Jy_p[i][j][k][pol] = tempy;
// 			}
// 		}
// 
// 
// #pragma omp for schedule(static)
// 		for (i = 1; i < txsize; i++) 
// 		for (j = 1; j < tysize; j++) 
// 		for (k = 0; k < tzsize; k++) 
// 		for (int pol = 0; pol < total_pol; pol++)
// 		{
// 			mz = Media_Ez[i][j][k];
// 			if (pMp[mz].TypeOfMaterials == 'D')
// 			{
// 				tempz = Jz[i][j][k][pol];
// 				Jz[i][j][k][pol] = pMp[mz].Ca[pol] * Jz[i][j][k][pol] + pMp[mz].Cb[pol] * Jz_p[i][j][k][pol] + pMp[mz].Cc[pol] * Ez[i][j][k] + pMp[mz].Cd[pol] * Ez_p[i][j][k] + pMp[mz].Ce[pol] * Ez_pp[i][j][k];
// 				Jz_p[i][j][k][pol] = tempz;
// 			}
// 		}
// 	}
// }