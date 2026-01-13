//---------------------------------------------------------------------------
#pragma hdrstop
#include "TFDTD_Boundary_TFSF.h"
#include <iostream>
using namespace std;
//---------------------------------------------------------------------------
__fastcall TFDTD_Boundary_TFSF::TFDTD_Boundary_TFSF(double ***Ex, double ***Ey, double ***Ez, double ***Hx, double ***Hy, double ***Hz)
{
	m_pi = acos(-1);
	eps_0 = 1e-9 / 36 / m_pi;
	mu_0 = 4 * m_pi * 1e-7;
	c_0 = 1.0 / sqrt(eps_0*mu_0);

	this->Ex = Ex; this->Hx = Hx;
	this->Ey = Ey; this->Hy = Hy;
	this->Ez = Ez; this->Hz = Hz;
}

TFDTD_Boundary_TFSF::TFDTD_Boundary_TFSF(struct_Field_carrier* pFc)
{
	this->pFc = pFc;
	m_pi = acos(-1);
	eps_0 = 1e-9 / 36 / m_pi;
	mu_0 = 4 * m_pi * 1e-7;
	c_0 = 1.0 / sqrt(eps_0 * mu_0);

	this->Ex = pFc->Ex_address; this->Hx = pFc->Hx_address;
	this->Ey = pFc->Ey_address; this->Hy = pFc->Hy_address;
	this->Ez = pFc->Ez_address; this->Hz = pFc->Hz_address;

	this->dx = pFc->dx;
	this->dy = pFc->dy;
	this->dz = pFc->dz;
	this->dt = pFc->dt;
	this->f0 = pFc->f0;
	this->txsize = pFc->txsize;
	this->tysize = pFc->tysize;
	this->tzsize = pFc->tzsize;
}

__fastcall TFDTD_Boundary_TFSF::~TFDTD_Boundary_TFSF()
{
	memfree(&Ex_xyN_correction, inc_plane_x);
	memfree(&Ex_xyP_correction, inc_plane_x);
	memfree(&Ex_xzN_correction, inc_plane_x);
	memfree(&Ex_xzP_correction, inc_plane_x);

	memfree(&Ey_xyN_correction, inc_plane_x);
	memfree(&Ey_xyP_correction, inc_plane_x);
	memfree(&Ey_yzN_correction, inc_plane_y);
	memfree(&Ey_yzP_correction, inc_plane_y);

	memfree(&Ez_xzN_correction, inc_plane_x);
	memfree(&Ez_xzP_correction, inc_plane_x);
	memfree(&Ez_yzN_correction, inc_plane_y);
	memfree(&Ez_yzP_correction, inc_plane_y);

	memfree(&Hx_xyN_correction, inc_plane_x);
	memfree(&Hx_xyP_correction, inc_plane_x);
	memfree(&Hx_xzN_correction, inc_plane_x);
	memfree(&Hx_xzP_correction, inc_plane_x);

	memfree(&Hy_xyN_correction, inc_plane_x);
	memfree(&Hy_xyP_correction, inc_plane_x);
	memfree(&Hy_yzN_correction, inc_plane_y);
	memfree(&Hy_yzP_correction, inc_plane_y);

	memfree(&Hz_xzN_correction, inc_plane_x);
	memfree(&Hz_xzP_correction, inc_plane_x);
	memfree(&Hz_yzN_correction, inc_plane_y);
	memfree(&Hz_yzP_correction, inc_plane_y);
}

//void TFDTD_Boundary_TFSF::insert_parameters(int txsize, int tysize, int tzsize, double dx, double dy, double dz, double dt, double f0, int Nstep)
//{
//	this->dx = dx;
//	this->dy = dy;
//	this->dz = dz;
//	this->dt = dt;
//	this->f0 = f0;
//	this->txsize = txsize;
//	this->tysize = tysize;
//	this->tzsize = tzsize;
//	Cedx = dt / eps_0 / dx;
//	Cedy = dt / eps_0 / dy;
//	Cedz = dt / eps_0 / dz;
//
//	Chdx = dt / mu_0 / dx;
//	Chdy = dt / mu_0 / dy;
//	Chdz = dt / mu_0 / dz;
//}

void TFDTD_Boundary_TFSF::insert_inc_angle(double input_Theta, double input_phi, double polar_psi)
{
	this->Input_Degree(input_phi, input_Theta, polar_psi);
}

void TFDTD_Boundary_TFSF::insert_TFSF_Box(int xx00, int xx11, int yy00, int yy11, int zz00, int zz11)
{
	Cedx1 = dt / eps_0 / pFc->E_dx[xx00];
	Cedy1 = dt / eps_0 / pFc->E_dy[yy00];
	Cedz1 = dt / eps_0 / pFc->E_dz[zz00];
	Chdx1 = dt / mu_0 / pFc->H_dx[xx00];
	Chdy1 = dt / mu_0 / pFc->H_dy[yy00];
	Chdz1 = dt / mu_0 / pFc->H_dz[zz00];

	Cedx2 = dt / eps_0 / pFc->E_dx[xx11];
	Cedy2 = dt / eps_0 / pFc->E_dy[yy11];
	Cedz2 = dt / eps_0 / pFc->E_dz[zz11];
	Chdx2 = dt / mu_0 / pFc->H_dx[xx11];
	Chdy2 = dt / mu_0 / pFc->H_dy[yy11];
	Chdz2 = dt / mu_0 / pFc->H_dz[zz11];

	this->Input_Boundary(xx00, xx11, yy00, yy11, zz00, zz11);
	this->Initial_condition();
	this->Define_Initial_Position(); // <- dds define
	this->Calcuate_Initial_Field();
}

// void TFDTD_Boundary_TFSF::Set_Pulse(double Pulse)
// {
// 	this->Pulse = Pulse;
// 	Electric_Field_update();
// 	Magnetic_Field_update();
// 	Calcuate_inc_Field();
// }

void TFDTD_Boundary_TFSF::DApply()
{
	if (B_Tx_ini != 0)
		Correction_X_Boundary_Exn();
	if (B_Tx_ter != txsize)
		Correction_X_Boundary_Exp();
	if (B_Ty_ini != 0)
		Correction_Y_Boundary_Eyn();
	if (B_Ty_ter != tysize)
		Correction_Y_Boundary_Eyp();
	if (B_Tz_ini != 0)
		Correction_Z_Boundary_Ezn();
	if (B_Tz_ter != tzsize)
		Correction_Z_Boundary_Ezp();
}

void TFDTD_Boundary_TFSF::HApply()
{
	if (B_Tx_ini != 0)
		Correction_X_Boundary_Hxn();
	if (B_Tx_ter != txsize)
		Correction_X_Boundary_Hxp();
	if (B_Ty_ini != 0)
		Correction_Y_Boundary_Hyn();
	if (B_Ty_ter != tysize)
		Correction_Y_Boundary_Hyp();
	if (B_Tz_ini != 0)
		Correction_Z_Boundary_Hzn();
	if (B_Tz_ter != tzsize)
		Correction_Z_Boundary_Hzp();
}

void TFDTD_Boundary_TFSF::Input_Degree(double input_phi, double input_theta, double polar_psi)
{
	this->input_phi = input_phi;
	this->input_theta = input_theta;
	Phi = input_phi * 2.0 * m_pi / 360;     //Radian
	Theta = input_theta * 2.0 * m_pi / 360; //Radian
	this->Polar_psi = polar_psi * 2.0 * m_pi / 360;

	kk_inc = abs((sin(Theta)*cos(Phi) + sin(Theta)*sin(Phi) + cos(Theta)));
}

void TFDTD_Boundary_TFSF::Get_Polarization(double polar_psi)
{
	this->Polar_psi = polar_psi * 2.0 * m_pi / 360;
}

void TFDTD_Boundary_TFSF::Input_Boundary(int xi, int xt, int yi, int yt, int zi, int zt)
{
	B_Tx_ini = xi;
	B_Tx_ter = xt;

	B_Ty_ini = yi;
	B_Ty_ter = yt;

	B_Tz_ini = zi;
	B_Tz_ter = zt;

	disd = 10;
}

void TFDTD_Boundary_TFSF::Define_dds(int xpos, int ypos, int zpos)
{
	dds = (double)sqrt(sin(Theta) * cos(Phi) * pFc->E_dx[xpos] * sin(Theta) * cos(Phi) * pFc->E_dx[xpos]
		+ sin(Theta) * sin(Phi) * pFc->E_dy[ypos] * sin(Theta) * sin(Phi) * pFc->E_dy[ypos]
		+ cos(Theta) * pFc->E_dz[zpos] * cos(Theta) * pFc->E_dz[zpos]);
}

void TFDTD_Boundary_TFSF::Define_Initial_Position()
{
	if (input_theta >= 0 && input_theta <  90)
	{
		if (input_phi >= 0 && input_phi <= 90)
		{
			DD_center = B_Tx_ini*sin(Theta)*cos(Phi) + B_Ty_ini*sin(Theta)*sin(Phi) + B_Tz_ini*cos(Theta);
			DD_center_H = (B_Tx_ini - 0.5 * sin(Theta) * cos(Phi)) * sin(Theta) * cos(Phi) + (B_Ty_ini - 0.5 * sin(Theta) * sin(Phi)) * sin(Theta) * sin(Phi) + (B_Tz_ini - 0.5 * cos(Theta)) * cos(Theta);
			Define_dds(B_Tx_ini, B_Ty_ini, B_Tz_ini);
		}
		else if (input_phi > 90 && input_phi <= 180)
		{
			DD_center = B_Tx_ter*sin(Theta)*cos(Phi) + B_Ty_ini*sin(Theta)*sin(Phi) + B_Tz_ini*cos(Theta);
			DD_center_H = (B_Tx_ter - 0.5 * sin(Theta) * cos(Phi)) * sin(Theta) * cos(Phi) + (B_Ty_ini - 0.5 * sin(Theta) * sin(Phi)) * sin(Theta) * sin(Phi) + (B_Tz_ini - 0.5 * cos(Theta)) * cos(Theta);
			Define_dds(B_Tx_ter, B_Ty_ini, B_Tz_ini);
		}
		else if (input_phi > 180 && input_phi <= 270)
		{
			DD_center = B_Tx_ter*sin(Theta)*cos(Phi) + B_Ty_ter*sin(Theta)*sin(Phi) + B_Tz_ini*cos(Theta);
			DD_center_H = (B_Tx_ter - 0.5 * sin(Theta) * cos(Phi)) * sin(Theta) * cos(Phi) + (B_Ty_ter - 0.5 * sin(Theta) * sin(Phi)) * sin(Theta) * sin(Phi) + (B_Tz_ini - 0.5 * cos(Theta)) * cos(Theta);
			Define_dds(B_Tx_ter, B_Ty_ter, B_Tz_ini);
		}
		else if (input_phi > 270 && input_phi < 359)
		{
			DD_center = B_Tx_ini*sin(Theta)*cos(Phi) + B_Ty_ter*sin(Theta)*sin(Phi) + B_Tz_ini*cos(Theta);
			DD_center_H = (B_Tx_ini - 0.5 * sin(Theta) * cos(Phi)) * sin(Theta) * cos(Phi) + (B_Ty_ter - 0.5 * sin(Theta) * sin(Phi)) * sin(Theta) * sin(Phi) + (B_Tz_ini - 0.5 * cos(Theta)) * cos(Theta);
			Define_dds(B_Tx_ini, B_Ty_ter, B_Tz_ini);
		}
	}
	else if (input_theta >= 90 && input_theta <= 180)
	{
		if (input_phi >= 0 && input_phi <= 90)
		{
			DD_center = B_Tx_ini*sin(Theta)*cos(Phi) + B_Ty_ini*sin(Theta)*sin(Phi) + B_Tz_ter*cos(Theta);
			DD_center_H = (B_Tx_ini - 0.5 * sin(Theta) * cos(Phi)) * sin(Theta) * cos(Phi) + (B_Ty_ini - 0.5 * sin(Theta) * sin(Phi)) * sin(Theta) * sin(Phi) + (B_Tz_ter - 0.5 * cos(Theta)) * cos(Theta);
			Define_dds(B_Tx_ini, B_Ty_ini, B_Tz_ter);
		}
		else if (input_phi > 90 && input_phi <= 180)
		{
			DD_center = B_Tx_ter*sin(Theta)*cos(Phi) + B_Ty_ini*sin(Theta)*sin(Phi) + B_Tz_ter*cos(Theta);
			DD_center_H = (B_Tx_ter - 0.5 * sin(Theta) * cos(Phi)) * sin(Theta) * cos(Phi) + (B_Ty_ini - 0.5 * sin(Theta) * sin(Phi)) * sin(Theta) * sin(Phi) + (B_Tz_ter - 0.5 * cos(Theta)) * cos(Theta);
			Define_dds(B_Tx_ter, B_Ty_ini, B_Tz_ter);
		}
		else if (input_phi > 180 && input_phi <= 270)
		{
			DD_center = B_Tx_ter*sin(Theta)*cos(Phi) + B_Ty_ter*sin(Theta)*sin(Phi) + B_Tz_ter*cos(Theta);
			DD_center_H = (B_Tx_ter - 0.5 * sin(Theta) * cos(Phi)) * sin(Theta) * cos(Phi) + (B_Ty_ter - 0.5 * sin(Theta) * sin(Phi)) * sin(Theta) * sin(Phi) + (B_Tz_ter - 0.5 * cos(Theta)) * cos(Theta);
			Define_dds(B_Tx_ter, B_Ty_ter, B_Tz_ter);
		}
		else if (input_phi > 270 && input_phi < 359)
		{
			DD_center = B_Tx_ini*sin(Theta)*cos(Phi) + B_Ty_ter*sin(Theta)*sin(Phi) + B_Tz_ter*cos(Theta);
			DD_center_H = (B_Tx_ini - 0.5 * sin(Theta) * cos(Phi)) * sin(Theta) * cos(Phi) + (B_Ty_ter - 0.5 * sin(Theta) * sin(Phi)) * sin(Theta) * sin(Phi) + (B_Tz_ter - 0.5 * cos(Theta)) * cos(Theta);
			Define_dds(B_Tx_ini, B_Ty_ter, B_Tz_ter);
		}
	}
}

void TFDTD_Boundary_TFSF::Initial_condition()
{
	this->Nststep = (int)(5*sqrt((double)(txsize*txsize + tysize*tysize + tzsize*tzsize + 0.5)));

	E_inc_oblic = Tptr1<double>(Nststep);
	H_inc_oblic = Tptr1<double>(Nststep);
	E_inc_oblic_mur = Tptr1<double>(2);

	inc_plane_x = B_Tx_ter - B_Tx_ini + 1;
	inc_plane_y = B_Ty_ter - B_Ty_ini + 1;
	inc_plane_z = B_Tz_ter - B_Tz_ini + 1;

	Ex_xyN_correction = Tptr2<double>(inc_plane_x, inc_plane_y);
	Ex_xyP_correction = Tptr2<double>(inc_plane_x, inc_plane_y);
	Ex_xzN_correction = Tptr2<double>(inc_plane_x, inc_plane_z);
	Ex_xzP_correction = Tptr2<double>(inc_plane_x, inc_plane_z);

	Ey_xyN_correction = Tptr2<double>(inc_plane_x, inc_plane_y);
	Ey_xyP_correction = Tptr2<double>(inc_plane_x, inc_plane_y);
	Ey_yzN_correction = Tptr2<double>(inc_plane_y, inc_plane_z);
	Ey_yzP_correction = Tptr2<double>(inc_plane_y, inc_plane_z);

	Ez_xzN_correction = Tptr2<double>(inc_plane_x, inc_plane_z);
	Ez_xzP_correction = Tptr2<double>(inc_plane_x, inc_plane_z);
	Ez_yzN_correction = Tptr2<double>(inc_plane_y, inc_plane_z);
	Ez_yzP_correction = Tptr2<double>(inc_plane_y, inc_plane_z);

	Hx_xyN_correction = Tptr2<double>(inc_plane_x, inc_plane_y);
	Hx_xyP_correction = Tptr2<double>(inc_plane_x, inc_plane_y);
	Hx_xzN_correction = Tptr2<double>(inc_plane_x, inc_plane_z);
	Hx_xzP_correction = Tptr2<double>(inc_plane_x, inc_plane_z);

	Hy_xyN_correction = Tptr2<double>(inc_plane_x, inc_plane_y);
	Hy_xyP_correction = Tptr2<double>(inc_plane_x, inc_plane_y);
	Hy_yzN_correction = Tptr2<double>(inc_plane_y, inc_plane_z);
	Hy_yzP_correction = Tptr2<double>(inc_plane_y, inc_plane_z);

	Hz_xzN_correction = Tptr2<double>(inc_plane_x, inc_plane_z);
	Hz_xzP_correction = Tptr2<double>(inc_plane_x, inc_plane_z);
	Hz_yzN_correction = Tptr2<double>(inc_plane_y, inc_plane_z);
	Hz_yzP_correction = Tptr2<double>(inc_plane_y, inc_plane_z);
}

void TFDTD_Boundary_TFSF::Calcuate_Initial_Field()
{
	///////////////////////////////////////////////////////////////////////////
	//	double dds = (double)sqrt(((double)((dx*dx + dy*dy + dz*dz))));//abs(sin(Theta)*cos(Phi)*dx + sin(Phi)*cos(Theta)*dy + cos(Theta)*dz);
	double S_CFS = c_0*dt / dds; //*sqrt(((double)((1 / (dx*dx) + 1 / (dy*dy) + 1 / (dz*dz)) / 3)));
								 //double S_CFS = c_0*dt/dds; //*sqrt(((double)((1 / (dx*dx) + 1 / (dy*dy) + 1 / (dz*dz)) / 3)));
	double N_lam = (c_0 / f0) / dds;
	double psi_dir = (double)sqrt((double)(sin(Theta)*cos(Phi)*sin(Theta)*cos(Phi))) +
		(double)sqrt((double)(sin(Theta)*sin(Phi)*sin(Theta)*sin(Phi))) +
		(double)sqrt((double)(cos(Theta)*cos(Theta)));
	v_nomal = m_pi*c_0 / (N_lam          *asin(1.0 / S_CFS 			*sin(m_pi*S_CFS / N_lam)));
	v_oblic = m_pi*c_0 / (N_lam / psi_dir*asin(1.0 / S_CFS * psi_dir*sin(m_pi*S_CFS / N_lam)));
	v_ratio = v_nomal / v_oblic;

	// Mur
	double rho = c_0*dt / dds;
	abs_Gamma_media = (1.0 - rho) / (1.0 + rho);

	//	double dds = dz;
	//	//double dds = 1.0/sqrt((double)(1.0/dx * 1.0/dx + 1.0/dy * 1.0/dy + 1.0/dz * 1.0/dz));
	//	//double dds   = 1.0/sqrt((double)(1.0/dx * 1.0/dx + 1.0/dy * 1.0/dy + 1.0/dz * 1.0/dz);
	//	double N_lam = (c_0 / f0) / dds;
	//	double S_CFS = c_0*dt / dds;
	//	double psi_dir = abs(sin(Theta)*cos(Phi) + sin(Theta)*sin(Phi) + cos(Theta));
	//	v_nomal = m_pi*c_0 / (N_lam*asin(1.0 / S_CFS*sin(m_pi*S_CFS / N_lam)));
	//	v_oblic = m_pi*c_0 / (N_lam*psi_dir*asin(1.0 / S_CFS / psi_dir*sin(m_pi*S_CFS / N_lam)));
	//	// n_delay = 1*dx/(v_oblic*dt);
	//	v_ratio = v_nomal / v_oblic;

	FILE *pfInc;
	pfInc = fopen("Incident_field_plane.dat", "w");
	fclose(pfInc);
	pf_debug1 = fopen("./Inc_Apply_to_En.txt", "w");
	pf_debug2 = fopen("./Inc_Apply_to_Hn.txt", "w");
	fclose(pf_debug1);
	fclose(pf_debug2);
}


void TFDTD_Boundary_TFSF::Electric_Field_update()
{
	FILE* pfInc;
	pfInc = fopen("Incident_field_plane.dat", "a");

	//double dds = sqrt((double)((dx*dx + dy*dy + dz*dz) / 3)); //1.0/sqrt((double)(1.0/dx/dx + 1.0/dy/dy + 1.0/dz/dz));
	// double dds = (double)sqrt(sin(Theta)*cos(Phi)*dx*sin(Theta)*cos(Phi)*dx
	// 					    + sin(Theta)*sin(Phi)*dy*sin(Theta)*sin(Phi)*dy
	// 					    + cos(Theta)*dz*cos(Theta)*dz);

	// Mur abs initial
	for (int ii = 0; ii < 2; ii++)
		E_inc_oblic_mur[ii] = E_inc_oblic[Nststep - 2 + ii];

	// normal fdtd
	for (int ii = 2; ii < Nststep - 1; ii++)
		E_inc_oblic[ii + 1] = E_inc_oblic[ii + 1] + dt / (v_ratio) / eps_0 * (H_inc_oblic[ii] - H_inc_oblic[ii + 1]) / dds;

	// Mur abs apply
	E_inc_oblic[Nststep - 1] = abs_Gamma_media * E_inc_oblic_mur[1] - abs_Gamma_media * E_inc_oblic[Nststep - 2] + E_inc_oblic_mur[0];

	E_inc_oblic[6] = Pulse;

	fprintf(pfInc, "%e\n", E_inc_oblic[7]);
	fclose(pfInc);
}

void TFDTD_Boundary_TFSF::Magnetic_Field_update()
{
	//double dds = sqrt((double)((dx*dx + dy*dy + dz*dz) / 3)); //1.0/sqrt((double)(1.0/dx/dx + 1.0/dy/dy + 1.0/dz/dz));
	// double dds = (double)sqrt(sin(Theta)*cos(Phi)*dx*sin(Theta)*cos(Phi)*dx
	// 					    + sin(Theta)*sin(Phi)*dy*sin(Theta)*sin(Phi)*dy
	// 					    + cos(Theta)*dz*cos(Theta)*dz);

	// Mur abs initial
	for (int ii = 2; ii < Nststep; ii++)
		H_inc_oblic[ii] = H_inc_oblic[ii] + dt / (v_ratio) / mu_0 * (E_inc_oblic[ii] - E_inc_oblic[ii + 1]) / dds;
}

void TFDTD_Boundary_TFSF::Calcuate_inc_Field()
{
	double DD = 0;       // Distance
	double DDPP = 0;     // Distance from source corner
	int IFIX_D = 0;   // Gause Distance
	double DDP = 0;      // Distance from Gause distance
	double E_inc_grid;
	double H_inc_grid;

	for (int ii = B_Tx_ini; ii <= B_Tx_ter; ii++)
		for (int jj = B_Ty_ini; jj <= B_Ty_ter; jj++)
		{
			if (ii <= B_Tx_ter - 1)
			{
				// Ex_inc (i+1/2, j, k0)
				DD = Calculate_Distance(ii + 0.5, jj, B_Tz_ini, Theta, Phi) - DD_center;
				DDPP = DD;
				IFIX_D = (int)(DDPP);
				DDP = DDPP - IFIX_D;
				E_inc_grid = (1 - DDP)*E_inc_oblic[IFIX_D + disd] + DDP * E_inc_oblic[IFIX_D + 1 + disd];
				Ex_xyN_correction[ii - B_Tx_ini][jj - B_Ty_ini] = E_inc_grid * Vector_component_Exinc(Polar_psi, Phi, Theta);
				
				// Ex_inc (i+1/2, j, k1)
				DD = Calculate_Distance(ii + 0.5, jj, B_Tz_ter, Theta, Phi) - DD_center;
				DDPP = DD;
				IFIX_D = (int)(DDPP);
				DDP = DDPP - IFIX_D;
				E_inc_grid = (1 - DDP)*E_inc_oblic[IFIX_D + disd] + DDP * E_inc_oblic[IFIX_D + 1 + disd];
				Ex_xyP_correction[ii - B_Tx_ini][jj - B_Ty_ini] = E_inc_grid * Vector_component_Exinc(Polar_psi, Phi, Theta);
			}
			if (jj <= B_Ty_ter - 1)
			{
				// Ey_inc (i, j+1/2, k0)
				DD = Calculate_Distance(ii, jj + 0.5, B_Tz_ini, Theta, Phi) - DD_center;
				DDPP = DD;
				IFIX_D = (int)(DDPP);
				DDP = DDPP - IFIX_D;
				E_inc_grid = (1 - DDP)*E_inc_oblic[IFIX_D + disd] + DDP * E_inc_oblic[IFIX_D + 1 + disd];
				Ey_xyN_correction[ii - B_Tx_ini][jj - B_Ty_ini] = E_inc_grid * Vector_component_Eyinc(Polar_psi, Phi, Theta);
				

				// Ey_inc (i, j+1/2, k0)
				DD = Calculate_Distance(ii, jj + 0.5, B_Tz_ter, Theta, Phi) - DD_center;
				DDPP = DD;
				IFIX_D = (int)(DDPP);
				DDP = DDPP - IFIX_D;
				E_inc_grid = (1 - DDP)*E_inc_oblic[IFIX_D + disd] + DDP * E_inc_oblic[IFIX_D + 1 + disd];
				Ey_xyP_correction[ii - B_Tx_ini][jj - B_Ty_ini] = E_inc_grid * Vector_component_Eyinc(Polar_psi, Phi, Theta);
				
			}
		}

	for (int jj = B_Ty_ini; jj <= B_Ty_ter; jj++)
		for (int kk = B_Tz_ini; kk <= B_Tz_ter; kk++)
		{
			if (jj <= B_Ty_ter - 1)
			{
				// Ex_inc (i0, j+1/2, k0)
				DD = Calculate_Distance(B_Tx_ini, jj + 0.5, kk, Theta, Phi) - DD_center;
				DDPP = DD;
				IFIX_D = (int)(DDPP);
				DDP = DDPP - IFIX_D;
				E_inc_grid = (1 - DDP)*E_inc_oblic[IFIX_D + disd] + DDP * E_inc_oblic[IFIX_D + 1 + disd];
				Ey_yzN_correction[jj - B_Ty_ini][kk - B_Tz_ini] = E_inc_grid * Vector_component_Eyinc(Polar_psi, Phi, Theta);

				// Ex_inc (10, j+1/2, k0)
				DD = Calculate_Distance(B_Tx_ter, jj + 0.5, kk, Theta, Phi) - DD_center;
				DDPP = DD;
				IFIX_D = (int)(DDPP);
				DDP = DDPP - IFIX_D;
				E_inc_grid = (1 - DDP)*E_inc_oblic[IFIX_D + disd] + DDP * E_inc_oblic[IFIX_D + 1 + disd];
				Ey_yzP_correction[jj - B_Ty_ini][kk - B_Tz_ini] = E_inc_grid * Vector_component_Eyinc(Polar_psi, Phi, Theta);
			}
			if (kk <= B_Tz_ter - 1)
			{
				// Ey_inc (i0, j, k+1/2)
				DD = Calculate_Distance(B_Tx_ini, jj, kk + 0.5, Theta, Phi) - DD_center;
				DDPP = DD;
				IFIX_D = (int)(DDPP);
				DDP = DDPP - IFIX_D;
				E_inc_grid = (1 - DDP)*E_inc_oblic[IFIX_D + disd] + DDP * E_inc_oblic[IFIX_D + 1 + disd];
				Ez_yzN_correction[jj - B_Ty_ini][kk - B_Tz_ini] = E_inc_grid * Vector_component_Ezinc(Polar_psi, Phi, Theta);

				// Ey_inc (i1, j, k+1/2)
				DD = Calculate_Distance(B_Tx_ter, jj, kk + 0.5, Theta, Phi) - DD_center;
				DDPP = DD;
				IFIX_D = (int)(DDPP);
				DDP = DDPP - IFIX_D;
				E_inc_grid = (1 - DDP)*E_inc_oblic[IFIX_D + disd] + DDP * E_inc_oblic[IFIX_D + 1 + disd];
				Ez_yzP_correction[jj - B_Ty_ini][kk - B_Tz_ini] = E_inc_grid * Vector_component_Ezinc(Polar_psi, Phi, Theta);
			}
		}

	for (int ii = B_Tx_ini; ii <= B_Tx_ter; ii++)
		for (int kk = B_Tz_ini; kk <= B_Tz_ter; kk++)
		{
			if (ii <= B_Tx_ter - 1)
			{
				// Ex_inc (i+1/2, j0, k)
				DD = Calculate_Distance(ii + 0.5, B_Ty_ini, kk, Theta, Phi) - DD_center;
				DDPP = DD;
				IFIX_D = (int)(DDPP);
				DDP = DDPP - IFIX_D;
				E_inc_grid = (1 - DDP)*E_inc_oblic[IFIX_D + disd] + DDP * E_inc_oblic[IFIX_D + 1 + disd];
				Ex_xzN_correction[ii - B_Tx_ini][kk - B_Tz_ini] = E_inc_grid  * Vector_component_Exinc(Polar_psi, Phi, Theta);

				// Ex_inc (i+1/2, j1, k)
				DD = Calculate_Distance(ii + 0.5, B_Ty_ter, kk, Theta, Phi) - DD_center;
				DDPP = DD;
				IFIX_D = (int)(DDPP);
				DDP = DDPP - IFIX_D;
				E_inc_grid = (1 - DDP)*E_inc_oblic[IFIX_D + disd] + DDP * E_inc_oblic[IFIX_D + 1 + disd];
				Ex_xzP_correction[ii - B_Tx_ini][kk - B_Tz_ini] = E_inc_grid * Vector_component_Exinc(Polar_psi, Phi, Theta);
			}
			if (kk <= B_Tz_ter - 1)
			{
				// Ex_inc (i, j0, k+1/2)
				DD = Calculate_Distance(ii, B_Ty_ini, kk + 0.5, Theta, Phi) - DD_center;
				DDPP = DD;
				IFIX_D = (int)(DDPP);
				DDP = DDPP - IFIX_D;
				E_inc_grid = (1 - DDP)*E_inc_oblic[IFIX_D + disd] + DDP * E_inc_oblic[IFIX_D + 1 + disd];
				Ez_xzN_correction[ii - B_Tx_ini][kk - B_Tz_ini] = E_inc_grid * Vector_component_Ezinc(Polar_psi, Phi, Theta);

				// Ex_inc (i, j1, k+1/2)
				DD = Calculate_Distance(ii, B_Ty_ter, kk + 0.5, Theta, Phi) - DD_center;
				DDPP = DD;
				IFIX_D = (int)(DDPP);
				DDP = DDPP - IFIX_D;
				E_inc_grid = (1 - DDP)*E_inc_oblic[IFIX_D + disd] + DDP * E_inc_oblic[IFIX_D + 1 + disd];
				Ez_xzP_correction[ii - B_Tx_ini][kk - B_Tz_ini] = E_inc_grid * Vector_component_Ezinc(Polar_psi, Phi, Theta);
			}
		}

	for (int jj = B_Ty_ini; jj <= B_Ty_ter; jj++)
		for (int kk = B_Tz_ini; kk <= B_Tz_ter; kk++)
		{
			if (jj <= B_Ty_ter - 1)
			{
				// Hz_inc (i0-1/2, j+1/2, k)
				DD = Calculate_Distance(B_Tx_ini - 0.5, jj + 0.5, kk, Theta, Phi) - DD_center_H;
				DDPP = DD;// + 0.5;
				IFIX_D = (int)(DDPP);
				DDP = DDPP - IFIX_D;
				H_inc_grid = (1 - DDP)*H_inc_oblic[IFIX_D - 1 + disd] + DDP * H_inc_oblic[IFIX_D + disd];
				Hz_yzN_correction[jj - B_Ty_ini][kk - B_Tz_ini] = H_inc_grid * Vector_component_Hzinc(Polar_psi, Phi, Theta);

				// Hz_inc (i1-1/2, j+1/2, k)
				DD = Calculate_Distance((double)B_Tx_ter - 0.5, jj + 0.5, kk, Theta, Phi) - DD_center_H;
				DDPP = DD;// + 0.5;
				IFIX_D = (int)(DDPP);
				DDP = DDPP - IFIX_D;
				H_inc_grid = (1 - DDP)*H_inc_oblic[IFIX_D - 1 + disd] + DDP * H_inc_oblic[IFIX_D + disd];
				Hz_yzP_correction[jj - B_Ty_ini][kk - B_Tz_ini] = H_inc_grid * Vector_component_Hzinc(Polar_psi, Phi, Theta);
			}
			if (kk <= B_Tz_ter - 1)
			{
				// Hy_inc (i0-1/2, j, k+1/2)
				DD = Calculate_Distance((double)B_Tx_ini - 0.5, jj, kk + 0.5, Theta, Phi) - DD_center_H;
				DDPP = DD;// + 0.5;
				IFIX_D = (int)(DDPP);
				DDP = DDPP - IFIX_D;
				H_inc_grid = (1 - DDP)*H_inc_oblic[IFIX_D - 1 + disd] + DDP * H_inc_oblic[IFIX_D + disd];
				Hy_yzN_correction[jj - B_Ty_ini][kk - B_Tz_ini] = H_inc_grid * Vector_component_Hyinc(Polar_psi, Phi, Theta);

				// Hy_inc (i1-1/2, j, k+1/2)
				DD = Calculate_Distance((double)B_Tx_ter - 0.5, jj, kk + 0.5, Theta, Phi) - DD_center_H;
				DDPP = DD;// + 0.5;
				IFIX_D = (int)(DDPP);
				DDP = DDPP - IFIX_D;
				H_inc_grid = (1 - DDP)*H_inc_oblic[IFIX_D - 1 + disd] + DDP * H_inc_oblic[IFIX_D + disd];
				Hy_yzP_correction[jj - B_Ty_ini][kk - B_Tz_ini] = H_inc_grid * Vector_component_Hyinc(Polar_psi, Phi, Theta);
			}
		}

	for (int ii = B_Tx_ini; ii <= B_Tx_ter; ii++)
		for (int kk = B_Tz_ini; kk <= B_Tz_ter; kk++)
		{
			if (kk <= B_Tz_ter - 1)
			{
				// Hx_inc (i, j0-1/2, k+1/2)
				DD = Calculate_Distance(ii, (double)B_Ty_ini - 0.5, kk + 0.5, Theta, Phi) - DD_center_H;
				DDPP = DD;// + 0.5;
				IFIX_D = (int)(DDPP);
				DDP = DDPP - IFIX_D;
				H_inc_grid = (1 - DDP)*H_inc_oblic[IFIX_D - 1 + disd] + DDP * H_inc_oblic[IFIX_D + disd];
				Hx_xzN_correction[ii - B_Tx_ini][kk - B_Tz_ini] = H_inc_grid * Vector_component_Hxinc(Polar_psi, Phi, Theta);

				// Hx_inc (i, j1-1/2, k+1/2)
				DD = Calculate_Distance(ii, (double)B_Ty_ter - 0.5, kk + 0.5, Theta, Phi) - DD_center_H;
				DDPP = DD;// + 0.5;
				IFIX_D = (int)(DDPP);
				DDP = DDPP - IFIX_D;
				H_inc_grid = (1 - DDP)*H_inc_oblic[IFIX_D - 1 + disd] + DDP * H_inc_oblic[IFIX_D + disd];
				Hx_xzP_correction[ii - B_Tx_ini][kk - B_Tz_ini] = H_inc_grid * Vector_component_Hxinc(Polar_psi, Phi, Theta);
			}
			if (ii <= B_Tx_ter - 1)
			{
				// Hz_inc (i+1/2, j0-1/2, k)
				DD = Calculate_Distance(ii + 0.5, (double)B_Ty_ini - 0.5, kk, Theta, Phi) - DD_center_H;
				DDPP = DD;// + 0.5;
				IFIX_D = (int)(DDPP);
				DDP = DDPP - IFIX_D;
				H_inc_grid = (1 - DDP)*H_inc_oblic[IFIX_D - 1 + disd] + DDP * H_inc_oblic[IFIX_D + disd];
				Hz_xzN_correction[ii - B_Tx_ini][kk - B_Tz_ini] = H_inc_grid * Vector_component_Hzinc(Polar_psi, Phi, Theta);

				// Hz_inc (i+1/2, j1-1/2, k)
				DD = Calculate_Distance(ii + 0.5, (double)B_Ty_ter - 0.5, kk, Theta, Phi) - DD_center_H;
				DDPP = DD;// + 0.5;
				IFIX_D = (int)(DDPP);
				DDP = DDPP - IFIX_D;
				H_inc_grid = (1 - DDP)*H_inc_oblic[IFIX_D - 1 + disd] + DDP * H_inc_oblic[IFIX_D + disd];
				Hz_xzP_correction[ii - B_Tx_ini][kk - B_Tz_ini] = H_inc_grid * Vector_component_Hzinc(Polar_psi, Phi, Theta);
			}
		}

	for (int ii = B_Tx_ini; ii <= B_Tx_ter; ii++)
		for (int jj = B_Ty_ini; jj <= B_Ty_ter; jj++)
		{
			if (jj <= B_Ty_ter - 1)
			{
				// Hx_inc (i, j+1/2, k0-1/2)
				DD = Calculate_Distance(ii, jj + 0.5, (double)B_Tz_ini - 0.5, Theta, Phi) - DD_center_H;
				DDPP = DD;// + 0.5;
				IFIX_D = (int)(DDPP);
				DDP = DDPP - IFIX_D;
				H_inc_grid = (1 - DDP)*H_inc_oblic[IFIX_D - 1 + disd] + DDP * H_inc_oblic[IFIX_D + disd];
				Hx_xyN_correction[ii - B_Tx_ini][jj - B_Ty_ini] = H_inc_grid * Vector_component_Hxinc(Polar_psi, Phi, Theta);

				// Hx_inc (i, j+1/2, k1-1/2)
				DD = Calculate_Distance(ii, jj + 0.5, (double)B_Tz_ter - 0.5, Theta, Phi) - DD_center_H;
				DDPP = DD;// + 0.5;
				IFIX_D = (int)(DDPP);
				DDP = DDPP - IFIX_D;
				H_inc_grid = (1 - DDP)*H_inc_oblic[IFIX_D - 1 + disd] + DDP * H_inc_oblic[IFIX_D + disd];
				Hx_xyP_correction[ii - B_Tx_ini][jj - B_Ty_ini] = H_inc_grid * Vector_component_Hxinc(Polar_psi, Phi, Theta);
			}
			if (ii <= B_Tx_ter - 1)
			{
				// Hy_inc  (i+1/2, j, k0-1/2)
				DD = Calculate_Distance(ii + 0.5, jj, (double)B_Tz_ini - 0.5, Theta, Phi) - DD_center_H;
				DDPP = DD;// + 0.5;
				IFIX_D = (int)(DDPP);
				DDP = DDPP - IFIX_D;
				H_inc_grid = (1 - DDP)*H_inc_oblic[IFIX_D - 1 + disd] + DDP * H_inc_oblic[IFIX_D + disd];
				Hy_xyN_correction[ii - B_Tx_ini][jj - B_Ty_ini] = H_inc_grid * Vector_component_Hyinc(Polar_psi, Phi, Theta);

				// Hy_inc  (i+1/2, j, k1-1/2)
				DD = Calculate_Distance(ii + 0.5, jj, (double)B_Tz_ter - 0.5, Theta, Phi) - DD_center_H;
				DDPP = DD;// + 0.5;
				IFIX_D = (int)(DDPP);
				DDP = DDPP - IFIX_D;
				H_inc_grid = (1 - DDP)*H_inc_oblic[IFIX_D - 1 + disd] + DDP * H_inc_oblic[IFIX_D + disd];
				Hy_xyP_correction[ii - B_Tx_ini][jj - B_Ty_ini] = H_inc_grid * Vector_component_Hyinc(Polar_psi, Phi, Theta);
			}
		}
}

///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////		E - Field 	Correction 			  /////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////

void TFDTD_Boundary_TFSF::Correction_X_Boundary_Exn()
{
	// x- correction
	for (int jj = B_Ty_ini; jj < B_Ty_ter; jj++)
		for (int kk = B_Tz_ini; kk <= B_Tz_ter; kk++)
			// Ey( i = i0; j = j0 + 1/2, ... , j1 - 1/2; k = k0, ... k1)
			// Ey|(n+1)(i0,j,k) = ~ + dt/eps0/dx * Hz_xyN_corr->(i0-1, j, k)
			Ey[B_Tx_ini][jj][kk] = Ey[B_Tx_ini][jj][kk] + Cedx1 * Hz_yzN_correction[jj - B_Ty_ini][kk - B_Tz_ini];
	//*(*(*(*EPy + B_Tx_ini) + jj) + kk) = *(*(*(*EPy + B_Tx_ini) + jj) + kk) + dt / eps_0 / dx * Hz_yzN_correction[jj - B_Ty_ini][kk - B_Tz_ini];


	for (int jj = B_Ty_ini; jj <= B_Ty_ter; jj++)
		for (int kk = B_Tz_ini; kk < B_Tz_ter; kk++)
			// Ez( i = i0; j = j0, ... , j1; k = k0+1/2, ... k1-1/2)
			// Ez|(n+1)(i0, j = j0, ... , j1, k =  k0 + 1/2, ... , k1 - 1/2
			Ez[B_Tx_ini][jj][kk] = Ez[B_Tx_ini][jj][kk] - Cedx1 * Hy_yzN_correction[jj - B_Ty_ini][kk - B_Tz_ini];
	//*(*(*(*EPz + B_Tx_ini) + jj) + kk) = *(*(*(*EPz + B_Tx_ini) + jj) + kk) - dt / eps_0 / dx * Hy_yzN_correction[jj - B_Ty_ini][kk - B_Tz_ini];
}

void TFDTD_Boundary_TFSF::Correction_X_Boundary_Exp()
{
	// x+ correction
	for (int jj = B_Ty_ini; jj < B_Ty_ter; jj++)
		for (int kk = B_Tz_ini; kk <= B_Tz_ter; kk++)
			// Ey( i = i1; j = j0 + 1/2, ... , j1 - 1/2; k = k0, ..., k1)
			// Ey|(n+1)(i1,j,k) = ~ + dt/eps0/dx * Hz_xyN_corr->(i0-1, j, k)
			Ey[B_Tx_ter][jj][kk] = Ey[B_Tx_ter][jj][kk] - Cedx2 * Hz_yzP_correction[jj - B_Ty_ini][kk - B_Tz_ini];
	//*(*(*(*EPy + B_Tx_ter) + jj) + kk) = *(*(*(*EPy + B_Tx_ter) + jj) + kk) - (dt / eps_0 / dx) * Hz_yzP_correction[jj - B_Ty_ini][kk - B_Tz_ini];

	// x+ correction
	for (int jj = B_Ty_ini; jj <= B_Ty_ter; jj++)
		for (int kk = B_Tz_ini; kk < B_Tz_ter; kk++)
			// Ez( i = i0; j = j0, ... , j1; k = k0+1/2, ... k1-1/2)
			// Ez|(n+1)(i1, j = j0, ... , j1, k =  k0 + 1/2, ... , k1 - 1/2
			Ez[B_Tx_ter][jj][kk] = Ez[B_Tx_ter][jj][kk] + Cedx2 * Hy_yzP_correction[jj - B_Ty_ini][kk - B_Tz_ini];
	//*(*(*(*EPz + B_Tx_ter) + jj) + kk) = *(*(*(*EPz + B_Tx_ter) + jj) + kk) + (dt / eps_0 / dx) * Hy_yzP_correction[jj - B_Ty_ini][kk - B_Tz_ini];
}

void TFDTD_Boundary_TFSF::Correction_Y_Boundary_Eyn()
{
	// y- correction
	for (int ii = B_Tx_ini; ii < B_Tx_ter; ii++)
		for (int kk = B_Tz_ini; kk <= B_Tz_ter; kk++)
		{
			// Ex( i = i0 ... i1; j = j0; k = k0+1/2, ... k1-1/2)
			// Ex|(n+1)(i0 + 1/2, ... , i1-1/2, j = j0, k =  k0, ... , k1)
			//*(*(*(*EPx + ii) + B_Ty_ini) + kk) = *(*(*(*EPx + ii) + B_Ty_ini) + kk) - dt / eps_0 / dy * Hz_xzN_correction[ii - B_Tx_ini][kk - B_Tz_ini];
			Ex[ii][B_Ty_ini][kk] = Ex[ii][B_Ty_ini][kk] - Cedy1 * Hz_xzN_correction[ii - B_Tx_ini][kk - B_Tz_ini];
		}

	for (int ii = B_Tx_ini; ii <= B_Tx_ter; ii++)
		for (int kk = B_Tz_ini; kk < B_Tz_ter; kk++)
		{
			// Ez(n+1)(i,j0,k) = ~ + dt/eps0/dy * Hx_zxN_corr->(i0-1, j, k)
			//*(*(*(*EPz + ii) + B_Ty_ini) + kk) = *(*(*(*EPz + ii) + B_Ty_ini) + kk) + dt / eps_0 / dy * Hx_xzN_correction[ii - B_Tx_ini][kk - B_Tz_ini];
			Ez[ii][B_Ty_ini][kk] = Ez[ii][B_Ty_ini][kk] + Cedy1 * Hx_xzN_correction[ii - B_Tx_ini][kk - B_Tz_ini];
		}


}

void TFDTD_Boundary_TFSF::Correction_Y_Boundary_Eyp()
{
	// y+ correction
	for (int ii = B_Tx_ini; ii < B_Tx_ter; ii++)
		for (int kk = B_Tz_ini; kk <= B_Tz_ter; kk++)
		{
			// Ex( i = i0 ... i1; j = j0; k = k0+1/2, ... k1-1/2)
			// Ex|(n+1)(i0 + 1/2, ... , i1-1/2, j = j0, k =  k0, ... , k1)
			//*(*(*(*EPx + ii) + B_Ty_ter) + kk) = *(*(*(*EPx + ii) + B_Ty_ter) + kk) + dt / eps_0 / dy * Hz_xzP_correction[ii - B_Tx_ini][kk - B_Tz_ini];
			Ex[ii][B_Ty_ter][kk] = Ex[ii][B_Ty_ter][kk] + Cedy2 * Hz_xzP_correction[ii - B_Tx_ini][kk - B_Tz_ini];
		}

	for (int ii = B_Tx_ini; ii <= B_Tx_ter; ii++)
		for (int kk = B_Tz_ini; kk < B_Tz_ter; kk++)
		{
			// Ez(n+1)(i,j0,k) = ~ + dt/eps0/dy * Hx_zxN_corr->(i0-1, j, k)
			//*(*(*(*EPz + ii) + B_Ty_ter) + kk) = *(*(*(*EPz + ii) + B_Ty_ter) + kk) - dt / eps_0 / dy * Hx_xzP_correction[ii - B_Tx_ini][kk - B_Tz_ini];
			Ez[ii][B_Ty_ter][kk] = Ez[ii][B_Ty_ter][kk] - Cedy2 * Hx_xzP_correction[ii - B_Tx_ini][kk - B_Tz_ini];
		}
}

void TFDTD_Boundary_TFSF::Correction_Z_Boundary_Ezn()
{
	// z- correction
	for (int ii = B_Tx_ini; ii < B_Tx_ter; ii++)
		for (int jj = B_Ty_ini; jj <= B_Ty_ter; jj++)
		{
			// Ex( i = i0 +1/2 ... i1 -1/2; j = j0 ... j1; k = k0)
			// Ex|(n+1)(i0 + 1/2, ... , i1-1/2, j = j0, k =  k0) = ~ - dt/eps0/dy * Hy_zxN_corr->(i, j, k0-1/2)
			//*(*(*(*EPx + ii) + jj) + B_Tz_ini) = *(*(*(*EPx + ii) + jj) + B_Tz_ini) + dt / eps_0 / dz * Hy_xyN_correction[ii - B_Tx_ini][jj - B_Ty_ini];
			Ex[ii][jj][B_Tz_ini] = Ex[ii][jj][B_Tz_ini] + Cedz1 * Hy_xyN_correction[ii - B_Tx_ini][jj - B_Ty_ini];
		}

	for (int ii = B_Tx_ini; ii <= B_Tx_ter; ii++)
		for (int jj = B_Ty_ini; jj < B_Ty_ter; jj++)
		{
			// Ez(n+1)(i0, ... , i1; j= j0+1/2, ... , j1 - 1/2, k=k0) = ~ + dt/eps0/dy * Hx_zxN_corr->(i, j, k0-1/2)
			//*(*(*(*EPy + ii) + jj) + B_Tz_ini) = *(*(*(*EPy + ii) + jj) + B_Tz_ini) - dt / eps_0 / dz * Hx_xyN_correction[ii - B_Tx_ini][jj - B_Ty_ini];
			Ey[ii][jj][B_Tz_ini] = Ey[ii][jj][B_Tz_ini] - Cedz1 * Hx_xyN_correction[ii - B_Tx_ini][jj - B_Ty_ini];
		}
	// debug code
	//        int ci = (B_Tx_ini + B_Tx_ter)/2;
	//        int cj = (B_Ty_ini + B_Ty_ter)/2;
	//        pf_debug1 = fopen("./Inc_Apply_to_En.txt", "a");
	//        fprintf(pf_debug1, "%le\t%le\n", dt / dz * Hx_xyN_correction[ci][cj], dt / dz * Hy_xyN_correction[ci][cj]);
	//        fclose(pf_debug1);
}

void TFDTD_Boundary_TFSF::Correction_Z_Boundary_Ezp()
{
	// z+ correction
	for (int ii = B_Tx_ini; ii < B_Tx_ter; ii++)
	for (int jj = B_Ty_ini; jj <= B_Ty_ter; jj++)
	{
		// Ex( i = i0 +1/2 ... i1 -1/2; j = j0 ... j1; k = k0)
		// Ex|(n+1)(i0 + 1/2, ... , i1-1/2, j = j0, k =  k0) = ~ - dt/eps0/dy * Hy_zxN_corr->(i, j, k0-1/2)
		//*(*(*(*EPx + ii) + jj) + B_Tz_ter) = *(*(*(*EPx + ii) + jj) + B_Tz_ter) - dt / eps_0 / dz * Hy_xyP_correction[ii - B_Tx_ini][jj - B_Ty_ini];
		Ex[ii][jj][B_Tz_ter] = Ex[ii][jj][B_Tz_ter] - Cedz2 * Hy_xyP_correction[ii - B_Tx_ini][jj - B_Ty_ini];
	}

	for (int ii = B_Tx_ini; ii <= B_Tx_ter; ii++)
	for (int jj = B_Ty_ini; jj < B_Ty_ter; jj++)
	{
		// Ez(n+1)(i0, ... , i1; j= j0+1/2, ... , j1 - 1/2, k=k0) = ~ + dt/eps0/dy * Hx_zxN_corr->(i, j, k0-1/2)
		//*(*(*(*EPy + ii) + jj) + B_Tz_ter) = *(*(*(*EPy + ii) + jj) + B_Tz_ter) + dt / eps_0 / dz * Hx_xyP_correction[ii - B_Tx_ini][jj - B_Ty_ini];
		Ey[ii][jj][B_Tz_ter] = Ey[ii][jj][B_Tz_ter] + Cedz2 * Hx_xyP_correction[ii - B_Tx_ini][jj - B_Ty_ini];
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////		H - Field 	Correction 			  ////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

void TFDTD_Boundary_TFSF::Correction_X_Boundary_Hxn()
{
	// x- correction i0 -1/2 face
	for (int jj = B_Ty_ini; jj <= B_Ty_ter; jj++)
		for (int kk = B_Tz_ini; kk <= B_Tz_ter; kk++)
		{

			if (jj < B_Ty_ter) // jj = j0+1/2 .. j1 - 1/2
			{
				// Hz(i = i0-1/2 ; j= j0 +1/2, ... , j1 -1/2 ; k =  k0, ... , k1 )
				// Hz|(n+1/2)(i0 -1/2,j,k) = ~ + dt/mu0/dx * Ey_xyN_corr->(i0, j, k)
				Hz[B_Tx_ini - 1][jj][kk] = Hz[B_Tx_ini - 1][jj][kk] + Chdx1 * Ey_yzN_correction[jj - B_Ty_ini][kk - B_Tz_ini];
				//				*(*(*(*HPz + B_Tx_ini) + jj) + kk) = *(*(*(*HPz + B_Tx_ini) + jj) + kk) + dt / mu_0 / dx * Ey_yzN_correction[jj - B_Ty_ini][kk - B_Tz_ini];
			}

			if (kk < B_Tz_ter) // kk = k0+1/2 ... k1-1/2
			{
				// Hy(i = i0-1/2 ; j= j0, ... , j1, k =  k0 + 1/2, ... , k1 - 1/2 )
				// Hy|(n+1/2)(i0 -1/2, j, k) = ~ -dt/mu0/dx * Ez_xyN_corr->(i0, j, k)
				Hy[B_Tx_ini - 1][jj][kk] = Hy[B_Tx_ini - 1][jj][kk] - Chdx1 * Ez_yzN_correction[jj - B_Ty_ini][kk - B_Tz_ini];
				//				*(*(*(*HPy + B_Tx_ini) + jj) + kk) = *(*(*(*HPy + B_Tx_ini) + jj) + kk) - dt / mu_0 / dx * Ez_yzN_correction[jj - B_Ty_ini][kk - B_Tz_ini];
			}
		}
}

void TFDTD_Boundary_TFSF::Correction_X_Boundary_Hxp()
{
	// x+ correction
	for (int jj = B_Ty_ini; jj <= B_Ty_ter; jj++)
		for (int kk = B_Tz_ini; kk <= B_Tz_ter; kk++)
		{
			if (jj < B_Ty_ter) // jj = j0+1/2 .. j1 - 1/2
			{
				// Hz(i = i1+1/2 ; j= j0 +1/2, ... , j1 -1/2 ; k =  k0, ... , k1 )
				// Hz|(n+1/2)(i1 +1/2,j,k) = ~ - dt/mu0/dx * Ey_xyP_corr->(i0, j, k)
				Hz[B_Tx_ter][jj][kk] = Hz[B_Tx_ter][jj][kk] - Chdx2 * Ey_yzP_correction[jj - B_Ty_ini][kk - B_Tz_ini];
				//				*(*(*(*HPz + B_Tx_ter-1) + jj) + kk) = *(*(*(*HPz + B_Tx_ter-1) + jj) + kk) - dt / mu_0 / dx * Ey_yzP_correction[jj - B_Ty_ini][kk - B_Tz_ini];
			}

			if (kk < B_Tz_ter) // kk = k0+1/2 ... k1-1/2
			{
				// Hy(i = i1+1/2 ; j= j0, ... , j1, k =  k0 + 1/2, ... , k1 - 1/2 )
				// Hy|(n+1/2)(i1 +1/2, j, k) = ~ -dt/mu0/dx * Ez_xyN_corr->(i0, j, k)
				Hy[B_Tx_ter][jj][kk] = Hy[B_Tx_ter][jj][kk] + Chdx2 * Ez_yzP_correction[jj - B_Ty_ini][kk - B_Tz_ini];
				//				*(*(*(*HPy + B_Tx_ter-1) + jj) + kk) = *(*(*(*HPy + B_Tx_ter-1) + jj) + kk) + dt / mu_0 / dx * Ez_yzP_correction[jj - B_Ty_ini][kk - B_Tz_ini];
			}
		}
}

void TFDTD_Boundary_TFSF::Correction_Y_Boundary_Hyn()
{
	// y- correction : j = j0 - 1/2
	for (int ii = B_Tx_ini; ii <= B_Tx_ter; ii++)
		for (int kk = B_Tz_ini; kk <= B_Tz_ter; kk++)
		{
			if (ii < B_Tx_ter)
			{
				// Hz( i = i0+1/2 ... i1 -1/2; j = j0- 1/2; k = k0, ... k1)
				// Hz|(n+1/2)(i; j = j0-1/2 ; k) = ~ - dt/mu0/dy * Ex
				Hz[ii][B_Ty_ini - 1][kk] = Hz[ii][B_Ty_ini - 1][kk] - Chdy1 * Ex_xzN_correction[ii - B_Tx_ini][kk - B_Tz_ini];
				//				*(*(*(*HPz + ii) + B_Ty_ini) + kk) = *(*(*(*HPz + ii) + B_Ty_ini) + kk) - dt / mu_0 / dy * Ex_xzN_correction[ii - B_Tx_ini][kk - B_Tz_ini];
			}
			if (kk < B_Tz_ter)
			{
				// Hx( i = i0, ... , i1; j = j0-1/2 ; k = k0+1/2, ..., k1 - 1/2)
				//			// Hx(n+1.2)(i,j0,k) = ~ + dt/eps0/dy * Hx_zxN_corr->(i0-1, j, k)
				Hx[ii][B_Ty_ini - 1][kk] = Hx[ii][B_Ty_ini - 1][kk] + Chdy1 * Ez_xzN_correction[ii - B_Tx_ini][kk - B_Tz_ini];
				//				*(*(*(*HPx + ii) + B_Ty_ini) + kk) = *(*(*(*HPx + ii) + B_Ty_ini) + kk) + dt / mu_0 / dy * Ez_xzN_correction[ii - B_Tx_ini][kk - B_Tz_ini];
			}
		}
}

void TFDTD_Boundary_TFSF::Correction_Y_Boundary_Hyp()
{
	// y+ correction
	for (int ii = B_Tx_ini; ii <= B_Tx_ter; ii++)
		for (int kk = B_Tz_ini; kk <= B_Tz_ter; kk++)
		{
			if (ii < B_Tx_ter)
			{
				// Hz( i = i0+1/2, ..., i1 -1/2 ; j = j1+ 1/2 ; k = k0, ... k1)
				// Hz|(n+1/2)(i; j = j1-1/2 ; k) = ~ + dt/mu0/dy * Ex
				Hz[ii][B_Ty_ter][kk] = Hz[ii][B_Ty_ter][kk] + Chdy2 * Ex_xzP_correction[ii - B_Tx_ini][kk - B_Tz_ini];
				//				*(*(*(*HPz + ii) + B_Ty_ter-1) + kk) = *(*(*(*HPz + ii) + B_Ty_ter-1) + kk) + dt / mu_0 / dy * Ex_xzP_correction[ii - B_Tx_ini][kk - B_Tz_ini];
			}
			if (kk < B_Tz_ter)
			{
				// Hx( i = i0, ... , i1 ; j = j1 + 1/2 ; k = k0+1/2, ..., k1 - 1/2)
				//			// Hx(n+1.2)(i,j1 - 1/2,k) = ~ - dt/eps0/dy * Hx_zxN_corr->(i0-1, j, k)
				Hx[ii][B_Ty_ter][kk] = Hx[ii][B_Ty_ter][kk] - Chdy2 * Ez_xzP_correction[ii - B_Tx_ini][kk - B_Tz_ini];
				//				*(*(*(*HPx + ii) + B_Ty_ter-1) + kk) = *(*(*(*HPx + ii) + B_Ty_ter-1) + kk) - dt / mu_0 / dy * Ez_xzP_correction[ii - B_Tx_ini][kk - B_Tz_ini];
			}
		}
}

void TFDTD_Boundary_TFSF::Correction_Z_Boundary_Hzn()
{
	// z- correction
	for (int ii = B_Tx_ini; ii <= B_Tx_ter; ii++)
		for (int jj = B_Ty_ini; jj <= B_Ty_ter; jj++)
		{
			if (ii < B_Tx_ter)
			{
				// Hy( i = i0 +1/2 ... i1 -1/2; j = j0, ..., j1; k = k0 - 1/2)
				// Hy|(n+1/2)(i, j, k0 -1/2) =  ~ + dt/mu0/dz * Ex
				Hy[ii][jj][B_Tz_ini - 1] = Hy[ii][jj][B_Tz_ini - 1] + Chdz1 * Ex_xyN_correction[ii - B_Tx_ini][jj - B_Ty_ini];
				//				*(*(*(*HPy + ii) + jj) + B_Tz_ini) = *(*(*(*HPy + ii) + jj) + B_Tz_ini) + dt / mu_0 / dz * Ex_xyN_correction[ii - B_Tx_ini][jj - B_Ty_ini];
			}
			if (jj < B_Ty_ter)
			{
				// Hx( i = i0, ..., i1; j = j0 + 1/2, ..., j1 -1/2; k = k0 - 1/2)
				// Hx(n+1)(i,j0,k) = ~ - dt/eps0/dy * Ey_zxN_corr->(i0-1, j, k)
				Hx[ii][jj][B_Tz_ini - 1] = Hx[ii][jj][B_Tz_ini - 1] - Chdz1 * Ey_xyN_correction[ii - B_Tx_ini][jj - B_Ty_ini];
				//				*(*(*(*HPx + ii) + jj) + B_Tz_ini) = *(*(*(*HPx + ii) + jj) + B_Tz_ini) - dt / mu_0 / dz * Ey_xyN_correction[ii - B_Tx_ini][jj - B_Ty_ini];
			}
		}
	// debug code
	//       int ci = (B_Tx_ini + B_Tx_ter)/2;
	//       int cj = (B_Ty_ini + B_Ty_ter)/2;
	//       pf_debug2 = fopen("./Inc_Apply_to_Hn.txt", "a");
	//       fprintf(pf_debug2, "%le\t%le\n",dt / mu_0 / dz *Ex_xyN_correction[ci][cj],dt / mu_0 / dz *Ey_xyN_correction[ci][cj]);
	//       fclose(pf_debug2);
}

void TFDTD_Boundary_TFSF::Correction_Z_Boundary_Hzp()
{
	// z+ correction
	for (int ii = B_Tx_ini; ii <= B_Tx_ter; ii++)
		for (int jj = B_Ty_ini; jj <= B_Ty_ter; jj++)
		{
			if (ii < B_Tx_ter)
			{
				// Hy( i = i0 +1/2 ... i1 -1/2; j = j0, ..., j1; k = k0 - 1/2)
				// Hy|(n+1/2)(i, j, k0 -1/2) =  ~ + dt/mu0/dz * Ex
				Hy[ii][jj][B_Tz_ter] = Hy[ii][jj][B_Tz_ter] - Chdz2 * Ex_xyP_correction[ii - B_Tx_ini][jj - B_Ty_ini];
				//				*(*(*(*HPy + ii) + jj) + B_Tz_ter-1) = *(*(*(*HPy + ii) + jj) + B_Tz_ter-1) - dt / mu_0 / dz * Ex_xyP_correction[ii - B_Tx_ini][jj - B_Ty_ini];
			}
			if (jj < B_Ty_ter)
			{
				// Hx( i = i0, ..., i1; j = j0 + 1/2, ..., j1 -1/2; k = k0 - 1/2)
				// Hx(n+1)(i,j0,k) = ~ - dt/eps0/dy * Ey_zxN_corr->(i0-1, j, k)
				Hx[ii][jj][B_Tz_ter] = Hx[ii][jj][B_Tz_ter] + Chdz2 * Ey_xyP_correction[ii - B_Tx_ini][jj - B_Ty_ini];
				//				*(*(*(*HPx + ii) + jj) + B_Tz_ter-1) = *(*(*(*HPx + ii) + jj) + B_Tz_ter-1) + dt / mu_0 / dz * Ey_xyP_correction[ii - B_Tx_ini][jj - B_Ty_ini];
			}
		}
}

void TFDTD_Boundary_TFSF::Calcuate_inc_Field_E_OMP()
{
	double DD = 0;       // Distance
	double DDPP = 0;     // Distance from source corner
	int IFIX_D = 0;   // Gause Distance
	double DDP = 0;      // Distance from Gause distance
	double E_inc_grid;
	double H_inc_grid;
	int ii, jj, kk;
#pragma omp parallel private(ii, jj, kk, DD, DDPP, IFIX_D, DDP, E_inc_grid, H_inc_grid)
	{
#pragma omp for schedule(dynamic)
		for (ii = B_Tx_ini; ii <= B_Tx_ter; ii++)
			for (jj = B_Ty_ini; jj <= B_Ty_ter; jj++)
			{
				if (ii <= B_Tx_ter - 1)
				{
					// Ex_inc (i+1/2, j, k0)
					DD = Calculate_Distance(ii + 0.5, jj, B_Tz_ini, Theta, Phi) - DD_center;
					DDPP = DD;
					IFIX_D = (int)(DDPP);
					DDP = DDPP - IFIX_D;
					E_inc_grid = (1 - DDP) * E_inc_oblic[IFIX_D + disd] + DDP * E_inc_oblic[IFIX_D + 1 + disd];
					Ex_xyN_correction[ii - B_Tx_ini][jj - B_Ty_ini] = E_inc_grid * Vector_component_Exinc(Polar_psi, Phi, Theta);

					// Ex_inc (i+1/2, j, k1)
					DD = Calculate_Distance(ii + 0.5, jj, B_Tz_ter, Theta, Phi) - DD_center;
					DDPP = DD;
					IFIX_D = (int)(DDPP);
					DDP = DDPP - IFIX_D;
					E_inc_grid = (1 - DDP) * E_inc_oblic[IFIX_D + disd] + DDP * E_inc_oblic[IFIX_D + 1 + disd];
					Ex_xyP_correction[ii - B_Tx_ini][jj - B_Ty_ini] = E_inc_grid * Vector_component_Exinc(Polar_psi, Phi, Theta);

				}
				if (jj <= B_Ty_ter - 1)
				{
					// Ey_inc (i, j+1/2, k0)
					DD = Calculate_Distance(ii, jj + 0.5, B_Tz_ini, Theta, Phi) - DD_center;
					DDPP = DD;
					IFIX_D = (int)(DDPP);
					DDP = DDPP - IFIX_D;
					E_inc_grid = (1 - DDP) * E_inc_oblic[IFIX_D + disd] + DDP * E_inc_oblic[IFIX_D + 1 + disd];
					Ey_xyN_correction[ii - B_Tx_ini][jj - B_Ty_ini] = E_inc_grid * Vector_component_Eyinc(Polar_psi, Phi, Theta);

					// Ey_inc (i, j+1/2, k0)
					DD = Calculate_Distance(ii, jj + 0.5, B_Tz_ter, Theta, Phi) - DD_center;
					DDPP = DD;
					IFIX_D = (int)(DDPP);
					DDP = DDPP - IFIX_D;
					E_inc_grid = (1 - DDP) * E_inc_oblic[IFIX_D + disd] + DDP * E_inc_oblic[IFIX_D + 1 + disd];
					Ey_xyP_correction[ii - B_Tx_ini][jj - B_Ty_ini] = E_inc_grid * Vector_component_Eyinc(Polar_psi, Phi, Theta);
				}
			}
#pragma omp for schedule(dynamic)
		for (jj = B_Ty_ini; jj <= B_Ty_ter; jj++)
			for (kk = B_Tz_ini; kk <= B_Tz_ter; kk++)
			{
				if (jj <= B_Ty_ter - 1)
				{
					// Ex_inc (i0, j+1/2, k0)
					DD = Calculate_Distance(B_Tx_ini, jj + 0.5, kk, Theta, Phi) - DD_center;
					DDPP = DD;
					IFIX_D = (int)(DDPP);
					DDP = DDPP - IFIX_D;
					E_inc_grid = (1 - DDP) * E_inc_oblic[IFIX_D + disd] + DDP * E_inc_oblic[IFIX_D + 1 + disd];
					Ey_yzN_correction[jj - B_Ty_ini][kk - B_Tz_ini] = E_inc_grid * Vector_component_Eyinc(Polar_psi, Phi, Theta);

					// Ex_inc (10, j+1/2, k0)
					DD = Calculate_Distance(B_Tx_ter, jj + 0.5, kk, Theta, Phi) - DD_center;
					DDPP = DD;
					IFIX_D = (int)(DDPP);
					DDP = DDPP - IFIX_D;
					E_inc_grid = (1 - DDP) * E_inc_oblic[IFIX_D + disd] + DDP * E_inc_oblic[IFIX_D + 1 + disd];
					Ey_yzP_correction[jj - B_Ty_ini][kk - B_Tz_ini] = E_inc_grid * Vector_component_Eyinc(Polar_psi, Phi, Theta);
				}
				if (kk <= B_Tz_ter - 1)
				{
					// Ey_inc (i0, j, k+1/2)
					DD = Calculate_Distance(B_Tx_ini, jj, kk + 0.5, Theta, Phi) - DD_center;
					DDPP = DD;
					IFIX_D = (int)(DDPP);
					DDP = DDPP - IFIX_D;
					E_inc_grid = (1 - DDP) * E_inc_oblic[IFIX_D + disd] + DDP * E_inc_oblic[IFIX_D + 1 + disd];
					Ez_yzN_correction[jj - B_Ty_ini][kk - B_Tz_ini] = E_inc_grid * Vector_component_Ezinc(Polar_psi, Phi, Theta);

					// Ey_inc (i1, j, k+1/2)
					DD = Calculate_Distance(B_Tx_ter, jj, kk + 0.5, Theta, Phi) - DD_center;
					DDPP = DD;
					IFIX_D = (int)(DDPP);
					DDP = DDPP - IFIX_D;
					E_inc_grid = (1 - DDP) * E_inc_oblic[IFIX_D + disd] + DDP * E_inc_oblic[IFIX_D + 1 + disd];
					Ez_yzP_correction[jj - B_Ty_ini][kk - B_Tz_ini] = E_inc_grid * Vector_component_Ezinc(Polar_psi, Phi, Theta);
				}
			}
#pragma omp for schedule(dynamic)
		for (ii = B_Tx_ini; ii <= B_Tx_ter; ii++)
			for (kk = B_Tz_ini; kk <= B_Tz_ter; kk++)
			{
				if (ii <= B_Tx_ter - 1)
				{
					// Ex_inc (i+1/2, j0, k)
					DD = Calculate_Distance(ii + 0.5, B_Ty_ini, kk, Theta, Phi) - DD_center;
					DDPP = DD;
					IFIX_D = (int)(DDPP);
					DDP = DDPP - IFIX_D;
					E_inc_grid = (1 - DDP) * E_inc_oblic[IFIX_D + disd] + DDP * E_inc_oblic[IFIX_D + 1 + disd];
					Ex_xzN_correction[ii - B_Tx_ini][kk - B_Tz_ini] = E_inc_grid * Vector_component_Exinc(Polar_psi, Phi, Theta);

					// Ex_inc (i+1/2, j1, k)
					DD = Calculate_Distance(ii + 0.5, B_Ty_ter, kk, Theta, Phi) - DD_center;
					DDPP = DD;
					IFIX_D = (int)(DDPP);
					DDP = DDPP - IFIX_D;
					E_inc_grid = (1 - DDP) * E_inc_oblic[IFIX_D + disd] + DDP * E_inc_oblic[IFIX_D + 1 + disd];
					Ex_xzP_correction[ii - B_Tx_ini][kk - B_Tz_ini] = E_inc_grid * Vector_component_Exinc(Polar_psi, Phi, Theta);
				}
				if (kk <= B_Tz_ter - 1)
				{
					// Ex_inc (i, j0, k+1/2)
					DD = Calculate_Distance(ii, B_Ty_ini, kk + 0.5, Theta, Phi) - DD_center;
					DDPP = DD;
					IFIX_D = (int)(DDPP);
					DDP = DDPP - IFIX_D;
					E_inc_grid = (1 - DDP) * E_inc_oblic[IFIX_D + disd] + DDP * E_inc_oblic[IFIX_D + 1 + disd];
					Ez_xzN_correction[ii - B_Tx_ini][kk - B_Tz_ini] = E_inc_grid * Vector_component_Ezinc(Polar_psi, Phi, Theta);

					// Ex_inc (i, j1, k+1/2)
					DD = Calculate_Distance(ii, B_Ty_ter, kk + 0.5, Theta, Phi) - DD_center;
					DDPP = DD;
					IFIX_D = (int)(DDPP);
					DDP = DDPP - IFIX_D;
					E_inc_grid = (1 - DDP) * E_inc_oblic[IFIX_D + disd] + DDP * E_inc_oblic[IFIX_D + 1 + disd];
					Ez_xzP_correction[ii - B_Tx_ini][kk - B_Tz_ini] = E_inc_grid * Vector_component_Ezinc(Polar_psi, Phi, Theta);
				}
			}
	}
}

void TFDTD_Boundary_TFSF::Calcuate_inc_Field_H_OMP()
{
			double DD = 0;       // Distance
			double DDPP = 0;     // Distance from source corner
			int IFIX_D = 0;   // Gause Distance
			double DDP = 0;      // Distance from Gause distance
			double E_inc_grid;
			double H_inc_grid;
			int ii, jj, kk;
#pragma omp parallel private(ii, jj, kk, DD, DDPP, IFIX_D, DDP, E_inc_grid, H_inc_grid)
	{
#pragma omp for schedule(dynamic)
		for (jj = B_Ty_ini; jj <= B_Ty_ter; jj++)
			for (kk = B_Tz_ini; kk <= B_Tz_ter; kk++)
			{
				if (jj <= B_Ty_ter - 1)
				{
					// Hz_inc (i0-1/2, j+1/2, k)
					DD = Calculate_Distance(B_Tx_ini - 0.5, jj + 0.5, kk, Theta, Phi) - DD_center_H;
					DDPP = DD;// + 0.5;
					IFIX_D = (int)(DDPP);
					DDP = DDPP - IFIX_D;
					H_inc_grid = (1 - DDP)*H_inc_oblic[IFIX_D - 1 + disd] + DDP * H_inc_oblic[IFIX_D + disd];
					Hz_yzN_correction[jj - B_Ty_ini][kk - B_Tz_ini] = H_inc_grid * Vector_component_Hzinc(Polar_psi, Phi, Theta);

					// Hz_inc (i1-1/2, j+1/2, k)
					DD = Calculate_Distance((double)B_Tx_ter + 0.5, jj + 0.5, kk, Theta, Phi) - DD_center_H;
					DDPP = DD;// + 0.5;
					IFIX_D = (int)(DDPP);
					DDP = DDPP - IFIX_D;
					H_inc_grid = (1 - DDP)*H_inc_oblic[IFIX_D - 1 + disd] + DDP * H_inc_oblic[IFIX_D + disd];
					Hz_yzP_correction[jj - B_Ty_ini][kk - B_Tz_ini] = H_inc_grid * Vector_component_Hzinc(Polar_psi, Phi, Theta);
				}
				if (kk <= B_Tz_ter - 1)
				{
					// Hy_inc (i0-1/2, j, k+1/2)
					DD = Calculate_Distance((double)B_Tx_ini - 0.5, jj, kk + 0.5, Theta, Phi) - DD_center_H;
					DDPP = DD;// + 0.5;
					IFIX_D = (int)(DDPP);
					DDP = DDPP - IFIX_D;
					H_inc_grid = (1 - DDP)*H_inc_oblic[IFIX_D - 1 + disd] + DDP * H_inc_oblic[IFIX_D + disd];
					Hy_yzN_correction[jj - B_Ty_ini][kk - B_Tz_ini] = H_inc_grid * Vector_component_Hyinc(Polar_psi, Phi, Theta);

					// Hy_inc (i1-1/2, j, k+1/2)
					DD = Calculate_Distance((double)B_Tx_ter + 0.5, jj, kk + 0.5, Theta, Phi) - DD_center_H;
					DDPP = DD;// + 0.5;
					IFIX_D = (int)(DDPP);
					DDP = DDPP - IFIX_D;
					H_inc_grid = (1 - DDP)*H_inc_oblic[IFIX_D - 1 + disd] + DDP * H_inc_oblic[IFIX_D + disd];
					Hy_yzP_correction[jj - B_Ty_ini][kk - B_Tz_ini] = H_inc_grid * Vector_component_Hyinc(Polar_psi, Phi, Theta);
				}
			}
#pragma omp for schedule(dynamic)
		for (ii = B_Tx_ini; ii <= B_Tx_ter; ii++)
			for (kk = B_Tz_ini; kk <= B_Tz_ter; kk++)
			{
				if (kk <= B_Tz_ter - 1)
				{
					// Hx_inc (i, j0-1/2, k+1/2)
					DD = Calculate_Distance(ii, (double)B_Ty_ini - 0.5, kk + 0.5, Theta, Phi) - DD_center_H;
					DDPP = DD;// + 0.5;
					IFIX_D = (int)(DDPP);
					DDP = DDPP - IFIX_D;
					H_inc_grid = (1 - DDP)*H_inc_oblic[IFIX_D - 1 + disd] + DDP * H_inc_oblic[IFIX_D + disd];
					Hx_xzN_correction[ii - B_Tx_ini][kk - B_Tz_ini] = H_inc_grid * Vector_component_Hxinc(Polar_psi, Phi, Theta);

					// Hx_inc (i, j1-1/2, k+1/2)
					DD = Calculate_Distance(ii, (double)B_Ty_ter + 0.5, kk + 0.5, Theta, Phi) - DD_center_H;
					DDPP = DD;// + 0.5;
					IFIX_D = (int)(DDPP);
					DDP = DDPP - IFIX_D;
					H_inc_grid = (1 - DDP)*H_inc_oblic[IFIX_D - 1 + disd] + DDP * H_inc_oblic[IFIX_D + disd];
					Hx_xzP_correction[ii - B_Tx_ini][kk - B_Tz_ini] = H_inc_grid * Vector_component_Hxinc(Polar_psi, Phi, Theta);
				}
				if (ii <= B_Tx_ter - 1)
				{
					// Hz_inc (i+1/2, j0-1/2, k)
					DD = Calculate_Distance(ii + 0.5, (double)B_Ty_ini - 0.5, kk, Theta, Phi) - DD_center_H;
					DDPP = DD;// + 0.5;
					IFIX_D = (int)(DDPP);
					DDP = DDPP - IFIX_D;
					H_inc_grid = (1 - DDP)*H_inc_oblic[IFIX_D - 1 + disd] + DDP * H_inc_oblic[IFIX_D + disd];
					Hz_xzN_correction[ii - B_Tx_ini][kk - B_Tz_ini] = H_inc_grid * Vector_component_Hzinc(Polar_psi, Phi, Theta);

					// Hz_inc (i+1/2, j1-1/2, k)
					DD = Calculate_Distance(ii + 0.5, (double)B_Ty_ter + 0.5, kk, Theta, Phi) - DD_center_H;
					DDPP = DD;// + 0.5;
					IFIX_D = (int)(DDPP);
					DDP = DDPP - IFIX_D;
					H_inc_grid = (1 - DDP)*H_inc_oblic[IFIX_D - 1 + disd] + DDP * H_inc_oblic[IFIX_D + disd];
					Hz_xzP_correction[ii - B_Tx_ini][kk - B_Tz_ini] = H_inc_grid * Vector_component_Hzinc(Polar_psi, Phi, Theta);
				}
			}
#pragma omp for schedule(dynamic)
		for (ii = B_Tx_ini; ii <= B_Tx_ter; ii++)
			for (jj = B_Ty_ini; jj <= B_Ty_ter; jj++)
			{
				if (jj <= B_Ty_ter - 1)
				{
					// Hx_inc (i, j+1/2, k0-1/2)
					DD = Calculate_Distance(ii, jj + 0.5, (double)B_Tz_ini - 0.5, Theta, Phi) - DD_center_H;
					DDPP = DD;// + 0.5;
					IFIX_D = (int)(DDPP);
					DDP = DDPP - IFIX_D;
					H_inc_grid = (1 - DDP)*H_inc_oblic[IFIX_D - 1 + disd] + DDP * H_inc_oblic[IFIX_D + disd];
					Hx_xyN_correction[ii - B_Tx_ini][jj - B_Ty_ini] = H_inc_grid * Vector_component_Hxinc(Polar_psi, Phi, Theta);

					// Hx_inc (i, j+1/2, k1-1/2)
					DD = Calculate_Distance(ii, jj + 0.5, (double)B_Tz_ter + 0.5, Theta, Phi) - DD_center_H;
					DDPP = DD;// + 0.5;
					IFIX_D = (int)(DDPP);
					DDP = DDPP - IFIX_D;
					H_inc_grid = (1 - DDP)*H_inc_oblic[IFIX_D - 1 + disd] + DDP * H_inc_oblic[IFIX_D + disd];
					Hx_xyP_correction[ii - B_Tx_ini][jj - B_Ty_ini] = H_inc_grid * Vector_component_Hxinc(Polar_psi, Phi, Theta);
				}
				if (ii <= B_Tx_ter - 1)
				{
					// Hy_inc  (i+1/2, j, k0-1/2)
					DD = Calculate_Distance(ii + 0.5, jj, (double)B_Tz_ini - 0.5, Theta, Phi) - DD_center_H;
					DDPP = DD;// + 0.5;
					IFIX_D = (int)(DDPP);
					DDP = DDPP - IFIX_D;
					H_inc_grid = (1 - DDP)*H_inc_oblic[IFIX_D - 1 + disd] + DDP * H_inc_oblic[IFIX_D + disd];
					Hy_xyN_correction[ii - B_Tx_ini][jj - B_Ty_ini] = H_inc_grid * Vector_component_Hyinc(Polar_psi, Phi, Theta);

					// Hy_inc  (i+1/2, j, k1-1/2)
					DD = Calculate_Distance(ii + 0.5, jj, (double)B_Tz_ter + 0.5, Theta, Phi) - DD_center_H;
					DDPP = DD;// + 0.5;
					IFIX_D = (int)(DDPP);
					DDP = DDPP - IFIX_D;
					H_inc_grid = (1 - DDP)*H_inc_oblic[IFIX_D - 1 + disd] + DDP * H_inc_oblic[IFIX_D + disd];
					Hy_xyP_correction[ii - B_Tx_ini][jj - B_Ty_ini] = H_inc_grid * Vector_component_Hyinc(Polar_psi, Phi, Theta);
				}
			}
	}
}

void TFDTD_Boundary_TFSF::Correction_X_Boundary_Exn_OMP()
{
	int jj, kk;
#pragma omp parallel private(jj, kk)
	{
#pragma omp for schedule(dynamic)
		// x- correction
		for (jj = B_Ty_ini; jj < B_Ty_ter; jj++)
			for (kk = B_Tz_ini; kk <= B_Tz_ter; kk++)
				// Ey( i = i0; j = j0 + 1/2, ... , j1 - 1/2; k = k0, ... k1)
				// Ey|(n+1)(i0,j,k) = ~ + dt/eps0/dx * Hz_xyN_corr->(i0-1, j, k)
				Ey[B_Tx_ini][jj][kk] = Ey[B_Tx_ini][jj][kk] + Cedx1 * Hz_yzN_correction[jj - B_Ty_ini][kk - B_Tz_ini];
		//*(*(*(*EPy + B_Tx_ini) + jj) + kk) = *(*(*(*EPy + B_Tx_ini) + jj) + kk) + dt / eps_0 / dx * Hz_yzN_correction[jj - B_Ty_ini][kk - B_Tz_ini];

#pragma omp for schedule(dynamic)
		for (jj = B_Ty_ini; jj <= B_Ty_ter; jj++)
			for (kk = B_Tz_ini; kk < B_Tz_ter; kk++)
				// Ez( i = i0; j = j0, ... , j1; k = k0+1/2, ... k1-1/2)
				// Ez|(n+1)(i0, j = j0, ... , j1, k =  k0 + 1/2, ... , k1 - 1/2
				Ez[B_Tx_ini][jj][kk] = Ez[B_Tx_ini][jj][kk] - Cedx1 * Hy_yzN_correction[jj - B_Ty_ini][kk - B_Tz_ini];
		//*(*(*(*EPz + B_Tx_ini) + jj) + kk) = *(*(*(*EPz + B_Tx_ini) + jj) + kk) - dt / eps_0 / dx * Hy_yzN_correction[jj - B_Ty_ini][kk - B_Tz_ini];
	}
}
void TFDTD_Boundary_TFSF::Correction_X_Boundary_Exp_OMP()
{
	int jj, kk;
#pragma omp parallel private(jj, kk)
	{
#pragma omp for schedule(dynamic)
		// x+ correction
		for (jj = B_Ty_ini; jj < B_Ty_ter; jj++)
			for (kk = B_Tz_ini; kk <= B_Tz_ter; kk++)
				// Ey( i = i1; j = j0 + 1/2, ... , j1 - 1/2; k = k0, ..., k1)
				// Ey|(n+1)(i1,j,k) = ~ + dt/eps0/dx * Hz_xyN_corr->(i0-1, j, k)
				Ey[B_Tx_ter][jj][kk] = Ey[B_Tx_ter][jj][kk] - Cedx2 * Hz_yzP_correction[jj - B_Ty_ini][kk - B_Tz_ini];
		//*(*(*(*EPy + B_Tx_ter) + jj) + kk) = *(*(*(*EPy + B_Tx_ter) + jj) + kk) - (dt / eps_0 / dx) * Hz_yzP_correction[jj - B_Ty_ini][kk - B_Tz_ini];
#pragma omp for schedule(dynamic)
		// x+ correction
		for (jj = B_Ty_ini; jj <= B_Ty_ter; jj++)
			for (kk = B_Tz_ini; kk < B_Tz_ter; kk++)
				// Ez( i = i0; j = j0, ... , j1; k = k0+1/2, ... k1-1/2)
				// Ez|(n+1)(i1, j = j0, ... , j1, k =  k0 + 1/2, ... , k1 - 1/2
				Ez[B_Tx_ter][jj][kk] = Ez[B_Tx_ter][jj][kk] + Cedx2 * Hy_yzP_correction[jj - B_Ty_ini][kk - B_Tz_ini];
		//*(*(*(*EPz + B_Tx_ter) + jj) + kk) = *(*(*(*EPz + B_Tx_ter) + jj) + kk) + (dt / eps_0 / dx) * Hy_yzP_correction[jj - B_Ty_ini][kk - B_Tz_ini];
	}
}
void TFDTD_Boundary_TFSF::Correction_Y_Boundary_Eyn_OMP()
{
	int ii, kk;
#pragma omp parallel private(ii, kk)
	{
#pragma omp for schedule(dynamic)
		// y- correction
		for (ii = B_Tx_ini; ii < B_Tx_ter; ii++)
			for (kk = B_Tz_ini; kk <= B_Tz_ter; kk++)
			{
				// Ex( i = i0 ... i1; j = j0; k = k0+1/2, ... k1-1/2)
				// Ex|(n+1)(i0 + 1/2, ... , i1-1/2, j = j0, k =  k0, ... , k1)
				//*(*(*(*EPx + ii) + B_Ty_ini) + kk) = *(*(*(*EPx + ii) + B_Ty_ini) + kk) - dt / eps_0 / dy * Hz_xzN_correction[ii - B_Tx_ini][kk - B_Tz_ini];
				Ex[ii][B_Ty_ini][kk] = Ex[ii][B_Ty_ini][kk] - Cedy1 * Hz_xzN_correction[ii - B_Tx_ini][kk - B_Tz_ini];
			}
#pragma omp for schedule(dynamic)
		for (ii = B_Tx_ini; ii <= B_Tx_ter; ii++)
			for (kk = B_Tz_ini; kk < B_Tz_ter; kk++)
			{
				// Ez(n+1)(i,j0,k) = ~ + dt/eps0/dy * Hx_zxN_corr->(i0-1, j, k)
				//*(*(*(*EPz + ii) + B_Ty_ini) + kk) = *(*(*(*EPz + ii) + B_Ty_ini) + kk) + dt / eps_0 / dy * Hx_xzN_correction[ii - B_Tx_ini][kk - B_Tz_ini];
				Ez[ii][B_Ty_ini][kk] = Ez[ii][B_Ty_ini][kk] + Cedy1 * Hx_xzN_correction[ii - B_Tx_ini][kk - B_Tz_ini];
			}
	}
}
void TFDTD_Boundary_TFSF::Correction_Y_Boundary_Eyp_OMP()
{
	int ii, kk;
#pragma omp parallel private(ii, kk)
	{
#pragma omp for schedule(dynamic)
		// y+ correction
		for (ii = B_Tx_ini; ii < B_Tx_ter; ii++)
			for (kk = B_Tz_ini; kk <= B_Tz_ter; kk++)
			{
				// Ex( i = i0 ... i1; j = j0; k = k0+1/2, ... k1-1/2)
				// Ex|(n+1)(i0 + 1/2, ... , i1-1/2, j = j0, k =  k0, ... , k1)
				//*(*(*(*EPx + ii) + B_Ty_ter) + kk) = *(*(*(*EPx + ii) + B_Ty_ter) + kk) + dt / eps_0 / dy * Hz_xzP_correction[ii - B_Tx_ini][kk - B_Tz_ini];
				Ex[ii][B_Ty_ter][kk] = Ex[ii][B_Ty_ter][kk] + Cedy2 * Hz_xzP_correction[ii - B_Tx_ini][kk - B_Tz_ini];
			}
#pragma omp for schedule(dynamic)
		for (ii = B_Tx_ini; ii <= B_Tx_ter; ii++)
			for (kk = B_Tz_ini; kk < B_Tz_ter; kk++)
			{
				// Ez(n+1)(i,j0,k) = ~ + dt/eps0/dy * Hx_zxN_corr->(i0-1, j, k)
				//*(*(*(*EPz + ii) + B_Ty_ter) + kk) = *(*(*(*EPz + ii) + B_Ty_ter) + kk) - dt / eps_0 / dy * Hx_xzP_correction[ii - B_Tx_ini][kk - B_Tz_ini];
				Ez[ii][B_Ty_ter][kk] = Ez[ii][B_Ty_ter][kk] - Cedy2 * Hx_xzP_correction[ii - B_Tx_ini][kk - B_Tz_ini];
			}
	}
}
void TFDTD_Boundary_TFSF::Correction_Z_Boundary_Ezn_OMP()
{
	int ii, jj;
#pragma omp parallel private(ii, jj)
	{
#pragma omp for schedule(dynamic)
		// z- correction
		for (ii = B_Tx_ini; ii < B_Tx_ter; ii++)
			for (jj = B_Ty_ini; jj <= B_Ty_ter; jj++)
			{
				// Ex( i = i0 +1/2 ... i1 -1/2; j = j0 ... j1; k = k0)
				// Ex|(n+1)(i0 + 1/2, ... , i1-1/2, j = j0, k =  k0) = ~ - dt/eps0/dy * Hy_zxN_corr->(i, j, k0-1/2)
				//*(*(*(*EPx + ii) + jj) + B_Tz_ini) = *(*(*(*EPx + ii) + jj) + B_Tz_ini) + dt / eps_0 / dz * Hy_xyN_correction[ii - B_Tx_ini][jj - B_Ty_ini];
				Ex[ii][jj][B_Tz_ini] = Ex[ii][jj][B_Tz_ini] + Cedz1 * Hy_xyN_correction[ii - B_Tx_ini][jj - B_Ty_ini];
			}
#pragma omp for schedule(dynamic)
		for (ii = B_Tx_ini; ii <= B_Tx_ter; ii++)
			for (jj = B_Ty_ini; jj < B_Ty_ter; jj++)
			{
				// Ez(n+1)(i0, ... , i1; j= j0+1/2, ... , j1 - 1/2, k=k0) = ~ + dt/eps0/dy * Hx_zxN_corr->(i, j, k0-1/2)
				//*(*(*(*EPy + ii) + jj) + B_Tz_ini) = *(*(*(*EPy + ii) + jj) + B_Tz_ini) - dt / eps_0 / dz * Hx_xyN_correction[ii - B_Tx_ini][jj - B_Ty_ini];
				Ey[ii][jj][B_Tz_ini] = Ey[ii][jj][B_Tz_ini] - Cedz1 * Hx_xyN_correction[ii - B_Tx_ini][jj - B_Ty_ini];
			}
	}
}
void TFDTD_Boundary_TFSF::Correction_Z_Boundary_Ezp_OMP()
{
	int ii, jj;
#pragma omp parallel private(ii, jj)
	{
#pragma omp for schedule(dynamic)
		// z+ correction
		for (ii = B_Tx_ini; ii < B_Tx_ter; ii++)
			for (jj = B_Ty_ini; jj <= B_Ty_ter; jj++)
			{
				// Ex( i = i0 +1/2 ... i1 -1/2; j = j0 ... j1; k = k0)
				// Ex|(n+1)(i0 + 1/2, ... , i1-1/2, j = j0, k =  k0) = ~ - dt/eps0/dy * Hy_zxN_corr->(i, j, k0-1/2)
				//*(*(*(*EPx + ii) + jj) + B_Tz_ter) = *(*(*(*EPx + ii) + jj) + B_Tz_ter) - dt / eps_0 / dz * Hy_xyP_correction[ii - B_Tx_ini][jj - B_Ty_ini];
				Ex[ii][jj][B_Tz_ter] = Ex[ii][jj][B_Tz_ter] - Cedz2 * Hy_xyP_correction[ii - B_Tx_ini][jj - B_Ty_ini];
			}
#pragma omp for schedule(dynamic)
		for (ii = B_Tx_ini; ii <= B_Tx_ter; ii++)
			for (jj = B_Ty_ini; jj < B_Ty_ter; jj++)
			{
				// Ez(n+1)(i0, ... , i1; j= j0+1/2, ... , j1 - 1/2, k=k0) = ~ + dt/eps0/dy * Hx_zxN_corr->(i, j, k0-1/2)
				//*(*(*(*EPy + ii) + jj) + B_Tz_ter) = *(*(*(*EPy + ii) + jj) + B_Tz_ter) + dt / eps_0 / dz * Hx_xyP_correction[ii - B_Tx_ini][jj - B_Ty_ini];
				Ey[ii][jj][B_Tz_ter] = Ey[ii][jj][B_Tz_ter] + Cedz2 * Hx_xyP_correction[ii - B_Tx_ini][jj - B_Ty_ini];
			}
		// cout << Dy[20][20][B_Tz_ter] << endl;
	}
}
void TFDTD_Boundary_TFSF::Correction_X_Boundary_Hxn_OMP()
{
	int kk, jj;
#pragma omp parallel private(kk, jj)
	{
#pragma omp for schedule(dynamic)
		// x- correction i0 -1/2 face
		for (jj = B_Ty_ini; jj <= B_Ty_ter; jj++)
			for (kk = B_Tz_ini; kk <= B_Tz_ter; kk++)
			{
				if (jj < B_Ty_ter) // jj = j0+1/2 .. j1 - 1/2
				{
					// Hz(i = i0-1/2 ; j= j0 +1/2, ... , j1 -1/2 ; k =  k0, ... , k1 )
					// Hz|(n+1/2)(i0 -1/2,j,k) = ~ + dt/mu0/dx * Ey_xyN_corr->(i0, j, k)
					Hz[B_Tx_ini - 1][jj][kk] = Hz[B_Tx_ini - 1][jj][kk] + Chdx1 * Ey_yzN_correction[jj - B_Ty_ini][kk - B_Tz_ini];
					//				*(*(*(*HPz + B_Tx_ini) + jj) + kk) = *(*(*(*HPz + B_Tx_ini) + jj) + kk) + dt / mu_0 / dx * Ey_yzN_correction[jj - B_Ty_ini][kk - B_Tz_ini];
				}

				if (kk < B_Tz_ter) // kk = k0+1/2 ... k1-1/2
				{
					// Hy(i = i0-1/2 ; j= j0, ... , j1, k =  k0 + 1/2, ... , k1 - 1/2 )
					// Hy|(n+1/2)(i0 -1/2, j, k) = ~ -dt/mu0/dx * Ez_xyN_corr->(i0, j, k)
					Hy[B_Tx_ini - 1][jj][kk] = Hy[B_Tx_ini - 1][jj][kk] - Chdx1 * Ez_yzN_correction[jj - B_Ty_ini][kk - B_Tz_ini];
					//				*(*(*(*HPy + B_Tx_ini) + jj) + kk) = *(*(*(*HPy + B_Tx_ini) + jj) + kk) - dt / mu_0 / dx * Ez_yzN_correction[jj - B_Ty_ini][kk - B_Tz_ini];
				}
			}
	}
}
void TFDTD_Boundary_TFSF::Correction_X_Boundary_Hxp_OMP()
{
	int kk, jj;
#pragma omp parallel private(kk, jj)
	{
#pragma omp for schedule(dynamic)

		// x+ correction
		for (jj = B_Ty_ini; jj <= B_Ty_ter; jj++)
			for (kk = B_Tz_ini; kk <= B_Tz_ter; kk++)
			{
				if (jj < B_Ty_ter) // jj = j0+1/2 .. j1 - 1/2
				{
					// Hz(i = i1+1/2 ; j= j0 +1/2, ... , j1 -1/2 ; k =  k0, ... , k1 )
					// Hz|(n+1/2)(i1 +1/2,j,k) = ~ - dt/mu0/dx * Ey_xyP_corr->(i0, j, k)
					Hz[B_Tx_ter][jj][kk] = Hz[B_Tx_ter][jj][kk] - Chdx2 * Ey_yzP_correction[jj - B_Ty_ini][kk - B_Tz_ini];
					//				*(*(*(*HPz + B_Tx_ter-1) + jj) + kk) = *(*(*(*HPz + B_Tx_ter-1) + jj) + kk) - dt / mu_0 / dx * Ey_yzP_correction[jj - B_Ty_ini][kk - B_Tz_ini];
				}

				if (kk < B_Tz_ter) // kk = k0+1/2 ... k1-1/2
				{
					// Hy(i = i1+1/2 ; j= j0, ... , j1, k =  k0 + 1/2, ... , k1 - 1/2 )
					// Hy|(n+1/2)(i1 +1/2, j, k) = ~ -dt/mu0/dx * Ez_xyN_corr->(i0, j, k)
					Hy[B_Tx_ter][jj][kk] = Hy[B_Tx_ter][jj][kk] + Chdx2 * Ez_yzP_correction[jj - B_Ty_ini][kk - B_Tz_ini];
					//				*(*(*(*HPy + B_Tx_ter-1) + jj) + kk) = *(*(*(*HPy + B_Tx_ter-1) + jj) + kk) + dt / mu_0 / dx * Ez_yzP_correction[jj - B_Ty_ini][kk - B_Tz_ini];
				}
			}
	}
}
void TFDTD_Boundary_TFSF::Correction_Y_Boundary_Hyn_OMP()
{
	int kk, ii;
#pragma omp parallel private(kk, ii)
	{
#pragma omp for schedule(dynamic)

		// y- correction : j = j0 - 1/2
		for (ii = B_Tx_ini; ii <= B_Tx_ter; ii++)
			for (kk = B_Tz_ini; kk <= B_Tz_ter; kk++)
			{
				if (ii < B_Tx_ter)
				{
					// Hz( i = i0+1/2 ... i1 -1/2; j = j0- 1/2; k = k0, ... k1)
					// Hz|(n+1/2)(i; j = j0-1/2 ; k) = ~ - dt/mu0/dy * Ex
					Hz[ii][B_Ty_ini - 1][kk] = Hz[ii][B_Ty_ini - 1][kk] - Chdy1 * Ex_xzN_correction[ii - B_Tx_ini][kk - B_Tz_ini];
					//				*(*(*(*HPz + ii) + B_Ty_ini) + kk) = *(*(*(*HPz + ii) + B_Ty_ini) + kk) - dt / mu_0 / dy * Ex_xzN_correction[ii - B_Tx_ini][kk - B_Tz_ini];
				}
				if (kk < B_Tz_ter)
				{
					// Hx( i = i0, ... , i1; j = j0-1/2 ; k = k0+1/2, ..., k1 - 1/2)
					//			// Hx(n+1.2)(i,j0,k) = ~ + dt/eps0/dy * Hx_zxN_corr->(i0-1, j, k)
					Hx[ii][B_Ty_ini - 1][kk] = Hx[ii][B_Ty_ini - 1][kk] + Chdy1 * Ez_xzN_correction[ii - B_Tx_ini][kk - B_Tz_ini];
					//				*(*(*(*HPx + ii) + B_Ty_ini) + kk) = *(*(*(*HPx + ii) + B_Ty_ini) + kk) + dt / mu_0 / dy * Ez_xzN_correction[ii - B_Tx_ini][kk - B_Tz_ini];
				}
			}
	}
}
void TFDTD_Boundary_TFSF::Correction_Y_Boundary_Hyp_OMP()
{
	int kk, ii;
#pragma omp parallel private(kk, ii)
	{
#pragma omp for schedule(dynamic)

		// y+ correction
		for (ii = B_Tx_ini; ii <= B_Tx_ter; ii++)
			for (kk = B_Tz_ini; kk <= B_Tz_ter; kk++)
			{
				if (ii < B_Tx_ter)
				{
					// Hz( i = i0+1/2, ..., i1 -1/2 ; j = j1+ 1/2 ; k = k0, ... k1)
					// Hz|(n+1/2)(i; j = j1-1/2 ; k) = ~ + dt/mu0/dy * Ex
					Hz[ii][B_Ty_ter][kk] = Hz[ii][B_Ty_ter][kk] + Chdy2 * Ex_xzP_correction[ii - B_Tx_ini][kk - B_Tz_ini];
					//				*(*(*(*HPz + ii) + B_Ty_ter-1) + kk) = *(*(*(*HPz + ii) + B_Ty_ter-1) + kk) + dt / mu_0 / dy * Ex_xzP_correction[ii - B_Tx_ini][kk - B_Tz_ini];
				}
				if (kk < B_Tz_ter)
				{
					// Hx( i = i0, ... , i1 ; j = j1 + 1/2 ; k = k0+1/2, ..., k1 - 1/2)
					//			// Hx(n+1.2)(i,j1 - 1/2,k) = ~ - dt/eps0/dy * Hx_zxN_corr->(i0-1, j, k)
					Hx[ii][B_Ty_ter][kk] = Hx[ii][B_Ty_ter][kk] - Chdy2 * Ez_xzP_correction[ii - B_Tx_ini][kk - B_Tz_ini];
					//				*(*(*(*HPx + ii) + B_Ty_ter-1) + kk) = *(*(*(*HPx + ii) + B_Ty_ter-1) + kk) - dt / mu_0 / dy * Ez_xzP_correction[ii - B_Tx_ini][kk - B_Tz_ini];
				}
			}
	}
}
void TFDTD_Boundary_TFSF::Correction_Z_Boundary_Hzn_OMP()
{
	int jj, ii;
#pragma omp parallel private(jj, ii)
	{
#pragma omp for schedule(dynamic)
		// z- correction
		for (ii = B_Tx_ini; ii <= B_Tx_ter; ii++)
			for (jj = B_Ty_ini; jj <= B_Ty_ter; jj++)
			{
				if (ii < B_Tx_ter)
				{
					// Hy( i = i0 +1/2 ... i1 -1/2; j = j0, ..., j1; k = k0 - 1/2)
					// Hy|(n+1/2)(i, j, k0 -1/2) =  ~ + dt/mu0/dz * Ex
					Hy[ii][jj][B_Tz_ini - 1] = Hy[ii][jj][B_Tz_ini - 1] + Chdz1 * Ex_xyN_correction[ii - B_Tx_ini][jj - B_Ty_ini];
					//				*(*(*(*HPy + ii) + jj) + B_Tz_ini) = *(*(*(*HPy + ii) + jj) + B_Tz_ini) + dt / mu_0 / dz * Ex_xyN_correction[ii - B_Tx_ini][jj - B_Ty_ini];
				}
				if (jj < B_Ty_ter)
				{
					// Hx( i = i0, ..., i1; j = j0 + 1/2, ..., j1 -1/2; k = k0 - 1/2)
					// Hx(n+1)(i,j0,k) = ~ - dt/eps0/dy * Ey_zxN_corr->(i0-1, j, k)
					Hx[ii][jj][B_Tz_ini - 1] = Hx[ii][jj][B_Tz_ini - 1] - Chdz1 * Ey_xyN_correction[ii - B_Tx_ini][jj - B_Ty_ini];
					//				*(*(*(*HPx + ii) + jj) + B_Tz_ini) = *(*(*(*HPx + ii) + jj) + B_Tz_ini) - dt / mu_0 / dz * Ey_xyN_correction[ii - B_Tx_ini][jj - B_Ty_ini];
				}
			}
	}
	// debug code
	//	int ci = (B_Tx_ini + B_Tx_ter) / 2;
	//	int cj = (B_Ty_ini + B_Ty_ter) / 2;
	//	pf_debug2 = fopen("./Inc_Apply_to_Hn.txt", "a");
	//	fprintf(pf_debug2, "%le\t%le\n", dt / mu_0 / dz *Ex_xyN_correction[ci][cj], dt / mu_0 / dz *Ey_xyN_correction[ci][cj]);
	//	fclose(pf_debug2);
}
void TFDTD_Boundary_TFSF::Correction_Z_Boundary_Hzp_OMP()
{
	int jj, ii;
#pragma omp parallel private(jj, ii)
	{
#pragma omp for schedule(dynamic)
		// z+ correction
		for (ii = B_Tx_ini; ii <= B_Tx_ter; ii++)
			for (jj = B_Ty_ini; jj <= B_Ty_ter; jj++)
			{
				if (ii < B_Tx_ter)
				{
					// Hy( i = i0 +1/2 ... i1 -1/2; j = j0, ..., j1; k = k0 - 1/2)
					// Hy|(n+1/2)(i, j, k0 -1/2) =  ~ + dt/mu0/dz * Ex
					Hy[ii][jj][B_Tz_ter] = Hy[ii][jj][B_Tz_ter] - Chdz2 * Ex_xyP_correction[ii - B_Tx_ini][jj - B_Ty_ini];
					//				*(*(*(*HPy + ii) + jj) + B_Tz_ter-1) = *(*(*(*HPy + ii) + jj) + B_Tz_ter-1) - dt / mu_0 / dz * Ex_xyP_correction[ii - B_Tx_ini][jj - B_Ty_ini];
				}
				if (jj < B_Ty_ter)
				{
					// Hx( i = i0, ..., i1; j = j0 + 1/2, ..., j1 -1/2; k = k0 - 1/2)
					// Hx(n+1)(i,j0,k) = ~ - dt/eps0/dy * Ey_zxN_corr->(i0-1, j, k)
					Hx[ii][jj][B_Tz_ter] = Hx[ii][jj][B_Tz_ter] + Chdz2 * Ey_xyP_correction[ii - B_Tx_ini][jj - B_Ty_ini];
					//				*(*(*(*HPx + ii) + jj) + B_Tz_ter-1) = *(*(*(*HPx + ii) + jj) + B_Tz_ter-1) + dt / mu_0 / dz * Ey_xyP_correction[ii - B_Tx_ini][jj - B_Ty_ini];
				}
			}
	}
}

void TFDTD_Boundary_TFSF::DApply_OMP()
{
	if (B_Tx_ini != 0)
		Correction_X_Boundary_Exn_OMP();
	if (B_Tx_ter != txsize)
		Correction_X_Boundary_Exp_OMP();
	if (B_Ty_ini != 0)
		Correction_Y_Boundary_Eyn_OMP();
	if (B_Ty_ter != tysize)
		Correction_Y_Boundary_Eyp_OMP();
	if (B_Tz_ini != 0)
		Correction_Z_Boundary_Ezn_OMP();
	if (B_Tz_ter != tzsize)
		Correction_Z_Boundary_Ezp_OMP();
}

void TFDTD_Boundary_TFSF::HApply_OMP()
{
	if (B_Tx_ini != 0)
		Correction_X_Boundary_Hxn_OMP();
	if (B_Tx_ter != txsize)
		Correction_X_Boundary_Hxp_OMP();
	if (B_Ty_ini != 0)
		Correction_Y_Boundary_Hyn_OMP();
	if (B_Ty_ter != tysize)
		Correction_Y_Boundary_Hyp_OMP();
	if (B_Tz_ini != 0)
		Correction_Z_Boundary_Hzn_OMP();
	if (B_Tz_ter != tzsize)
		Correction_Z_Boundary_Hzp_OMP();
}

void TFDTD_Boundary_TFSF::TFSF_Electric_Field_OMP(double Pulse)
{
	this->Pulse = Pulse;
	Electric_Field_update();
	Calcuate_inc_Field_E_OMP();
}

void TFDTD_Boundary_TFSF::TFSF_Magnetic_Field_OMP()
{
	Magnetic_Field_update();
	Calcuate_inc_Field_H_OMP();
}