//---------------------------------------------------------------------------
#ifndef TFDTD_BOUNDARY_TFSFH
#define TFDTD_BOUNDARY_TFSFH
//---------------------------------------------------------------------------
#define _USE_MATH_DEFINES
#include <iostream>
#include <math.h>
#include "TFunctions.h"
#include <stdio.h>
#include <stdlib.h>

class TFDTD_Boundary_TFSF
{
public:
	__fastcall TFDTD_Boundary_TFSF() {};
	 TFDTD_Boundary_TFSF(struct_Field_carrier *pFc);
	__fastcall TFDTD_Boundary_TFSF(double ***Ex, double ***Ey, double ***Ez, double ***Hx, double ***Hy, double ***Hz);
	__fastcall ~TFDTD_Boundary_TFSF();
protected:
	double m_pi;
public:
	//void insert_parameters(int txsize/*txsize*/, int tysize /*tysize*/, int tzsize /*tzsize*/, double dx/*dx*/, double dy /*dy*/, double dz/*dz*/, double dt/*dt*/, double f0/*f0*/, int Nstep /*Nstep*/, int xx00/*x0plane*/, int xx11/*x0plane*/, int yy00/*x0plane*/, int yy11/*x0plane*/, int zz00/*x0plane*/, int zz11/*x0plane*/, double input_theta, double input_phi, double polar_psi);
	//void insert_parameters(int txsize/*txsize*/, int tysize /*tysize*/, int tzsize/*tzsize*/, double dx/*dx*/, double dy /*dy*/, double dz/*dz*/, double dt/*dt*/, int f0/*f0*/, int Nstep /*Nstep*/, int xx00/*x0plane*/, int xx11/*x0plane*/, int yy00/*x0plane*/, int yy11/*x0plane*/, int zz00/*x0plane*/, int zz11/*x0plane*/, double input_phi, double input_theta, double polar_psi);
	void insert_parameters(int txsize, int tysize, int tzsize, double dx, double dy, double dz, double dt, double f0, int Nstep);
	void insert_inc_angle(double input_Theta, double input_phi, double polar_psi);
	void insert_TFSF_Box(int xx00, int xx11, int yy00, int yy11, int zz00, int zz11);	
	void Set_Pulse(double Pulse);
	void Get_Polarization(double polar_psi);
	void DApply();
	void HApply();

protected: // Correction fuctions
	void Correction_X_Boundary_Exn();
	void Correction_X_Boundary_Exp();
	void Correction_Y_Boundary_Eyn();
	void Correction_Y_Boundary_Eyp();
	void Correction_Z_Boundary_Ezn();
	void Correction_Z_Boundary_Ezp();

	void Correction_X_Boundary_Hxn();
	void Correction_X_Boundary_Hxp();
	void Correction_Y_Boundary_Hyn();
	void Correction_Y_Boundary_Hyp();
	void Correction_Z_Boundary_Hzn();
	void Correction_Z_Boundary_Hzp();

	void Initial_condition();
private:
	//basic functions;
	virtual void Calcuate_inc_Field();
	void Calcuate_Initial_Field();

	//	void Initial_condition(int Nststep);
	void Define_Initial_Position();
	void Input_Boundary(int xi, int xt, int yi, int yt, int zi, int zt);
	void Input_Degree(double input_phi, double input_theta, double polar_psi);
	void Electric_Field_update();
	void Magnetic_Field_update();
	// velocity of direction of planewave
	double v_nomal, v_oblic, v_ratio;
	double* E_inc_oblic;
	double* H_inc_oblic;
	// Mur
	double abs_Gamma_media;
	double* E_inc_oblic_mur;

	// Boundary of Planewave
	double input_phi, input_theta;
	double Polar_psi;
	int disd;

	// Direction of Planewave (rad)
	double kk_inc, Theta, Phi;
	double DD_center;
	double DD_center_H;

	// 1-D FDTD spatial size
	int Nststep;

	// Pulse amp
	double Pulse;

protected:
	// basic parameters
	int txsize, tysize, tzsize;
	double *dx, *dy, *dz, dt, f0;
	double dds;
	double eps_0;
	double mu_0;
	double c_0;
	double Cedx1, Cedy1, Cedz1;
	double Cedx2, Cedy2, Cedz2;

	double Chdx1, Chdy1, Chdz1;
	double Chdx2, Chdy2, Chdz2;

	// Boundary of Planewave
	int B_Tx_ini, B_Tx_ter;
	int B_Ty_ini, B_Ty_ter;
	int B_Tz_ini, B_Tz_ter;
	int inc_plane_x, inc_plane_y, inc_plane_z;

protected:
	struct_Field_carrier* pFc;
	void Define_dds(int xpos, int ypos, int zpos);
	// Field Pointer         alloc X
	double ***Ex, ***Ey, ***Ez;
	double ***Hx, ***Hy, ***Hz;
	// Field correction variables
	double **Ex_xyN_correction;
	double **Ex_xyP_correction;
	double **Ex_xzN_correction;
	double **Ex_xzP_correction;

	double **Ey_xyN_correction;
	double **Ey_xyP_correction;
	double **Ey_yzN_correction;
	double **Ey_yzP_correction;

	double **Ez_xzN_correction;
	double **Ez_xzP_correction;
	double **Ez_yzN_correction;
	double **Ez_yzP_correction;

	double **Hx_xyN_correction;
	double **Hx_xyP_correction;
	double **Hx_xzN_correction;
	double **Hx_xzP_correction;

	double **Hy_xyN_correction;
	double **Hy_xyP_correction;
	double **Hy_yzN_correction;
	double **Hy_yzP_correction;

	double **Hz_xzN_correction;
	double **Hz_xzP_correction;
	double **Hz_yzN_correction;
	double **Hz_yzP_correction;



protected: // debug
	FILE *pf_debug1, *pf_debug2;
	// OMP
	void Calcuate_inc_Field_E_OMP();
	void Calcuate_inc_Field_H_OMP();
	void Correction_X_Boundary_Exn_OMP();
	void Correction_X_Boundary_Exp_OMP();
	void Correction_Y_Boundary_Eyn_OMP();
	void Correction_Y_Boundary_Eyp_OMP();
	void Correction_Z_Boundary_Ezn_OMP();
	void Correction_Z_Boundary_Ezp_OMP();

	void Correction_X_Boundary_Hxn_OMP();
	void Correction_X_Boundary_Hxp_OMP();
	void Correction_Y_Boundary_Hyn_OMP();
	void Correction_Y_Boundary_Hyp_OMP();
	void Correction_Z_Boundary_Hzn_OMP();
	void Correction_Z_Boundary_Hzp_OMP();
public:
	void DApply_OMP();
	void HApply_OMP();
	void TFSF_Electric_Field_OMP(double Pulse);
	void TFSF_Magnetic_Field_OMP();
};

inline double Vector_component_Hxinc(double psi, double phi, double theta);
inline double Vector_component_Hyinc(double psi, double phi, double theta);
inline double Vector_component_Hzinc(double psi, double phi, double theta);
inline double Vector_component_Exinc(double psi, double phi, double theta);
inline double Vector_component_Eyinc(double psi, double phi, double theta);
inline double Vector_component_Ezinc(double psi, double phi, double theta);
template <typename T1, typename T2, typename T3, typename T4> double Calculate_Distance(T1 xx, T2 yy, T3 zz, T4 Theta, T4 Phi);

inline double Vector_component_Hxinc(double psi, double phi, double theta)
{
	double corrections = sin(psi)*sin(phi) + cos(psi)*cos(theta)*cos(phi);
	return corrections;
}
inline double Vector_component_Hyinc(double psi, double phi, double theta)
{
	double corrections = -sin(psi)*cos(phi) + cos(psi)*cos(theta)*sin(phi);
	return corrections;
}
inline double Vector_component_Hzinc(double psi, double phi, double theta)
{
	double corrections = -cos(psi)*sin(theta);
	return corrections;
}
inline double Vector_component_Exinc(double psi, double phi, double theta)
{
	double corrections = cos(psi)*sin(phi) - sin(psi)*cos(theta)*cos(phi);
	return corrections;
}
inline double Vector_component_Eyinc(double psi, double phi, double theta)
{
	double corrections = -cos(psi)*cos(phi) - sin(psi)*cos(theta)*sin(phi);
	return corrections;
}
inline double Vector_component_Ezinc(double psi, double phi, double theta)
{
	double corrections = sin(psi)*sin(theta);
	return corrections;
}

template <typename T1, typename T2, typename T3, typename T4> double Calculate_Distance(T1 xx, T2 yy, T3 zz, T4 Theta, T4 Phi)
{
	double DD = (xx)*sin(Theta)*cos(Phi) + (yy)*sin(Theta)*sin(Phi) + (zz)*cos(Theta);
	return DD;
}

#endif
