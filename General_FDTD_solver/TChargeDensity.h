#pragma once
#include "TSaveFDTD_P_pros_Comm.h"
#include "TEnvelop_Extract.h"

typedef struct chargeD_EJ
{
	double*** Ex_tt1, *** Ey_tt1, *** Ez_tt1;
	double*** Jx_tt1, *** Jy_tt1, *** Jz_tt1;
	double*** Ex_tt2, *** Ey_tt2, *** Ez_tt2;
	double*** Jx_tt2, *** Jy_tt2, *** Jz_tt2;

	complex<double>*** Ex_ri;
	complex<double>*** Ey_ri;
	complex<double>*** Ez_ri;
	complex<double>*** Jx_ri;
	complex<double>*** Jy_ri;
	complex<double>*** Jz_ri;

}chargeD_EJ;

class TSaveChargeDensity : public TSaveFDTD_P_pros_Comm
{
public: 
	TSaveChargeDensity();
	~TSaveChargeDensity() {};
	void Prepare2save();
	void Save(int Nststep_iter);
private:
	void Read_save_info();

	char Fname[512];
	double omega;
	complex <double> jomega;

	complex <double> *** Dx, *** Dy, *** Dz, *** pCharge;

// Time to Frequency
	chargeD_EJ pchargeDEJ;
	void assign_electric_flux_D();
	void Field_3D_ElectrifFluxD(int ss);
	void Calcuate_Amp_Phase();
	complex<double> AP_to_RI(double A, double P);

	double Complex_to_Tiem(complex<double> pFields, int ss);
	void Charge_Density();
	void Save_Charge_Density();

	void Field_3D_ChargeD_EJ_tt1(int Nststep_iter);
	void Field_3D_ChargeD_EJ_tt2(int Nststep_iter);

	char SaveMod[256];
};

