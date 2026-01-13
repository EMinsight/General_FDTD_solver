#include "TChargeDensity.h"

TSaveChargeDensity::TSaveChargeDensity()
{
	_mkdir("result\\Chargedensity");
}

void TSaveChargeDensity::Prepare2save()
{
	Read_save_info();
	assign_electric_flux_D();
	// Assing_save_length();
}

void TSaveChargeDensity::Read_save_info()
{
	ifstream ifs(JsonFileName);
	IStreamWrapper isw(ifs);
	Jdocument.ParseStream(isw);
	
	strcpy(SaveMod, Jdocument["SIMULATION_ENVIORONMENT"]["Chargedensity"]["SaveMod"].GetString());
	double* temp_coordx, * temp_coordy, * temp_coordz;
	double maxx, maxy, maxz;

	maxx = txsize + 1;
	maxy = tysize + 1;
	maxz = tzsize + 1;
	temp_coordx = Tptr1<double>(maxx);
	temp_coordy = Tptr1<double>(maxy);
	temp_coordz = Tptr1<double>(maxz);
	for (int ii = 0; ii < maxx; ii++)
		temp_coordx[ii] = Jdocument["MESH"]["CellLocation"]["xcoordinate"][ii].GetDouble();

	for (int ii = 0; ii < maxy; ii++)
		temp_coordy[ii] = Jdocument["MESH"]["CellLocation"]["ycoordinate"][ii].GetDouble();

	for (int ii = 0; ii < maxz; ii++)
		temp_coordz[ii] = Jdocument["MESH"]["CellLocation"]["zcoordinate"][ii].GetDouble();

	double temp = 0;

	pchargeDEJ.Ex_tt1 = Tptr3<double>(txsize, tysize + 1, tzsize + 1);
	pchargeDEJ.Ey_tt1 = Tptr3<double>(txsize + 1, tysize, tzsize + 1);
	pchargeDEJ.Ez_tt1 = Tptr3<double>(txsize + 1, tysize + 1, tzsize);
	pchargeDEJ.Jx_tt1 = Tptr3<double>(txsize, tysize + 1, tzsize + 1);
	pchargeDEJ.Jy_tt1 = Tptr3<double>(txsize + 1, tysize, tzsize + 1);
	pchargeDEJ.Jz_tt1 = Tptr3<double>(txsize + 1, tysize + 1, tzsize);

	pchargeDEJ.Ex_tt2 = Tptr3<double>(txsize, tysize + 1, tzsize + 1);
	pchargeDEJ.Ey_tt2 = Tptr3<double>(txsize + 1, tysize, tzsize + 1);
	pchargeDEJ.Ez_tt2 = Tptr3<double>(txsize + 1, tysize + 1, tzsize);
	pchargeDEJ.Jx_tt2 = Tptr3<double>(txsize, tysize + 1, tzsize + 1);
	pchargeDEJ.Jy_tt2 = Tptr3<double>(txsize + 1, tysize, tzsize + 1);
	pchargeDEJ.Jz_tt2 = Tptr3<double>(txsize + 1, tysize + 1, tzsize);

	memfree(&temp_coordx);
	memfree(&temp_coordy);
	memfree(&temp_coordz);
	ifs.close();

	set_steady_start(0.9);
	set_steady_end(0.91);

}

void TSaveChargeDensity::assign_electric_flux_D()
{
	Dx = Tptr3<complex<double>>(txsize,   tysize+1, tzsize+1);
	Dy = Tptr3<complex<double>>(txsize+1, tysize,   tzsize+1);
	Dz = Tptr3<complex<double>>(txsize+1, tysize+1, tzsize  );
	pCharge = Tptr3<complex<double>>(txsize, tysize, tzsize);
	omega = 2 * pi * f0;
	jomega = complex<double>(0, omega);
}

void TSaveChargeDensity::Field_3D_ChargeD_EJ_tt1(int Nststep_iter)
{
	for (int ii = 0; ii < txsize; ii++)
		for (int jj = 0; jj < tysize; jj++)
			for (int kk = 0; kk < tzsize; kk++)
			{
				Voxel_face_average(ii, jj, kk, "Field_3D_ChargeD_EJ_tt1");
				pchargeDEJ.Ex_tt1[ii][jj][kk] = ExV;
				pchargeDEJ.Ey_tt1[ii][jj][kk] = EyV;
				pchargeDEJ.Ez_tt1[ii][jj][kk] = EzV;
				pchargeDEJ.Jx_tt1[ii][jj][kk] = JxV;
				pchargeDEJ.Jy_tt1[ii][jj][kk] = JyV;
				pchargeDEJ.Jz_tt1[ii][jj][kk] = JzV;
			}
}

void TSaveChargeDensity::Field_3D_ChargeD_EJ_tt2(int Nststep_iter)
{
	for (int ii = 0; ii < txsize; ii++)
		for (int jj = 0; jj < tysize; jj++)
			for (int kk = 0; kk < tzsize; kk++)
			{
				Voxel_face_average(ii, jj, kk, "Field_3D_ChargeD_EJ_tt2");
				pchargeDEJ.Ex_tt2[ii][jj][kk] = ExV;
				pchargeDEJ.Ey_tt2[ii][jj][kk] = EyV;
				pchargeDEJ.Ez_tt2[ii][jj][kk] = EzV;
				pchargeDEJ.Jx_tt2[ii][jj][kk] = JxV;
				pchargeDEJ.Jy_tt2[ii][jj][kk] = JyV;
				pchargeDEJ.Jz_tt2[ii][jj][kk] = JzV;
			}
}

void TSaveChargeDensity::Field_3D_ElectrifFluxD(int ss) // only harmonic.
{
	double Ex_diff, Ey_diff, Ez_diff, Jx_diff, Jy_diff, Jz_diff;
	for (int ii = 0; ii < txsize; ii++)
	for (int jj = 0; jj < tysize; jj++)
	for (int kk = 0; kk < tzsize; kk++)
	{
		Dx[ii][jj][kk] =  (eps_0 * pMp->Epsr * pchargeDEJ.Ex_ri[ii][jj][kk] + pchargeDEJ.Jx_ri[ii][jj][kk] / jomega);
		Dy[ii][jj][kk] =  (eps_0 * pMp->Epsr * pchargeDEJ.Ey_ri[ii][jj][kk] + pchargeDEJ.Jy_ri[ii][jj][kk] / jomega);
		Dz[ii][jj][kk] =  (eps_0 * pMp->Epsr * pchargeDEJ.Ez_ri[ii][jj][kk] + pchargeDEJ.Jz_ri[ii][jj][kk] / jomega);
	}
}



double TSaveChargeDensity::Complex_to_Tiem(complex<double> pFields, int ss)
{
	complex<double> temp = pFields * complex<double>(cos(omega * ss * dt), sin(omega * ss * dt));
	return temp.real();
}

void TSaveChargeDensity::Charge_Density()
{
	// div D
	for (int ii = 0; ii < txsize-1; ii++)
	for (int jj = 0; jj < tysize-1; jj++)
	for (int kk = 0; kk < tzsize-1; kk++)
	{
		pCharge[ii][jj][kk] = (Dx[ii + 1][jj	 ][kk	 ]-Dx[ii][jj][kk]) / pFc->E_dx[ii] +
							  (Dy[ii	 ][jj + 1][kk	 ]-Dy[ii][jj][kk]) / pFc->E_dy[jj] +
							  (Dz[ii	 ][jj	 ][kk + 1]-Dz[ii][jj][kk]) / pFc->E_dz[kk];
	}
}


void TSaveChargeDensity::Save(int Nststep_iter)
{
	if (save_sta == Nststep_iter && strcmp(SaveMod, "Frequency") == 0)
	{
		tt1 = Nststep_iter;
		Field_3D_ChargeD_EJ_tt1(Nststep_iter);
	}
	if (save_end == Nststep_iter && strcmp(SaveMod, "Frequency") == 0)
	{
		tt2 = Nststep_iter;
		Field_3D_ChargeD_EJ_tt2(Nststep_iter);
		// E, J -> frequency
		Calcuate_Amp_Phase();
		// D(w; x,y,z) <-calculate
		Field_3D_ElectrifFluxD(Nststep_iter- save_sta);

		Charge_Density();
		Save_Charge_Density();
	}
}

void TSaveChargeDensity::Save_Charge_Density()
{
	sprintf(Fname, "result\\Chargedensity\\chargedensity.dat");
	FILE* pf = fopen(Fname, "wb");
	for (int ii = 0; ii < txsize; ii++)
	for (int jj = 0; jj < tysize; jj++)
		fwrite(pCharge[ii][jj], sizeof(complex<double>), tzsize, pf);

	fclose(pf);
	memfree(&pCharge, txsize, tysize);
}


void TSaveChargeDensity::Calcuate_Amp_Phase()
{
	double f1, f2;
	pchargeDEJ.Ex_ri = Tptr3<complex<double>>(txsize, tysize, tzsize);	
	for (int ii = 0; ii < txsize; ii++)
	for (int jj = 0; jj < tysize; jj++)
	for (int kk = 0; kk < tzsize; kk++)
	{
		f1 = pchargeDEJ.Ex_tt1[ii][jj][kk];
		f2 = pchargeDEJ.Ex_tt2[ii][jj][kk];

		//Phase = atan2(cos(2 * pi*f0*dt*tt2) - f2 / f1, sin(2 * pi*f0*tt2));
		Phase = atan2((f2 * sin(2 * pi * f0 * dt * tt1) - f1 * sin(2 * pi * f0 * dt * tt2)),
			(f1 * cos(2 * pi * f0 * dt * tt2) - f2 * cos(2 * pi * f0 * dt * tt1)));
		Amp = f1 / sin(2 * pi * f0 * dt * tt1 + Phase);

		Amp_Phase_Tunning();
		pchargeDEJ.Ex_ri[ii][jj][kk] = AP_to_RI(Amp, Phase);
	}
	memfree(&pchargeDEJ.Ex_tt1, txsize, tysize);
	memfree(&pchargeDEJ.Ex_tt2, txsize, tysize);

	pchargeDEJ.Ey_ri = Tptr3<complex<double>>(txsize, tysize, tzsize);
	for (int ii = 0; ii < txsize; ii++)
		for (int jj = 0; jj < tysize; jj++)
			for (int kk = 0; kk < tzsize; kk++)
			{
				f1 = pchargeDEJ.Ey_tt1[ii][jj][kk];
				f2 = pchargeDEJ.Ey_tt2[ii][jj][kk];

				//Phase = atan2(cos(2 * pi*f0*dt*tt2) - f2 / f1, sin(2 * pi*f0*tt2));
				Phase = atan2((f2 * sin(2 * pi * f0 * dt * tt1) - f1 * sin(2 * pi * f0 * dt * tt2)),
					(f1 * cos(2 * pi * f0 * dt * tt2) - f2 * cos(2 * pi * f0 * dt * tt1)));
				Amp = f1 / sin(2 * pi * f0 * dt * tt1 + Phase);

				Amp_Phase_Tunning();
				pchargeDEJ.Ey_ri[ii][jj][kk] = AP_to_RI(Amp, Phase);
			}
	memfree(&pchargeDEJ.Ey_tt1, txsize, tysize);
	memfree(&pchargeDEJ.Ey_tt2, txsize, tysize);

	pchargeDEJ.Ez_ri = Tptr3<complex<double>>(txsize, tysize, tzsize);
	for (int ii = 0; ii < txsize; ii++)
		for (int jj = 0; jj < tysize; jj++)
			for (int kk = 0; kk < tzsize; kk++)
			{
				f1 = pchargeDEJ.Ez_tt1[ii][jj][kk];
				f2 = pchargeDEJ.Ez_tt2[ii][jj][kk];

				//Phase = atan2(cos(2 * pi*f0*dt*tt2) - f2 / f1, sin(2 * pi*f0*tt2));
				Phase = atan2((f2 * sin(2 * pi * f0 * dt * tt1) - f1 * sin(2 * pi * f0 * dt * tt2)),
					(f1 * cos(2 * pi * f0 * dt * tt2) - f2 * cos(2 * pi * f0 * dt * tt1)));
				Amp = f1 / sin(2 * pi * f0 * dt * tt1 + Phase);

				Amp_Phase_Tunning();
				pchargeDEJ.Ez_ri[ii][jj][kk] = AP_to_RI(Amp, Phase);
			}
	memfree(&pchargeDEJ.Ez_tt1, txsize, tysize);
	memfree(&pchargeDEJ.Ez_tt2, txsize, tysize);

	pchargeDEJ.Jx_ri = Tptr3<complex<double>>(txsize, tysize, tzsize);
	for (int ii = 0; ii < txsize; ii++)
		for (int jj = 0; jj < tysize; jj++)
			for (int kk = 0; kk < tzsize; kk++)
			{
				f1 = pchargeDEJ.Jx_tt1[ii][jj][kk];
				f2 = pchargeDEJ.Jx_tt2[ii][jj][kk];

				//Phase = atan2(cos(2 * pi*f0*dt*tt2) - f2 / f1, sin(2 * pi*f0*tt2));
				Phase = atan2((f2 * sin(2 * pi * f0 * dt * tt1) - f1 * sin(2 * pi * f0 * dt * tt2)),
					(f1 * cos(2 * pi * f0 * dt * tt2) - f2 * cos(2 * pi * f0 * dt * tt1)));
				Amp = f1 / sin(2 * pi * f0 * dt * tt1 + Phase);

				Amp_Phase_Tunning();
				pchargeDEJ.Jx_ri[ii][jj][kk] = AP_to_RI(Amp, Phase);
			}
	memfree(&pchargeDEJ.Jx_tt1, txsize, tysize);
	memfree(&pchargeDEJ.Jx_tt2, txsize, tysize);

	pchargeDEJ.Jy_ri = Tptr3<complex<double>>(txsize, tysize, tzsize);
	for (int ii = 0; ii < txsize; ii++)
		for (int jj = 0; jj < tysize; jj++)
			for (int kk = 0; kk < tzsize; kk++)
			{
				f1 = pchargeDEJ.Jy_tt1[ii][jj][kk];
				f2 = pchargeDEJ.Jy_tt2[ii][jj][kk];

				//Phase = atan2(cos(2 * pi*f0*dt*tt2) - f2 / f1, sin(2 * pi*f0*tt2));
				Phase = atan2((f2 * sin(2 * pi * f0 * dt * tt1) - f1 * sin(2 * pi * f0 * dt * tt2)),
					(f1 * cos(2 * pi * f0 * dt * tt2) - f2 * cos(2 * pi * f0 * dt * tt1)));
				Amp = f1 / sin(2 * pi * f0 * dt * tt1 + Phase);

				Amp_Phase_Tunning();
				pchargeDEJ.Jy_ri[ii][jj][kk] = AP_to_RI(Amp, Phase);
			}
	memfree(&pchargeDEJ.Jy_tt1, txsize, tysize);
	memfree(&pchargeDEJ.Jy_tt2, txsize, tysize);

	pchargeDEJ.Jz_ri = Tptr3<complex<double>>(txsize, tysize, tzsize);
	for (int ii = 0; ii < txsize; ii++)
		for (int jj = 0; jj < tysize; jj++)
			for (int kk = 0; kk < tzsize; kk++)
			{
				f1 = pchargeDEJ.Jz_tt1[ii][jj][kk];
				f2 = pchargeDEJ.Jz_tt2[ii][jj][kk];

				//Phase = atan2(cos(2 * pi*f0*dt*tt2) - f2 / f1, sin(2 * pi*f0*tt2));
				Phase = atan2((f2 * sin(2 * pi * f0 * dt * tt1) - f1 * sin(2 * pi * f0 * dt * tt2)),
					(f1 * cos(2 * pi * f0 * dt * tt2) - f2 * cos(2 * pi * f0 * dt * tt1)));
				Amp = f1 / sin(2 * pi * f0 * dt * tt1 + Phase);

				Amp_Phase_Tunning();
				pchargeDEJ.Jz_ri[ii][jj][kk] = AP_to_RI(Amp, Phase);
			}
	memfree(&pchargeDEJ.Jz_tt1, txsize, tysize);
	memfree(&pchargeDEJ.Jz_tt2, txsize, tysize);
}

complex<double> TSaveChargeDensity::AP_to_RI(double A, double P)
{
	complex<double> temp = A * exp(complex<double>(0, 1) * P);
	return temp;
}