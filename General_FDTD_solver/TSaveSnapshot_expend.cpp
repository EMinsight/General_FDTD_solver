#include "TSaveSnapshot_expend.h"
#define CurrentDensityOn 1
#define ElectricFieldsOn 1

TSaveSnapshot_expend::TSaveSnapshot_expend(struct_Field_carrier* sFc, Material_properties* pMp)
{
	pi = acos(-1);
	this->EOx = sFc->Ex_address;
	this->EOy = sFc->Ey_address;
	this->EOz = sFc->Ez_address;

	this->HOx = sFc->Hx_address;
	this->HOy = sFc->Hy_address;
	this->HOz = sFc->Hz_address;

	this->Media = sFc->media;
	this->Media_Ex = sFc->Media_Ex;
	this->Media_Ey = sFc->Media_Ey;
	this->Media_Ez = sFc->Media_Ez;

	this->pMp = pMp;
}

void TSaveSnapshot_expend::set_Frequency_range(double f0, double sfre, double tfre, int sampling)
{
	this->f0 = f0;
	this->sfre = sfre;
	this->tfre = tfre;
	this->sampling = sampling;
	this->df = round((tfre - sfre) / (double)sampling);
	Exctract_Metal_surface();
	MemoryAssign_for_FieldDistribution();
}

void TSaveSnapshot_expend::MemoryAssign_for_FieldDistribution()
{
	// 주파수 별로 sopds 개수 할당,
	// 각 주파수의 평면값 할당.
	// 공통 : 주파수, 평면.
	sopds = Tptr2<struct_Obs_Plane_datasave>(sampling, no_planes);
	for (int ff = 0; ff < sampling; ff++)
	{
		for (int pp = 0; pp < no_planes; pp++)
		{
			switch (sopd[pp].s_plane)
			{
			case 1:
				sopds[ff][pp].Ex_save_plane1 = Tptr2<double>(txsize, tysize);
				sopds[ff][pp].Ey_save_plane1 = Tptr2<double>(txsize, tysize);
				sopds[ff][pp].Ez_save_plane1 = Tptr2<double>(txsize, tysize);
				sopds[ff][pp].Ex_save_plane2 = Tptr2<double>(txsize, tysize);
				sopds[ff][pp].Ey_save_plane2 = Tptr2<double>(txsize, tysize);
				sopds[ff][pp].Ez_save_plane2 = Tptr2<double>(txsize, tysize);

				sopds[ff][pp].AmpEx = Tptr2<double>(txsize, tysize);
				sopds[ff][pp].AmpEy = Tptr2<double>(txsize, tysize);
				sopds[ff][pp].AmpEz = Tptr2<double>(txsize, tysize);
				sopds[ff][pp].PhaseEx = Tptr2<double>(txsize, tysize);
				sopds[ff][pp].PhaseEy = Tptr2<double>(txsize, tysize);
				sopds[ff][pp].PhaseEz = Tptr2<double>(txsize, tysize);
				break;
			case 2:
				sopds[ff][pp].Ex_save_plane1 = Tptr2<double>(tysize, tzsize);
				sopds[ff][pp].Ey_save_plane1 = Tptr2<double>(tysize, tzsize);
				sopds[ff][pp].Ez_save_plane1 = Tptr2<double>(tysize, tzsize);
				sopds[ff][pp].Ex_save_plane2 = Tptr2<double>(tysize, tzsize);
				sopds[ff][pp].Ey_save_plane2 = Tptr2<double>(tysize, tzsize);
				sopds[ff][pp].Ez_save_plane2 = Tptr2<double>(tysize, tzsize);
				sopds[ff][pp].AmpEx = Tptr2<double>(tysize, tzsize);
				sopds[ff][pp].AmpEy = Tptr2<double>(tysize, tzsize);
				sopds[ff][pp].AmpEz = Tptr2<double>(tysize, tzsize);
				sopds[ff][pp].PhaseEx = Tptr2<double>(tysize, tzsize);
				sopds[ff][pp].PhaseEy = Tptr2<double>(tysize, tzsize);
				sopds[ff][pp].PhaseEz = Tptr2<double>(tysize, tzsize);

				break;
			case 3:
				sopds[ff][pp].Ex_save_plane1 = Tptr2<double>(txsize, tzsize);
				sopds[ff][pp].Ey_save_plane1 = Tptr2<double>(txsize, tzsize);
				sopds[ff][pp].Ez_save_plane1 = Tptr2<double>(txsize, tzsize);
				sopds[ff][pp].Ex_save_plane2 = Tptr2<double>(txsize, tzsize);
				sopds[ff][pp].Ey_save_plane2 = Tptr2<double>(txsize, tzsize);
				sopds[ff][pp].Ez_save_plane2 = Tptr2<double>(txsize, tzsize);
				sopds[ff][pp].AmpEx = Tptr2<double>(txsize, tzsize);
				sopds[ff][pp].AmpEy = Tptr2<double>(txsize, tzsize);
				sopds[ff][pp].AmpEz = Tptr2<double>(txsize, tzsize);
				sopds[ff][pp].PhaseEx = Tptr2<double>(txsize, tzsize);
 				sopds[ff][pp].PhaseEy = Tptr2<double>(txsize, tzsize);
				sopds[ff][pp].PhaseEz = Tptr2<double>(txsize, tzsize);
				break;
			}
		}
	}
	pSave_Temp_re = Tptr1<double>(sampling);
	pSave_Temp_im = Tptr1<double>(sampling);
	pFourier = new TFourierTransform();
	pFourier->setParameters(f0, f0, sampling, dt);
	
	
	// Current density
	if (CurrentDensityOn)
	{
		Jup = Tptr2<struct_Current_density>(txsize, tysize);
		Jdown = Tptr2<struct_Current_density>(txsize, tysize);
	}

}

void TSaveSnapshot_expend::Calcuate_Amp_PhaseXY(int ff, int pp)
{
	double f1, f2;
	for (int ii = 0; ii < txsize; ii++)
		for (int jj = 0; jj < tysize; jj++)
		{
			f1 = sopds[ff][pp].Ex_save_plane1[ii][jj];
			f2 = sopds[ff][pp].Ex_save_plane2[ii][jj];
			//Phase = atan2(cos(2 * pi*f0*dt*tt2) - f2 / f1, sin(2 * pi*f0*tt2));
			Phase = atan2((f2 * sin(2 * pi * f0 * dt * tt1) - f1 * sin(2 * pi * f0 * dt * tt2)),
				(f1 * cos(2 * pi * f0 * dt * tt2) - f2 * cos(2 * pi * f0 * dt * tt1)));
			Amp = f1 / sin(2 * pi * f0 * dt * tt1 + Phase);
			Amp_Phase_Tunning();
			sopds[ff][pp].AmpEx[ii][jj] = Amp;
			sopds[ff][pp].PhaseEx[ii][jj] = Phase;


			f1 = sopds[ff][pp].Ey_save_plane1[ii][jj];
			f2 = sopds[ff][pp].Ey_save_plane2[ii][jj];
			//Phase = atan2(cos(2 * pi*f0*dt*tt2) - f2 / f1, sin(2 * pi*f0*tt2));
			Phase = atan2((f2 * sin(2 * pi * f0 * dt * tt1) - f1 * sin(2 * pi * f0 * dt * tt2)),
				(f1 * cos(2 * pi * f0 * dt * tt2) - f2 * cos(2 * pi * f0 * dt * tt1)));
			Amp = f1 / sin(2 * pi * f0 * dt * tt1 + Phase);
			Amp_Phase_Tunning();
			sopds[ff][pp].AmpEy[ii][jj] = Amp;
			sopds[ff][pp].PhaseEy[ii][jj] = Phase;


			f1 = sopds[ff][pp].Ez_save_plane1[ii][jj];
			f2 = sopds[ff][pp].Ez_save_plane2[ii][jj];
			//Phase = atan2(cos(2 * pi*f0*dt*tt2) - f2 / f1, sin(2 * pi*f0*tt2));
			Phase = atan2((f2 * sin(2 * pi * f0 * dt * tt1) - f1 * sin(2 * pi * f0 * dt * tt2)),
				(f1 * cos(2 * pi * f0 * dt * tt2) - f2 * cos(2 * pi * f0 * dt * tt1)));
			Amp = f1 / sin(2 * pi * f0 * dt * tt1 + Phase);
			Amp_Phase_Tunning();
			sopds[ff][pp].AmpEz[ii][jj] = Amp;
			sopds[ff][pp].PhaseEz[ii][jj] = Phase;
		}
}
void TSaveSnapshot_expend::Calcuate_Amp_PhaseYZ(int ff, int pp)
{
	double f1, f2;
	for (int jj = 0; jj < tysize; jj++)
		for (int kk = 0; kk < tzsize; kk++)
		{
			f1 = sopds[ff][pp].Ex_save_plane1[jj][kk];
			f2 = sopds[ff][pp].Ex_save_plane2[jj][kk];
			//Phase = atan2(cos(2 * pi*f0*dt*tt2) - f2 / f1, sin(2 * pi*f0*tt2));
			Phase = atan2((f2 * sin(2 * pi * f0 * dt * tt1) - f1 * sin(2 * pi * f0 * dt * tt2)),
				(f1 * cos(2 * pi * f0 * dt * tt2) - f2 * cos(2 * pi * f0 * dt * tt1)));
			Amp = f1 / sin(2 * pi * f0 * dt * tt1 + Phase);
			Amp_Phase_Tunning();
			sopds[ff][pp].AmpEx[jj][kk] = Amp;
			sopds[ff][pp].PhaseEx[jj][kk] = Phase;


			f1 = sopds[ff][pp].Ey_save_plane1[jj][kk];
			f2 = sopds[ff][pp].Ey_save_plane2[jj][kk];
			//Phase = atan2(cos(2 * pi*f0*dt*tt2) - f2 / f1, sin(2 * pi*f0*tt2));
			Phase = atan2((f2 * sin(2 * pi * f0 * dt * tt1) - f1 * sin(2 * pi * f0 * dt * tt2)),
				(f1 * cos(2 * pi * f0 * dt * tt2) - f2 * cos(2 * pi * f0 * dt * tt1)));
			Amp = f1 / sin(2 * pi * f0 * dt * tt1 + Phase);
			Amp_Phase_Tunning();
			sopds[ff][pp].AmpEy[jj][kk] = Amp;
			sopds[ff][pp].PhaseEy[jj][kk] = Phase;


			f1 = sopds[ff][pp].Ez_save_plane1[jj][kk];
			f2 = sopds[ff][pp].Ez_save_plane2[jj][kk];
			//Phase = atan2(cos(2 * pi*f0*dt*tt2) - f2 / f1, sin(2 * pi*f0*tt2));
			Phase = atan2((f2 * sin(2 * pi * f0 * dt * tt1) - f1 * sin(2 * pi * f0 * dt * tt2)),
				(f1 * cos(2 * pi * f0 * dt * tt2) - f2 * cos(2 * pi * f0 * dt * tt1)));
			Amp = f1 / sin(2 * pi * f0 * dt * tt1 + Phase);
			Amp_Phase_Tunning();
			sopds[ff][pp].AmpEz[jj][kk] = Amp;
			sopds[ff][pp].PhaseEz[jj][kk] = Phase;
		}
}
void TSaveSnapshot_expend::Calcuate_Amp_PhaseXZ(int ff, int pp)
{
	double f1, f2;
	for (int ii = 0; ii < txsize; ii++)
		for (int kk = 0; kk < tzsize; kk++)
		{
			f1 = sopds[ff][pp].Ex_save_plane1[ii][kk];
			f2 = sopds[ff][pp].Ex_save_plane2[ii][kk];
			//Phase = atan2(cos(2 * pi*f0*dt*tt2) - f2 / f1, sin(2 * pi*f0*tt2));
			Phase = atan2((f2 * sin(2 * pi * f0 * dt * tt1) - f1 * sin(2 * pi * f0 * dt * tt2)),
				(f1 * cos(2 * pi * f0 * dt * tt2) - f2 * cos(2 * pi * f0 * dt * tt1)));
			Amp = f1 / sin(2 * pi * f0 * dt * tt1 + Phase);
			Amp_Phase_Tunning();
			sopds[ff][pp].AmpEx[ii][kk] = Amp;
			sopds[ff][pp].PhaseEx[ii][kk] = Phase;


			f1 = sopds[ff][pp].Ey_save_plane1[ii][kk];
			f2 = sopds[ff][pp].Ey_save_plane2[ii][kk];
			//Phase = atan2(cos(2 * pi*f0*dt*tt2) - f2 / f1, sin(2 * pi*f0*tt2));
			Phase = atan2((f2 * sin(2 * pi * f0 * dt * tt1) - f1 * sin(2 * pi * f0 * dt * tt2)),
				(f1 * cos(2 * pi * f0 * dt * tt2) - f2 * cos(2 * pi * f0 * dt * tt1)));
			Amp = f1 / sin(2 * pi * f0 * dt * tt1 + Phase);
			Amp_Phase_Tunning();
			sopds[ff][pp].AmpEy[ii][kk] = Amp;
			sopds[ff][pp].PhaseEy[ii][kk] = Phase;


			f1 = sopds[ff][pp].Ez_save_plane1[ii][kk];
			f2 = sopds[ff][pp].Ez_save_plane2[ii][kk];
			//Phase = atan2(cos(2 * pi*f0*dt*tt2) - f2 / f1, sin(2 * pi*f0*tt2));
			Phase = atan2((f2 * sin(2 * pi * f0 * dt * tt1) - f1 * sin(2 * pi * f0 * dt * tt2)),
				(f1 * cos(2 * pi * f0 * dt * tt2) - f2 * cos(2 * pi * f0 * dt * tt1)));
			Amp = f1 / sin(2 * pi * f0 * dt * tt1 + Phase);
			Amp_Phase_Tunning();
			sopds[ff][pp].AmpEz[ii][kk] = Amp;
			sopds[ff][pp].PhaseEz[ii][kk] = Phase;
		}
}

void TSaveSnapshot_expend::Calcuate_Amp_Phase()
{
	for (int ff = 0; ff < sampling; ff++)
		for (int pp = 0; pp < no_planes; pp++)
		{
			switch (sopd[pp].s_plane)
			{
			case 1:
				Calcuate_Amp_PhaseXY(ff, pp);
				break;
			case 2:
				Calcuate_Amp_PhaseYZ(ff, pp);
				break;
			case 3:
				Calcuate_Amp_PhaseXZ(ff, pp);
				break;
			}
		}
}

void TSaveSnapshot_expend::store_snapshot_Time(int Nststep_iter)
{
	int steady_state_offset = 0;
	this->Nststep_iter = Nststep_iter;
	complex <double> temp;
	double f1, f2;
	if (Nststep_iter == steady_state_start_discrete)
	{
		tt1 = Nststep_iter;
		for (int ff = 0; ff < sampling; ff++)
			for (int pp = 0; pp < no_planes; pp++)
			{
				switch (sopd[pp].s_plane)
				{
				case 1:
					SaveFieldXY(ff, pp, 1);
					break;
				case 2:
					SaveFieldYZ(ff, pp, 1);
					break;
				case 3:
					SaveFieldXZ(ff, pp, 1);
					break;
				}
			}
	}
	if (Nststep_iter == round(steady_state_start_discrete + 0.1 * T0_discrete))
	{
		// optical chirality 계산 가능.
		tt2 = Nststep_iter - tt1;
		for (int ff = 0; ff < sampling; ff++)
			for (int pp = 0; pp < no_planes; pp++)
			{
				switch (sopd[pp].s_plane)
				{
				case 1:
					SaveFieldXY(ff, pp, 1);
					break;
				case 2:
					SaveFieldYZ(ff, pp, 1);
					break;
				case 3:
					SaveFieldXZ(ff, pp, 1);
					break;
				}
			}

		Calcuate_Amp_Phase();
		save_FieldDistribution();
	}
}

void TSaveSnapshot_expend::store_current_Time(int Nststep_iter)
{
	int steady_state_offset = 0;
	this->Nststep_iter = Nststep_iter;
	complex <double> temp;
	double f1, f2;
	if (Nststep_iter == steady_state_start_discrete)
	{
		tt1 = Nststep_iter;
		for (int ii = 0; ii < txsize; ii++)
			for (int jj = 0; jj < tysize; jj++)
			{
				Single_Freq_Amp_Phase_J1(ii, jj);
			}
	}
	if (Nststep_iter == round(steady_state_start_discrete + 0.1 * T0_discrete))
	{
		// optical chirality 계산 가능.
		tt2 = Nststep_iter;
		for (int ii = 0; ii < txsize; ii++)
			for (int jj = 0; jj < tysize; jj++)
			{
				Single_Freq_Amp_Phase_J2(ii, jj);
			}

		save_nXh();
		for (int ii = 0; ii < txsize; ii++)
			for (int jj = 0; jj < tysize; jj++)
			{
				for (int kk = 0; kk < 3; kk++)
				{
					f1 = Jup[ii][jj].J1[kk];
					f2 = Jup[ii][jj].J2[kk];
					//Phase = atan2(cos(2 * pi*f0*dt*tt2) - f2 / f1, sin(2 * pi*f0*tt2));
					Phase = atan2((f2 * sin(2 * pi * f0 * dt * tt1) - f1 * sin(2 * pi * f0 * dt * tt2)),
						(f1 * cos(2 * pi * f0 * dt * tt2) - f2 * cos(2 * pi * f0 * dt * tt1)));
					Amp = f1 / sin(2 * pi * f0 * dt * tt1 + Phase);
					Amp_Phase_Tunning();
					Jup[ii][jj].Amp[kk] = Amp;
					Jup[ii][jj].Phase[kk] = Phase;
				}
				for (int kk = 0; kk < 3; kk++)
				{
					f1 = Jdown[ii][jj].J1[kk];
					f2 = Jdown[ii][jj].J2[kk];
					//Phase = atan2(cos(2 * pi*f0*dt*tt2) - f2 / f1, sin(2 * pi*f0*tt2));
					Phase = atan2((f2 * sin(2 * pi * f0 * dt * tt1) - f1 * sin(2 * pi * f0 * dt * tt2)),
						(f1 * cos(2 * pi * f0 * dt * tt2) - f2 * cos(2 * pi * f0 * dt * tt1)));
					Amp = f1 / sin(2 * pi * f0 * dt * tt1 + Phase);
					Amp_Phase_Tunning();
					Jdown[ii][jj].Amp[kk] = Amp;
					Jdown[ii][jj].Phase[kk] = Phase;
				}
			}
		save_Current_density();
	}
}

void TSaveSnapshot_expend::Exctract_Metal_surface()
{
	int Air = 0;
	int PDMS = 2;
	int Au = 3;
	int LL = 0;
	int Water = 1;
	

	//Get_Pulse_infomation();
	//Layer = no_composit - 2 + 3;
	// Num. Layer check;
	for (int ii = 0; ii < txsize; ii++)
	for (int jj = 0; jj < tysize; jj++)
	for (int kk = 0; kk < tzsize; kk++)
	{
		__try
		{
			if (Media[ii][jj][kk] - Media[ii - 1][jj][kk])
				;
			if (Media[ii][jj][kk] - Media[ii + 1][jj][kk])
				;
		}
		__except(NULL)
		{
			continue;
		}
		
		if (Media[ii][jj][kk])
			;


		// Media_Ex[ii][jj][kk];
		// Media_Ey[ii][jj][kk];
		// Media_Ez[ii][jj][kk];
	}

	

	pExtract_current_distribution_plane    = Tptr3<int>(2, txsize, tysize);
	pExtract_current_distribution_envelope = Tptr2<int>(txsize, tysize);

	for (int ii = 0; ii < txsize; ii++)
		for (int jj = 0; jj < tysize; jj++)
		{
			for (int kk = 0; kk < tzsize - 1; kk++)
			{
				if (Media[ii][jj][kk] != Media[ii][jj][kk + 1])
				{
					pExtract_current_distribution_envelope[ii][jj] = kk;
					break;
				}
			}
		}

	// 2차원 변환을 할 평면 일기 (= composit과 동일 위치)
	// Method 1.0 : PDMS와 Air, Au의 접합면 검출. (z 방향으로만)
	for (int ii = 0; ii < txsize; ii++)
		for (int jj = 0; jj < tysize; jj++)
		{
			LL = 0;
			for (int kk = 0; kk < tzsize - 1; kk++)
			{
				if (Media[ii][jj][kk] == PDMS || Media[ii][jj][kk] == Air)
				{
					// 카이랄 경계를 envelop.
					if (Media[ii][jj][kk - 1] == Au || Media[ii][jj][kk - 1] == Water)
						//if (Media[ii][jj][kk + 1] == Au || Media[ii][jj][kk + 1] == Air || Media[ii][jj][kk + 1] == Water)
					{
						pExtract_current_distribution_plane[LL][ii][jj] = kk;
					}
				}
				else
				{
					continue;
				}
			}
		}
	for (int ii = 0; ii < txsize; ii++)
		for (int jj = 0; jj < tysize; jj++)
		{
			LL = 1;
			for (int kk = 0; kk < tzsize - 1; kk++)
			{
				if (Media[ii][jj][kk] == PDMS || Media[ii][jj][kk] == Air)
				{
					// 카이랄 경계를 envelop.
					//if (Media[ii][jj][kk - 1] == Au || Media[ii][jj][kk - 1] == Air) // || Media[ii][jj][kk + 1] == Water)
					if (Media[ii][jj][kk + 1] == Au || Media[ii][jj][kk + 1] == Water)
					{
						pExtract_current_distribution_plane[LL][ii][jj] = kk;
					}
				}
				else
				{
					continue;
				}
			}
		}

	char LayerName[256];
	for (int LL = 0; LL < 2; LL++)
	{
		sprintf(LayerName, "./2D_current[%d].dat", LL);
		FILE* p2D_optical_chirality = fopen(LayerName, "wb");
		for (int ii = 0; ii < txsize; ii++)
			fwrite(pExtract_current_distribution_plane[LL][ii], tysize, sizeof(int), p2D_optical_chirality);
		fclose(p2D_optical_chirality);
	}
}

void TSaveSnapshot_expend::Single_Freq_Amp_Phase_J1(int xx, int yy)
{
	int zz = 0, nn = 0;
	//#pragma omp parallel private(zz, nn, ff)
	zz = pExtract_current_distribution_plane[nn][xx][yy];
	Jup[xx][yy].Hx1 = (HOx[xx][yy][zz]);//+HOx[xx + 1][yy][zz]) / 2.0;
	Jup[xx][yy].Hy1 = (HOy[xx][yy][zz]);//+HOy[xx][yy + 1][zz]) / 2.0;
	Jup[xx][yy].Hz1 = (HOz[xx][yy][zz]); //+HOz[xx][yy][zz + 1]) / 2.0;
	if (xx + 1 < txsize)
	{
		Jup[xx][yy].Hydx1 = HOy[xx + 1][yy][zz];
		Jup[xx][yy].Hzdx1 = HOz[xx + 1][yy][zz];
	}
	if (yy + 1 < tysize)
	{
		Jup[xx][yy].Hxdy1 = HOx[xx][yy + 1][zz];
		Jup[xx][yy].Hzdy1 = HOz[xx][yy + 1][zz];
	}
	if (zz + 1 < tzsize)
	{
		Jup[xx][yy].Hxdz1 = HOx[xx][yy][zz + 1];
		Jup[xx][yy].Hydz1 = HOy[xx][yy][zz + 1];
	}

	zz = pExtract_current_distribution_plane[1][xx][yy];
	Jdown[xx][yy].Hx1 = (HOx[xx][yy][zz]);
	Jdown[xx][yy].Hy1 = (HOy[xx][yy][zz]);
	Jdown[xx][yy].Hz1 = (HOz[xx][yy][zz]); //+ HOz[xx][yy][zz + 1]) / 2.0;
	if (xx + 1 < txsize)
	{
		Jdown[xx][yy].Hydx1 = HOy[xx + 1][yy][zz];
		Jdown[xx][yy].Hzdx1 = HOz[xx + 1][yy][zz];
	}
	if (yy + 1 < tysize)
	{
		Jdown[xx][yy].Hxdy1 = HOx[xx][yy + 1][zz];
		Jdown[xx][yy].Hzdy1 = HOz[xx][yy + 1][zz];
	}
	if (zz + 1 < tzsize)
	{
		Jdown[xx][yy].Hxdz1 = HOx[xx][yy][zz + 1];
		Jdown[xx][yy].Hydz1 = HOy[xx][yy][zz + 1];
	}
}
void TSaveSnapshot_expend::Single_Freq_Amp_Phase_J2(int xx, int yy)
{
	int zz = 0, nn = 0;
	//#pragma omp parallel private(zz, nn, ff)
	zz = pExtract_current_distribution_plane[nn][xx][yy];
	Jup[xx][yy].Hx2 = (HOx[xx][yy][zz]);// + HOx[xx + 1][yy][zz]) / 2.0;
	Jup[xx][yy].Hy2 = (HOy[xx][yy][zz]);// + HOy[xx][yy + 1][zz]) / 2.0;
	Jup[xx][yy].Hz2 = (HOz[xx][yy][zz]);// + HOz[xx][yy][zz + 1]) / 2.0;
	if (xx + 1 < txsize)
	{
		Jup[xx][yy].Hydx2 = HOy[xx + 1][yy][zz];
		Jup[xx][yy].Hzdx2 = HOz[xx + 1][yy][zz];
	}
	if (yy + 1 < tysize)
	{
		Jup[xx][yy].Hxdy2 = HOx[xx][yy + 1][zz];
		Jup[xx][yy].Hzdy2 = HOz[xx][yy + 1][zz];
	}
	if (zz + 1 < tzsize)
	{
		Jup[xx][yy].Hxdz2 = HOx[xx][yy][zz + 1];
		Jup[xx][yy].Hydz2 = HOy[xx][yy][zz + 1];
	}

	zz = pExtract_current_distribution_plane[1][xx][yy];
	Jdown[xx][yy].Hx2 = (HOx[xx][yy][zz]);//+ HOx[xx + 1][yy][zz]) / 2.0;
	Jdown[xx][yy].Hy2 = (HOy[xx][yy][zz]);//+ HOy[xx][yy + 1][zz]) / 2.0;
	Jdown[xx][yy].Hz2 = (HOz[xx][yy][zz]);//+ HOz[xx][yy][zz + 1]) / 2.0;
	if (xx + 1 < txsize)
	{
		Jdown[xx][yy].Hydx2 = HOy[xx + 1][yy][zz];
		Jdown[xx][yy].Hzdx2 = HOz[xx + 1][yy][zz];
	}
	if (yy + 1 < tysize)
	{
		Jdown[xx][yy].Hxdy2 = HOx[xx][yy + 1][zz];
		Jdown[xx][yy].Hzdy2 = HOz[xx][yy + 1][zz];
	}
	if (zz + 1 < tzsize)
	{
		Jdown[xx][yy].Hxdz2 = HOx[xx][yy][zz + 1];
		Jdown[xx][yy].Hydz2 = HOy[xx][yy][zz + 1];
	}
}

void TSaveSnapshot_expend::SaveFieldXY(int ff, int pp, int t1t2)
{
	if (t1t2 == 1)
		for (int ii = 0; ii < txsize; ii++)
			for (int jj = 0; jj < tysize; jj++)
			{
				sopds[ff][pp].Ex_save_plane1[ii][jj] = (EOx[ii][jj][sopd[pp].r_th] + EOx[ii][jj + 1][sopd[pp].r_th] +
					EOx[ii][jj][sopd[pp].r_th + 1] + EOx[ii][jj + 1][sopd[pp].r_th + 1]) / 4.0;
				sopds[ff][pp].Ey_save_plane1[ii][jj] = (EOy[ii][jj][sopd[pp].r_th] + EOy[ii + 1][jj][sopd[pp].r_th] +
					EOy[ii][jj][sopd[pp].r_th + 1] + EOy[ii + 1][jj][sopd[pp].r_th + 1]) / 4.0;
				sopds[ff][pp].Ez_save_plane1[ii][jj] = (EOz[ii][jj][sopd[pp].r_th] + EOz[ii + 1][jj][sopd[pp].r_th] +
					EOz[ii][jj + 1][sopd[pp].r_th + 1] + EOz[ii + 1][jj + 1][sopd[pp].r_th]) / 4.0;
			}

	if (t1t2 == 2)
		for (int ii = 0; ii < txsize; ii++)
			for (int jj = 0; jj < tysize; jj++)
			{
				sopds[ff][pp].Ex_save_plane2[ii][jj] = (EOx[ii][jj][sopd[pp].r_th] + EOx[ii][jj + 1][sopd[pp].r_th] +
					EOx[ii][jj][sopd[pp].r_th + 1] + EOx[ii][jj + 1][sopd[pp].r_th + 1]) / 4.0;
				sopds[ff][pp].Ey_save_plane2[ii][jj] = (EOy[ii][jj][sopd[pp].r_th] + EOy[ii + 1][jj][sopd[pp].r_th] +
					EOy[ii][jj][sopd[pp].r_th + 1] + EOy[ii + 1][jj][sopd[pp].r_th + 1]) / 4.0;
				sopds[ff][pp].Ez_save_plane2[ii][jj] = (EOz[ii][jj][sopd[pp].r_th] + EOz[ii + 1][jj][sopd[pp].r_th] +
					EOz[ii][jj + 1][sopd[pp].r_th + 1] + EOz[ii + 1][jj + 1][sopd[pp].r_th]) / 4.0;
			}
}

void TSaveSnapshot_expend::SaveFieldYZ(int ff, int pp, int t1t2)
{
	if (t1t2 == 1)
		for (int jj = 0; jj < tysize; jj++)
			for (int kk = 0; kk < tzsize; kk++)
			{
				sopds[ff][pp].Ex_save_plane1[jj][kk] = (EOx[sopd[pp].r_th][jj][kk] + EOx[sopd[pp].r_th][jj + 1][kk] +
					EOx[sopd[pp].r_th][jj][kk + 1] + EOx[sopd[pp].r_th][jj + 1][kk + 1]) / 4.0;
				sopds[ff][pp].Ey_save_plane1[jj][kk] = (EOy[sopd[pp].r_th][jj][kk] + EOy[sopd[pp].r_th + 1][jj][kk] +
					EOy[sopd[pp].r_th][jj][kk + 1] + EOy[sopd[pp].r_th + 1][jj][kk + 1]) / 4.0;
				sopds[ff][pp].Ez_save_plane1[jj][kk] = (EOz[sopd[pp].r_th][jj][kk] + EOz[sopd[pp].r_th + 1][jj][kk] +
					EOz[sopd[pp].r_th][jj][kk + 1] + EOz[sopd[pp].r_th + 1][jj + 1][kk]) / 4.0;
			}
	if (t1t2 == 2)
		for (int jj = 0; jj < tysize; jj++)
			for (int kk = 0; kk < tzsize; kk++)
			{
				sopds[ff][pp].Ex_save_plane2[jj][kk] = (EOx[sopd[pp].r_th][jj][kk] + EOx[sopd[pp].r_th][jj + 1][kk] +
					EOx[sopd[pp].r_th][jj][kk + 1] + EOx[sopd[pp].r_th][jj + 1][kk + 1]) / 4.0;
				sopds[ff][pp].Ey_save_plane2[jj][kk] = (EOy[sopd[pp].r_th][jj][kk] + EOy[sopd[pp].r_th + 1][jj][kk] +
					EOy[sopd[pp].r_th][jj][kk + 1] + EOy[sopd[pp].r_th + 1][jj][kk + 1]) / 4.0;
				sopds[ff][pp].Ez_save_plane2[jj][kk] = (EOz[sopd[pp].r_th][jj][kk] + EOz[sopd[pp].r_th + 1][jj][kk] +
					EOz[sopd[pp].r_th][jj][kk + 1] + EOz[sopd[pp].r_th + 1][jj + 1][kk]) / 4.0;
			}
}

void TSaveSnapshot_expend::SaveFieldXZ(int ff, int pp, int t1t2)
{
	if (t1t2 == 1)
		for (int ii = 0; ii < txsize; ii++)
			for (int kk = 0; kk < tzsize; kk++)
			{
				sopds[ff][pp].Ex_save_plane1[ii][kk] = (EOx[ii][sopd[pp].r_th][kk] + EOx[ii][sopd[pp].r_th + 1][kk] +
					EOx[ii][sopd[pp].r_th][kk + 1] + EOx[ii][sopd[pp].r_th + 1][kk + 1]) / 4.0;
				sopds[ff][pp].Ey_save_plane1[ii][kk] = (EOy[ii][sopd[pp].r_th][kk] + EOy[ii + 1][sopd[pp].r_th][kk] +
					EOy[ii][sopd[pp].r_th][kk + 1] + EOy[ii + 1][sopd[pp].r_th][kk + 1]) / 4.0;
				sopds[ff][pp].Ez_save_plane1[ii][kk] = (EOz[ii][sopd[pp].r_th][kk] + EOz[ii + 1][sopd[pp].r_th][kk] +
					EOz[ii][sopd[pp].r_th + 1][kk + 1] + EOz[ii + 1][sopd[pp].r_th + 1][kk]) / 4.0;
			}
	if (t1t2 == 2)
		for (int ii = 0; ii < txsize; ii++)
			for (int kk = 0; kk < tzsize; kk++)
			{
				sopds[ff][pp].Ex_save_plane2[ii][kk] = (EOx[ii][sopd[pp].r_th][kk] + EOx[ii][sopd[pp].r_th + 1][kk] +
					EOx[ii][sopd[pp].r_th][kk + 1] + EOx[ii][sopd[pp].r_th + 1][kk + 1]) / 4.0;
				sopds[ff][pp].Ey_save_plane2[ii][kk] = (EOy[ii][sopd[pp].r_th][kk] + EOy[ii + 1][sopd[pp].r_th][kk] +
					EOy[ii][sopd[pp].r_th][kk + 1] + EOy[ii + 1][sopd[pp].r_th][kk + 1]) / 4.0;
				sopds[ff][pp].Ez_save_plane2[ii][kk] = (EOz[ii][sopd[pp].r_th][kk] + EOz[ii + 1][sopd[pp].r_th][kk] +
					EOz[ii][sopd[pp].r_th + 1][kk + 1] + EOz[ii + 1][sopd[pp].r_th + 1][kk]) / 4.0;
			}
}

void TSaveSnapshot_expend::Voxel_center_average(int xx, int yy, int zz)
{
	Ex_of = (EOx[xx][yy][zz] +
		EOx[xx][yy + 1][zz] +
		EOx[xx][yy][zz + 1] +
		EOx[xx][yy + 1][zz + 1]) / 4.0;
	Ey_of = (EOy[xx][yy][zz] +
		EOy[xx + 1][yy][zz] +
		EOy[xx][yy][zz + 1] +
		EOy[xx + 1][yy][zz + 1]) / 4.0;
	Ez_of = (EOz[xx][yy][zz] +
		EOz[xx + 1][yy][zz] +
		EOz[xx][yy + 1][zz] +
		EOz[xx + 1][yy + 1][zz]) / 4.0;

	Hx_of = (HOx[xx][yy][zz] + HOx[xx + 1][yy][zz]) / 2.0;
	Hy_of = (HOy[xx][yy][zz] + HOy[xx][yy + 1][zz]) / 2.0;
	Hz_of = (HOz[xx][yy][zz] + HOz[xx][yy][zz + 1]) / 2.0;
}

void TSaveSnapshot_expend::save_FieldDistribution_make_saveFolder(int kk)
{
	switch (sopd[kk].s_plane)
	{
	case 1:
		sprintf(Fname_save_field_Path, "./result/Field_distribution_XY_[%04d]", sopd[kk].r_th);
		_mkdir(Fname_save_field_Path);
		break;
	case 2:
		sprintf(Fname_save_field_Path, "./result/Field_distribution_YZ_[%04d]", sopd[kk].r_th);
		_mkdir(Fname_save_field_Path);
		break;
	case 3:
		sprintf(Fname_save_field_Path, "./result/Field_distribution_XZ_[%04d]", sopd[kk].r_th);
		_mkdir(Fname_save_field_Path);
		break;
	}
}

void TSaveSnapshot_expend::save_FieldDistribution_assign_file_pointer(int ff)
{
	strcpy(Fname_save_field_Pathcpy_amp, Fname_save_field_Path);
	strcpy(Fname_save_field_Pathcpy_phase, Fname_save_field_Path);
	sprintf(Fname_save_field_Name_temp, "/Ex_Amp[%1.2le].dat", sfre + ff * df);
	strcpy(Fname_save_field_Amp, strcat(Fname_save_field_Pathcpy_amp, Fname_save_field_Name_temp));
	sprintf(Fname_save_field_Name_temp, "/Ex_Phase[%1.2le].dat", sfre + ff * df);
	strcpy(Fname_save_field_Phase, strcat(Fname_save_field_Pathcpy_phase, Fname_save_field_Name_temp));
	pFEx_Amp = fopen(Fname_save_field_Amp, "wb");
	pFEx_Phase = fopen(Fname_save_field_Phase, "wb");

	strcpy(Fname_save_field_Pathcpy_amp, Fname_save_field_Path);
	strcpy(Fname_save_field_Pathcpy_phase, Fname_save_field_Path);
	sprintf(Fname_save_field_Name_temp, "/Ey_Amp[%1.2le].dat", sfre + ff * df);
	strcpy(Fname_save_field_Amp, strcat(Fname_save_field_Pathcpy_amp, Fname_save_field_Name_temp));
	sprintf(Fname_save_field_Name_temp, "/Ey_Phase[%1.2le].dat", sfre + ff * df);
	strcpy(Fname_save_field_Phase, strcat(Fname_save_field_Pathcpy_phase, Fname_save_field_Name_temp));
	pFEy_Amp = fopen(Fname_save_field_Amp, "wb");
	pFEy_Phase = fopen(Fname_save_field_Phase, "wb");

	strcpy(Fname_save_field_Pathcpy_amp, Fname_save_field_Path);
	strcpy(Fname_save_field_Pathcpy_phase, Fname_save_field_Path);
	sprintf(Fname_save_field_Name_temp, "/Ez_Amp[%1.2le].dat", sfre + ff * df);
	strcpy(Fname_save_field_Amp, strcat(Fname_save_field_Pathcpy_amp, Fname_save_field_Name_temp));
	sprintf(Fname_save_field_Name_temp, "/Ez_Phase[%1.2le].dat", sfre + ff * df);
	strcpy(Fname_save_field_Phase, strcat(Fname_save_field_Pathcpy_phase, Fname_save_field_Name_temp));
	pFEz_Amp = fopen(Fname_save_field_Amp, "wb");
	pFEz_Phase = fopen(Fname_save_field_Phase, "wb");
}

void TSaveSnapshot_expend::save_FieldDistribution()
{
	for (int pp = 0; pp < no_planes; pp++)
	{
		save_FieldDistribution_make_saveFolder(pp);
		for (int ff = 0; ff < sampling; ff++)
		{
			save_FieldDistribution_assign_file_pointer(ff);

			switch (sopd[pp].s_plane)
			{
			case 1: // xy plane
				pSave_ExAmp = Tptr1<double>(tysize); pSave_ExPhase = Tptr1<double>(tysize);
				pSave_EyAmp = Tptr1<double>(tysize); pSave_EyPhase = Tptr1<double>(tysize);
				pSave_EzAmp = Tptr1<double>(tysize); pSave_EzPhase = Tptr1<double>(tysize);

				for (int ii = 0; ii < txsize; ii++)
				{
					for (int jj = 0; jj < tysize; jj++)
					{
						pSave_ExAmp[jj] = sopds[ff][pp].AmpEx[ii][jj];
						pSave_ExPhase[jj] = sopds[ff][pp].PhaseEx[ii][jj];
						pSave_EyAmp[jj] = sopds[ff][pp].AmpEy[ii][jj];
						pSave_EyPhase[jj] = sopds[ff][pp].PhaseEy[ii][jj];
						pSave_EzAmp[jj] = sopds[ff][pp].AmpEz[ii][jj];
						pSave_EzPhase[jj] = sopds[ff][pp].PhaseEz[ii][jj];
					}
					fwrite(pSave_ExAmp, tysize, sizeof(double), pFEx_Amp);
					fwrite(pSave_ExPhase, tysize, sizeof(double), pFEx_Phase);
					fwrite(pSave_EyAmp, tysize, sizeof(double), pFEy_Amp);
					fwrite(pSave_EyPhase, tysize, sizeof(double), pFEy_Phase);
					fwrite(pSave_EzAmp, tysize, sizeof(double), pFEz_Amp);
					fwrite(pSave_EzPhase, tysize, sizeof(double), pFEz_Phase);
				}

				memfree(&pSave_ExAmp); memfree(&pSave_ExPhase);
				memfree(&pSave_EyAmp); memfree(&pSave_EyPhase);
				memfree(&pSave_EzAmp); memfree(&pSave_EzPhase);
				fclose(pFEx_Amp); fclose(pFEx_Phase);
				fclose(pFEy_Amp); fclose(pFEy_Phase);
				fclose(pFEz_Amp); fclose(pFEz_Phase);

				break;
			case 2:// yz plane
				pSave_ExAmp = Tptr1<double>(tzsize); pSave_ExPhase = Tptr1<double>(tzsize);
				pSave_EyAmp = Tptr1<double>(tzsize); pSave_EyPhase = Tptr1<double>(tzsize);
				pSave_EzAmp = Tptr1<double>(tzsize); pSave_EzPhase = Tptr1<double>(tzsize);

				for (int jj = 0; jj < tysize; jj++)
				{
					for (int kk = 0; kk < tzsize; kk++)
					{
						pSave_ExAmp[kk] = sopds[ff][pp].AmpEx[jj][kk];
						pSave_ExPhase[kk] = sopds[ff][pp].PhaseEx[jj][kk];
						pSave_EyAmp[kk] = sopds[ff][pp].AmpEy[jj][kk];
						pSave_EyPhase[kk] = sopds[ff][pp].PhaseEy[jj][kk];
						pSave_EzAmp[kk] = sopds[ff][pp].AmpEz[jj][kk];
						pSave_EzPhase[kk] = sopds[ff][pp].PhaseEz[jj][kk];
					}
					fwrite(pSave_ExAmp, tzsize, sizeof(double), pFEx_Amp);
					fwrite(pSave_ExPhase, tzsize, sizeof(double), pFEx_Phase);
					fwrite(pSave_EyAmp, tzsize, sizeof(double), pFEy_Amp);
					fwrite(pSave_EyPhase, tzsize, sizeof(double), pFEy_Phase);
					fwrite(pSave_EzAmp, tzsize, sizeof(double), pFEz_Amp);
					fwrite(pSave_EzPhase, tzsize, sizeof(double), pFEz_Phase);
				}

				memfree(&pSave_ExAmp); memfree(&pSave_ExPhase);
				memfree(&pSave_EyAmp); memfree(&pSave_EyPhase);
				memfree(&pSave_EzAmp); memfree(&pSave_EzPhase);
				fclose(pFEx_Amp); fclose(pFEx_Phase);
				fclose(pFEy_Amp); fclose(pFEy_Phase);
				fclose(pFEz_Amp); fclose(pFEz_Phase);
				break;
			case 3:// xz plane
				pSave_ExAmp = Tptr1<double>(tzsize); pSave_ExPhase = Tptr1<double>(tzsize);
				pSave_EyAmp = Tptr1<double>(tzsize); pSave_EyPhase = Tptr1<double>(tzsize);
				pSave_EzAmp = Tptr1<double>(tzsize); pSave_EzPhase = Tptr1<double>(tzsize);

				for (int ii = 0; ii < txsize; ii++)
				{
					for (int kk = 0; kk < tzsize; kk++)
					{
						pSave_ExAmp[kk] = sopds[ff][pp].AmpEx[ii][kk];
						pSave_ExPhase[kk] = sopds[ff][pp].PhaseEx[ii][kk];
						pSave_EyAmp[kk] = sopds[ff][pp].AmpEy[ii][kk];
						pSave_EyPhase[kk] = sopds[ff][pp].PhaseEy[ii][kk];
						pSave_EzAmp[kk] = sopds[ff][pp].AmpEz[ii][kk];
						pSave_EzPhase[kk] = sopds[ff][pp].PhaseEz[ii][kk];
					}
					fwrite(pSave_ExAmp, tzsize, sizeof(double), pFEx_Amp);
					fwrite(pSave_ExPhase, tzsize, sizeof(double), pFEx_Phase);
					fwrite(pSave_EyAmp, tzsize, sizeof(double), pFEy_Amp);
					fwrite(pSave_EyPhase, tzsize, sizeof(double), pFEy_Phase);
					fwrite(pSave_EzAmp, tzsize, sizeof(double), pFEz_Amp);
					fwrite(pSave_EzPhase, tzsize, sizeof(double), pFEz_Phase);
				}

				memfree(&pSave_ExAmp); memfree(&pSave_ExPhase);
				memfree(&pSave_EyAmp); memfree(&pSave_EyPhase);
				memfree(&pSave_EzAmp); memfree(&pSave_EzPhase);
				fclose(pFEx_Amp); fclose(pFEx_Phase);
				fclose(pFEy_Amp); fclose(pFEy_Phase);
				fclose(pFEz_Amp); fclose(pFEz_Phase);
				break;
			}
		}
	}
}

void TSaveSnapshot_expend::save_Current_density()
{
	double* pSave_Amp = Tptr1<double>(tysize);
	double* pSave_Phase = Tptr1<double>(tysize);
	for (int kk = 0; kk < 2; kk++)
	{
		for (int dir = 0; dir < 3; dir++)
		{
			switch (dir)
			{
			case 0:
				sprintf(Fname_Current_density_Path, "./result/Current_distribution_[%04d]_Jx", kk);
				_mkdir(Fname_Current_density_Path);
				break;
			case 1:
				sprintf(Fname_Current_density_Path, "./result/Current_distribution_[%04d]_Jy", kk);
				_mkdir(Fname_Current_density_Path);
				break;
			case 2:
				sprintf(Fname_Current_density_Path, "./result/Current_distribution_[%04d]_Jz", kk);
				_mkdir(Fname_Current_density_Path);
				break;
			}

			for (int ff = 0; ff < sampling; ff++)
			{
				strcpy(Fname_Current_density_Pathcpy, Fname_Current_density_Path);
				sprintf(Fname_Current_density_Name_temp, "/Current_distribution_Amp[%1.2le].dat", sfre + ff * df);
				strcpy(Fname_Current_density_Amp, strcat(Fname_Current_density_Pathcpy, Fname_Current_density_Name_temp));
				strcpy(Fname_Current_density_Pathcpy, Fname_Current_density_Path);
				sprintf(Fname_Current_density_Name_temp, "/Current_distribution_Phase[%1.2le].dat", sfre + ff * df);
				strcpy(Fname_Current_density_Phase, strcat(Fname_Current_density_Pathcpy, Fname_Current_density_Name_temp));
				pFname_Current_density_Amp = fopen(Fname_Current_density_Amp, "wb");
				pFname_Current_density_Phase = fopen(Fname_Current_density_Phase, "wb");
				for (int ii = 0; ii < txsize; ii++)
				{
					if (kk == 0)
						for (int jj = 0; jj < tysize; jj++)
						{
							pSave_Amp[jj] = Jup[ii][jj].Amp[dir];
							pSave_Phase[jj] = Jup[ii][jj].Phase[dir];
						}
					if (kk == 1)
						for (int jj = 0; jj < tysize; jj++)
						{
							pSave_Amp[jj] = Jdown[ii][jj].Amp[dir];
							pSave_Phase[jj] = Jdown[ii][jj].Phase[dir];
						}
					fwrite(pSave_Amp, tysize, sizeof(double), pFname_Current_density_Amp);
					fwrite(pSave_Phase, tysize, sizeof(double), pFname_Current_density_Phase);
				}
				fclose(pFname_Current_density_Amp);
				fclose(pFname_Current_density_Phase);
			}
		}
	}
}
void TSaveSnapshot_expend::Amp_Phase_Tunning()
{
	if (Amp < 0 && Phase >= 0)
	{
		Amp = -Amp;
		Phase = Phase - pi;
	}
	else if (Amp < 0 && Phase < 0)
	{
		Amp = -Amp;
		Phase = Phase + pi;
	}
}
void TSaveSnapshot_expend::save_nXh()
{
	for (int xx = 0; xx < txsize - 1; xx++)
		for (int yy = 0; yy < tysize - 1; yy++)
		{
			Jup[xx][yy].J1[0] = (Jup[xx][yy].Hy1 + Jup[xx][yy].Hzdy1 - Jup[xx][yy].Hydz1 - Jup[xx][yy].Hz1);
			Jup[xx][yy].J1[1] = (Jup[xx][yy].Hz1 + Jup[xx][yy].Hxdz1 - Jup[xx][yy].Hzdx1 - Jup[xx][yy].Hx1);
			Jup[xx][yy].J1[2] = (Jup[xx][yy].Hx1 + Jup[xx][yy].Hydx1 - Jup[xx][yy].Hxdy1 - Jup[xx][yy].Hy1);
		}
	for (int xx = 0; xx < txsize - 1; xx++)
		for (int yy = 0; yy < tysize - 1; yy++)
		{
			Jup[xx][yy].J2[0] = (Jup[xx][yy].Hy2 + Jup[xx][yy].Hzdy2 - Jup[xx][yy].Hydz2 - Jup[xx][yy].Hz2);
			Jup[xx][yy].J2[1] = (Jup[xx][yy].Hz2 + Jup[xx][yy].Hxdz2 - Jup[xx][yy].Hzdx2 - Jup[xx][yy].Hx2);
			Jup[xx][yy].J2[2] = (Jup[xx][yy].Hx2 + Jup[xx][yy].Hydx2 - Jup[xx][yy].Hxdy2 - Jup[xx][yy].Hy2);
		}
	for (int xx = 0; xx < txsize - 1; xx++)
		for (int yy = 0; yy < tysize - 1; yy++)
		{
			Jdown[xx][yy].J1[0] = (Jdown[xx][yy].Hy1 + Jdown[xx][yy].Hzdy1 - Jdown[xx][yy].Hydz1 - Jdown[xx][yy].Hz1);
			Jdown[xx][yy].J1[1] = (Jdown[xx][yy].Hz1 + Jdown[xx][yy].Hxdz1 - Jdown[xx][yy].Hzdx1 - Jdown[xx][yy].Hx1);
			Jdown[xx][yy].J1[2] = (Jdown[xx][yy].Hx1 + Jdown[xx][yy].Hydx1 - Jdown[xx][yy].Hxdy1 - Jdown[xx][yy].Hy1);
		}
	for (int xx = 0; xx < txsize - 1; xx++)
		for (int yy = 0; yy < tysize - 1; yy++)
		{
			Jdown[xx][yy].J2[0] = (Jdown[xx][yy].Hy2 + Jdown[xx][yy].Hzdy2 - Jdown[xx][yy].Hydz2 - Jdown[xx][yy].Hz2);
			Jdown[xx][yy].J2[1] = (Jdown[xx][yy].Hz2 + Jdown[xx][yy].Hxdz2 - Jdown[xx][yy].Hzdx2 - Jdown[xx][yy].Hx2);
			Jdown[xx][yy].J2[2] = (Jdown[xx][yy].Hx2 + Jdown[xx][yy].Hydx2 - Jdown[xx][yy].Hxdy2 - Jdown[xx][yy].Hy2);
		}
	//for (int xx = 0; xx < txsize; xx++)
	//	for (int yy = 0; yy < tysize; yy++)
	//		for (int zz = 0; zz < tzsize; zz++)
	//		{
	//			Jdown[xx][yy].J[0] = 0.5*(Jup[xx][yy] + HOz[xx][yy + 1][zz] - HOy[xx][yy][zz] + HOy[xx][yy][zz + 1]);
	//			Jdown[xx][yy].J[1] = 0.5*(Jup[xx][yy] + HOx[xx][yy][zz + 1] - HOz[xx][yy][zz] + HOz[xx + 1][yy][zz]);
	//			Jdown[xx][yy].J[2] = 0.5*(Jup[xx][yy] + HOy[xx + 1][yy][zz] - HOx[xx][yy][zz] + HOx[xx][yy + 1][zz]);
	//		}
}

void TSaveSnapshot_expend::set_steady_start(double steady_state_start)
{
	this->steady_state_start_discrete = round(steady_state_start / dt);
}

void TSaveSnapshot_expend::set_T0(double T0)
{
	this->T0 = T0;
	T0_discrete = round(T0 / dt);
	if (T0_discrete > 100)
		T0_stepSize = 1; // round(T0_discrete / 100.0);
	else
		T0_stepSize = 1;
}