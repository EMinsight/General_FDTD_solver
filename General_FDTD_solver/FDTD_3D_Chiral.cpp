#include "FDTD_3D_Chiral.h"

TChiral::TChiral()
{
	__try
	{
		_mkdir("result");
		_mkdir("result\\media");
		_mkdir("result\\observation");
		_mkdir("result\\distribution");
		_mkdir("result\\composit");
		_mkdir("result\\optical_chirality");
		_mkdir("result\\optical_chirality_frequency");
		strcpy(JsonFileName, "./IEM_solver_input_format.json");
		eps0 = 1e-9 / 36 / m_pi;
	}
	__except (EXCEPTION_EXECUTE_HANDLER)
	{
	}

};
TChiral::~TChiral() {};

void TChiral::Initiallized_Fields()
{
	MEMORYSTATUSEX memInfo;
	memInfo.dwLength = sizeof(MEMORYSTATUSEX);
	GlobalMemoryStatusEx(&memInfo);
	DWORDLONG totalPhysMem = memInfo.ullTotalPhys;
	DWORDLONG availPhyMem = memInfo.ullAvailPhys;
	DWORDLONG physMemUsed = totalPhysMem - availPhyMem;

	if (availPhyMem < txsize * (tysize + 1) * (tzsize + 1) * sizeof(double))
	{
		printf("메모리가 부족합니다.\n");
		return;
	}
	else
		Ex = Tptr3<double>(txsize, tysize + 1, tzsize + 1);

	GlobalMemoryStatusEx(&memInfo);
	availPhyMem = memInfo.ullAvailPhys;
	physMemUsed = totalPhysMem - availPhyMem;
	if (availPhyMem < (txsize + 1) * (tysize) * (tzsize + 1) * sizeof(double))
	{
		memfree(&Ex, txsize, tysize + 1);
		printf("메모리가 부족합니다.\n");
		return;
	}
	else
		Ey = Tptr3<double>(txsize + 1, tysize, tzsize + 1);

	GlobalMemoryStatusEx(&memInfo);
	availPhyMem = memInfo.ullAvailPhys;
	physMemUsed = totalPhysMem - availPhyMem;
	if (availPhyMem < (txsize + 1) * (tysize + 1) * (tzsize) * sizeof(double))
	{
		memfree(&Ex, txsize, tysize + 1);
		memfree(&Ey, txsize + 1, tysize);
		printf("메모리가 부족합니다.\n");
		return;
	}
	else
		Ez = Tptr3<double>(txsize + 1, tysize + 1, tzsize);

	GlobalMemoryStatusEx(&memInfo);
	availPhyMem = memInfo.ullAvailPhys;
	physMemUsed = totalPhysMem - availPhyMem;
	if (availPhyMem < txsize * (tysize + 1) * (tzsize + 1) * sizeof(double))
	{
		memfree(&Ex, txsize, tysize + 1);
		memfree(&Ey, txsize + 1, tysize);
		memfree(&Ez, txsize + 1, tysize + 1);
		printf("메모리가 부족합니다.\n");
		return;
	}
	else
		Jx = Tptr4<double>(txsize, tysize + 1, tzsize + 1, total_pol);

	GlobalMemoryStatusEx(&memInfo);
	availPhyMem = memInfo.ullAvailPhys;
	physMemUsed = totalPhysMem - availPhyMem;
	if (availPhyMem < (txsize + 1) * (tysize) * (tzsize + 1) * sizeof(double))
	{
		memfree(&Ex, txsize, tysize + 1);
		memfree(&Ey, txsize + 1, tysize);
		memfree(&Ez, txsize + 1, tysize + 1);
		memfree(&Jx, txsize, tysize + 1, tzsize);
		printf("메모리가 부족합니다.\n");
		return;
	}
	else
		Jy = Tptr4<double>(txsize + 1, tysize, tzsize + 1, total_pol);

	GlobalMemoryStatusEx(&memInfo);
	availPhyMem = memInfo.ullAvailPhys;
	physMemUsed = totalPhysMem - availPhyMem;
	if (availPhyMem < (txsize + 1) * (tysize + 1) * (tzsize) * sizeof(double))
	{
		memfree(&Ex, txsize, tysize + 1);
		memfree(&Ey, txsize + 1, tysize);
		memfree(&Ez, txsize + 1, tysize + 1);
		memfree(&Jx, txsize, tysize + 1, tzsize);
		memfree(&Jy, txsize + 1, tysize, tzsize);
		printf("메모리가 부족합니다.\n");
		return;
	}
	else
		Jz = Tptr4<double>(txsize + 1, tysize + 1, tzsize, total_pol);

	GlobalMemoryStatusEx(&memInfo);
	availPhyMem = memInfo.ullAvailPhys;
	physMemUsed = totalPhysMem - availPhyMem;
	if (availPhyMem < (txsize + 1) * (tysize) * (tzsize) * sizeof(double))
	{
		memfree(&Ex, txsize, tysize + 1);
		memfree(&Ey, txsize + 1, tysize);
		memfree(&Ez, txsize + 1, tysize + 1);
		memfree(&Jx, txsize, tysize + 1, tzsize);
		memfree(&Jy, txsize + 1, tysize, tzsize);
		memfree(&Jz, txsize + 1, tysize + 1, tzsize);
		printf("메모리가 부족합니다.\n");
		return;
	}
	else
		Hx = Tptr3<double>(txsize + 1, tysize, tzsize);

	GlobalMemoryStatusEx(&memInfo);
	availPhyMem = memInfo.ullAvailPhys;
	physMemUsed = totalPhysMem - availPhyMem;
	if (availPhyMem < (txsize) * (tysize + 1) * (tzsize) * sizeof(double))
	{
		memfree(&Ex, txsize, tysize + 1);
		memfree(&Ey, txsize + 1, tysize);
		memfree(&Ez, txsize + 1, tysize + 1);
		memfree(&Jx, txsize, tysize + 1, tzsize);
		memfree(&Jy, txsize + 1, tysize, tzsize);
		memfree(&Jz, txsize + 1, tysize + 1, tzsize);
		memfree(&Hx, txsize + 1, tysize);
		printf("메모리가 부족합니다.\n");
		return;
	}
	else
		Hy = Tptr3<double>(txsize, tysize + 1, tzsize);

	GlobalMemoryStatusEx(&memInfo);
	availPhyMem = memInfo.ullAvailPhys;
	physMemUsed = totalPhysMem - availPhyMem;
	if (availPhyMem < (txsize) * (tysize) * (tzsize + 1) * sizeof(double))
	{
		memfree(&Ex, txsize, tysize + 1);
		memfree(&Ey, txsize + 1, tysize);
		memfree(&Ez, txsize + 1, tysize + 1);
		memfree(&Jx, txsize, tysize + 1, tzsize);
		memfree(&Jy, txsize + 1, tysize, tzsize);
		memfree(&Jz, txsize + 1, tysize + 1, tzsize);
		memfree(&Hx, txsize + 1, tysize);
		memfree(&Hx, txsize, tysize + 1);
		printf("메모리가 부족합니다.\n");
		return;
	}
	else
		Hz = Tptr3<double>(txsize, tysize, tzsize + 1);

	sFc.Ex_address = Ex;
	sFc.Ey_address = Ey;
	sFc.Ez_address = Ez;
	sFc.Hx_address = Hx;
	sFc.Hy_address = Hy;
	sFc.Hz_address = Hz;

}

void TChiral::Run_Initiate()
{
	// 5T0 ~ 7T0 , end 8T0
	NSTEPS = Source1->getNSTEPS();
	printf("NSTEPS = %d\n", NSTEPS);
	pChiral->Field_Initiate(Ex, Ey, Ez, Jx, Jy, Jz, Hx, Hy, Hz);
	media = pChiral->get_media();
	media_Ex = pChiral->get_media_Ex();
	media_Ey = pChiral->get_media_Ey();
	media_Ez = pChiral->get_media_Ez();

	pBC = new FDTD_Boundary_Condition(Ex, Ey, Ez, Hx, Hy, Hz);
	
	// Snapshot
	nosnap = NSTEPS / 1000;
	pSnap = new TSaveSnapshot_expend(Ex, Ey, Ez, Hx, Hy, Hz);
	pOptical = new Toptical_chirality(Ex, Ey, Ez, Hx, Hy, Hz);

	pSnap->Read_save_info();
	pOptical->Read_save_info();
	pSnap->set_TimeDiff_NstepMax(dt, NSTEPS);
	pSnap->set_spatial_and_dt(min1D(dx), min1D(dy), min1D(dz), dt);
	pSnap->Set_max_iter();

	if (nosnap > 0 && nosnap <= NSTEPS)
	{
		pSnap->get_mesh_length(txsize, tysize, tzsize);
		snapshot_iter = (int)((double)NSTEPS / (double)nosnap + 0.5);	// snapshot no. iteration
		snap_index = 0;
	}
	double T0 = 1 / f0;
	pOptical->get_mesh_length(txsize, tysize, tzsize);
	pOptical->set_spatial_and_dt(min1D(dx), min1D(dy), min1D(dz), dt);
	pSnap->set_Media(media);
	pOptical->set_Media(media);

	//pOptical->Exctract_chiral_surface();
	pSnap->set_T0(T0);
	pOptical->set_T0(T0);
	// test할때는 1T ~ 1.2T 로 충분.
	pSnap->set_steady_start(3 * T0);
	pOptical->set_steady_start(3 * T0);

}

void TChiral::TFDTD_main_routine()
{
	Field_data_assign_Generic();
	Initiallized_Fields();
	Source_infomation_upload();
	Run_Initiate();
	Run_Chiral();

	getchar();
	return;
}

void TChiral::Run_Chiral()
{
	double Pulse1, Pulse2;
	//pbProgress->Max = NSTEPS;
	int nth = (int)log10((double)NSTEPS);
	int Sx0, Sy0, Sz0, Sx1, Sy1, Sz1;

	if (SourceType == 3)
	{
		pSnap->get_Obs_point(&Sx0, &Sy0, &Sz0, &Sx1, &Sy1, &Sz1);
	}
	else
	{
		pSnap->get_Obs_point(&Sx0, &Sy0, &Sz0, &Sx1, &Sy1, &Sz1);
		//Sx0 = sst[0].Sx1;
		//Sy0 = sst[0].Sy1; 
		//Sz0 = sst[0].Sz1; 
	}

	// pSnap->get_Obs_point(&Sx0, &Sy0, &Sz0, &Sx1, &Sy1, &Sz1);
	// int ***medi = pGeneric->get_media();
	// pSnap->set_Media(medi);
	// pSnap->Initiate_optical_option();
	int ii, jj;
	int time1, time2;
	int nThread = omp_get_max_threads();
	//int nlx = round(txsize / nThread);
	//int nly = round(tysize / nThread);

	int n_Thread = omp_get_num_procs();
	printf("num_proc = %d\n", n_Thread);
	n_Thread = omp_get_max_threads();
	printf("max_thread = %d\n", n_Thread);
	omp_set_num_threads(round(n_Thread * 0.25));
	BW = 0;
	double sfre = f0 - BW / 10; // 640 nm
	double tfre = f0 + BW / 10;
	double sampling = 1;
	pSnap->set_Frequency_range(f0, sfre, tfre, sampling);
	pOptical->Initiate_optical_option_Frequency(f0, sfre, tfre, sampling);
	for (ss = 0; ss <= NSTEPS; ss++)
	{
		//	pbProgress->Position= ss;
		//	Pulse - point
		//	time1 = clock();
		pBC->Boundary_Condition_PMLF(); // PML
		Source1->Inject_SourceE(ss);

		pChiral->PBC_J_field_update();
		pChiral->FDTD_PBC_Chiral_Save_cal();
		
		pBC->Boundary_Condition_PBC();  // PBC
		
		pChiral->FDTD_PBC_Chiral_E_cal();
		Source1->Inject_SourceH(ss);
		Source1->Inject_Source_add(1);

		pChiral->FDTD_PBC_Chiral_H_cal();
		Source1->Inject_Source_add(2);

		pBC->Boundary_Condition_PMLG();


		pSnap->Set_Current_iter(ss);
		pSnap->Field_Obs_Points();
		pSnap->Save_Composition();

		//pSnap->store_snapshot_Time(ss);
		pSnap->store_current_Time(ss);
		pSnap->store_snapshot_Time(ss);
		pOptical->save_optical_chirality_Time(ss);
		// pSnap
		// 10T save snapshot
		// pSnap->Field_Obs_Plane(ss);
		if (ss % 1000 == 0)
			pSnap->Field_Obs_Plane(ss);
		else if (ss % (int)(0.01 * NSTEPS) == 0)
		{
			printf("%06.2lf[%%]\t%2.2le\t%2.2le\t%2.2le\t%2.2le\t%2.2le\t%2.2le\t%2.2le\n", ((double)ss / (double)NSTEPS) * 100, Pulse1,
				Ex[Sx0][Sy0][Sz0], Ey[Sx0][Sy0][Sz0], Ez[Sx0][Sy0][Sz0],
				Ex[Sx1][Sy1][Sz1], Ey[Sx1][Sy1][Sz1], Ez[Sx1][Sy1][Sz1]);
		}
		// if (ss == NSTEPS - 1)
		// {
		// 	Peroid = Source1->get_Peroid();
		// 	NSTEPS = (int)(((1 / 1.0 / f0) * Peroid) / dt); // for wideband
		// }

		//   pbProgress->Position = ss;
	} // *** Leap Frog Time Marching END ***//
	pOptical->Calculate_optical_chirality_Frequency();
	pOptical->save_optical_chirality_Frequency();

	delete pGeneric;
	delete pSnap;
}

void TChiral::Source_infomation_upload()
{
	Source1 = new TSource_Excitation(Ex, Ey, Ez, Hx, Hy, Hz);
	// snap_expand class함수로, snap class에서 처리
	// 또한, f0, BW, sampling 등은 json파일에서 획득
	Source1->Set_dsdt(dx, dy, dz, dt);
	Source1->Set_Pulse_parameters();

	f0 = Source1->get_f0();
	BW = Source1->get_BW();

	number_of_source = Source1->Get_Number_of_Source();
	sst = Source1->Get_source_data();
	SourceType = Source1->Get_sourceType();

	switch (SourceType)
	{
	case 0:
		printf("Resistive Volate Source");
		for (int i = 0; i < number_of_source; i++)
		{
			printf("%d\t%d\t%d\n", sst[i].Sx1, sst[i].Sy1, sst[i].Sz1);
			printf("%d\t%d\t%d\n", sst[i].Sx2, sst[i].Sy2, sst[i].Sz2);
		}
		break;
	case 1:
		printf("Current Source");
		for (int i = 0; i < number_of_source; i++)
		{
			printf("%d\t%d\t%d\n", sst[i].Sx1, sst[i].Sy1, sst[i].Sz1);
			printf("%d\t%d\t%d\n", sst[i].Sx2, sst[i].Sy2, sst[i].Sz2);
		}
		break;
	case 2:
		printf("Point Source");
		for (int i = 0; i < number_of_source; i++)
		{
			printf("%d\t%d\t%d\n", sst[i].Sx1, sst[i].Sy1, sst[i].Sz1);
			printf("%d\t%d\t%d\n", sst[i].Sx2, sst[i].Sy2, sst[i].Sz2);
		}
		break;
	case 3:
		printf("Planewave Source\n");
		break;
	}
}

void TChiral::Field_data_assign_Generic()
{
	pChiral = new TFDTD_Chiral_Contents();
	pChiral->Initial_Setting();
	//pChiral->Generic_dispersive_media_Initiate();

	txsize = pChiral->get_txsize();
	tysize = pChiral->get_tysize();
	tzsize = pChiral->get_tzsize();

	dx = pChiral->get_dx();
	dy = pChiral->get_dy();
	dz = pChiral->get_dz();
	dt = pChiral->get_dt();
	total_pol = pChiral->get_total_pol();
}
