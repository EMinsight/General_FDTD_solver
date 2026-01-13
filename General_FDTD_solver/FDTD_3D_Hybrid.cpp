<<<<<<<< HEAD:General_FDTD_solver/FDTD_Terminal.cpp
#include "FDTD_Terminal.h"
========
#include "FDTD_3D_Hybrid.h"
>>>>>>>> Hybrid_2024.0:General_FDTD_solver/FDTD_3D_Hybrid.cpp

Tmain::Tmain()
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
Tmain::~Tmain() {};

void Tmain::Initiallized_Fields()
{
	MEMORYSTATUSEX memInfo;
	memInfo.dwLength = sizeof(MEMORYSTATUSEX);
	GlobalMemoryStatusEx(&memInfo);
	DWORDLONG totalPhysMem = memInfo.ullTotalPhys;
	DWORDLONG availPhyMem = memInfo.ullAvailPhys;
	DWORDLONG physMemUsed = totalPhysMem - availPhyMem;

	if (availPhyMem < txsize*(tysize + 1)*(tzsize + 1) * sizeof(double))
	{
		printf("메모리가 부족합니다.\n");
		return;
	}
	else
		Ex = Tptr3<double>(txsize, tysize + 1, tzsize + 1);

	GlobalMemoryStatusEx(&memInfo);
	availPhyMem = memInfo.ullAvailPhys;
	physMemUsed = totalPhysMem - availPhyMem;
	if (availPhyMem < (txsize + 1)*(tysize)*(tzsize + 1) * sizeof(double))
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
	if (availPhyMem < (txsize + 1)*(tysize + 1)*(tzsize) * sizeof(double))
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
	if (availPhyMem < txsize*(tysize + 1)*(tzsize + 1) * sizeof(double))
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
	if (availPhyMem < (txsize + 1)*(tysize)*(tzsize + 1) * sizeof(double))
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
	if (availPhyMem < (txsize + 1)*(tysize + 1)*(tzsize) * sizeof(double))
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
	if (availPhyMem < (txsize + 1)*(tysize)*(tzsize) * sizeof(double))
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
	if (availPhyMem < (txsize)*(tysize + 1)*(tzsize) * sizeof(double))
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
	if (availPhyMem < (txsize)*(tysize)*(tzsize + 1) * sizeof(double))
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
}

<<<<<<<< HEAD:General_FDTD_solver/FDTD_Terminal.cpp
void Tmain::TFDTD_main_routine()
{
	// step1
	Assign_Generic_with_Get_mesh_parameters(); // Json read & update
	// step2
	Initiallized_Fields(); // Assign Ex, Ey, Ez, Jx, Jy, Jz, Hx, Hy, and Hz
	pGeneric->Field_Initiate(Ex, Ey, Ez, Jx, Jy, Jz, Hx, Hy, Hz);

	Source_infomation_upload();
	Run_Initiate();
	Time_marching();
========
void Tmain::Run_Initiate()
{
	// 5T0 ~ 7T0 , end 8T0
	NSTEPS = (int)(((1 / 1.0 / f0) * 10) / dt); // for wideband
	printf("NSTEPS = %d\n", NSTEPS);
	pGeneric->Field_Initiate(Ex, Ey, Ez, Jx, Jy, Jz, Hx, Hy, Hz);

	// pGeneric->set_PBC_Screwed_Parameters();
	media = pGeneric->get_media();
	media_Ex = pGeneric->get_media_Ex();
	media_Ey = pGeneric->get_media_Ey();
	media_Ez = pGeneric->get_media_Ez();

	// PML load
	PML1 = new TFDTD_2PML(Ex, Ey, Ez, Hx, Hy, Hz);
	PML1->PML_Initial_Setting(PMLX, PMLY, PMLZ);
	PML1->Set_Media(media_Ex, media_Ey, media_Ez);
	
	int IsPbuffer = 0;
	if (IsPbuffer == 1) // planewave - 
	{
		PbufferX = 0; 
		PbufferY = 0;
		PbufferZ = 10;
		// Plane wave
		planewave = new TFDTD_Boundary_TFSF * [number_of_source];
		for (int ii = 0; ii < number_of_source; ii++)
			planewave[ii] = new TFDTD_Boundary_TFSF(Ex, Ey, Ez, Hx, Hy, Hz);

		for (int ii = 0; ii < number_of_source; ii++)
		{
			//planewave2 = new TFDTD_Boundary_TFSF(Ex, Ey, Ez, Hx, Hy, Hz);
			int pmlwbuffer_x1, pmlwbuffer_x2, pmlwbuffer_y1, pmlwbuffer_y2, pmlwbuffer_z1, pmlwbuffer_z2;
			pmlwbuffer_x1 = sst[0].Sx1; // PMLX + PbufferX;
			pmlwbuffer_x2 = sst[0].Sx2; // txsize; //txsize-(PMLX+PbufferX);

			pmlwbuffer_y1 = sst[0].Sy1; // PMLY + PbufferY;
			pmlwbuffer_y2 = sst[0].Sy2; // tysize; //tysize-(PMLY+PbufferY); 

			pmlwbuffer_z1 = sst[0].Sz1; // 0; // PMLZ + PbufferZ;
			pmlwbuffer_z2 = sst[0].Sz2; // planewave_pos; // -(PMLZ + PbufferZ);
			
			planewave[ii]->insert_parameters(txsize, tysize, tzsize, min1D(dx), min1D(dy), dz[pmlwbuffer_z2], dt, f0, NSTEPS);
			planewave[ii]->insert_inc_angle(PW_Theta[ii], PW_Fhi[ii], PW_Polar[ii]);
			planewave[ii]->insert_TFSF_Box(pmlwbuffer_x1, pmlwbuffer_x2, pmlwbuffer_y1, pmlwbuffer_y2, pmlwbuffer_z1, pmlwbuffer_z2);

		}
	}

	// Snapshot
	nosnap = NSTEPS / 13; 
	pSnap = new TSaveSnapshot_expend(Ex, Ey, Ez, Hx, Hy, Hz);
	pSnap->Read_save_info();
	pSnap->set_TimeDiff_NstepMax(dt, NSTEPS);
	pSnap->set_spatial_and_dt( min1D(dx),  min1D(dy),  min1D(dz),  dt);
	pSnap->Set_max_iter();

	if (nosnap > 0 && nosnap <= NSTEPS)
	{
		pSnap->get_mesh_length(txsize, tysize, tzsize);
		snapshot_iter = (int)((double)NSTEPS / (double)nosnap + 0.5);	// snapshot no. iteration
		snap_index = 0;
	}
	double T0 = 1 / f0;
	pSnap->set_Media(media);
	// json -> sin, pulse
	// sine: 0.9T ~ 1T 
	// pulse: 0T ~ 1T 에 따라 아래 수치 변경
	pSnap->set_T0(T0);
	pSnap->set_steady_start(0.9 * T0);
}

void Tmain::TFDTD_main_routine()
{
	Field_data_assign_Generic(); // Json read
	Initiallized_Fields();
	Source_infomation_upload_Hybrid();
	Run_Initiate();
	Run_Hybrid();
>>>>>>>> Hybrid_2024.0:General_FDTD_solver/FDTD_3D_Hybrid.cpp

	getchar();
	return;
}

<<<<<<<< HEAD:General_FDTD_solver/FDTD_Terminal.cpp
void Tmain::Time_marching()
========
void Tmain::TFDTD_view_source() {}
void Tmain::Run_Hybrid()
>>>>>>>> Hybrid_2024.0:General_FDTD_solver/FDTD_3D_Hybrid.cpp
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
	}
	int ii, jj;
	int time1, time2;
	int nThread = omp_get_max_threads();
	
	// 아주 나중에 omp 고도화
	//int nlx = round(txsize / nThread);
	//int nly = round(tysize / nThread);
	
	int n_Thread = omp_get_num_procs();
	printf("num_proc = %d\n", n_Thread);
	n_Thread = omp_get_max_threads();
	printf("max_thread = %d\n", n_Thread);
<<<<<<<< HEAD:General_FDTD_solver/FDTD_Terminal.cpp
	omp_set_num_threads(round(n_Thread*1));

	// pSnap->set_Frequency_range(f0, sfre, tfre, sampling);
	// pOptical->Initiate_optical_option_Frequency(f0, sfre, tfre, sampling);
	for (ss = 0; ss <= NSTEPS; ss++)
	{
		//	pbProgress->Position= ss;
		//	Pulse - point
		//	time1 = clock();
		// for (int ii = 0; ii < number_of_source; ii++)
		// {
		// for (int ii = 0; ii < number_of_source; ii++)
		// 	Pulse1 = Source1->Get_Pulse(ss, Phase[ii]); // acos(-1) / 2.0);
		// 	planewave[ii]->Set_Pulse_OMP(Pulse1);
		// }
		// Pulse2 = Source1->Get_Pulse(ss, 0); 
		//planewave2->Set_Pulse_OMP(Pulse2);
		
		// Time marching 내부에 source를 업데이트를 하기 위해
		// source의 형태에 따라 edge, planewave, user define 등 으로 분배하는 cpp 파일을 추가.
		Source1->Inject_Source(ss);
========
	// update 필요: json에서 thread의 개수 혹은 비율을 받아와서 입력
	omp_set_num_threads(round(n_Thread*0.25));
	double sfre = f0 - BW / 2;
	double tfre = f0 + BW / 2;
	double sampling = 3;
	pSnap->set_Frequency_range(f0, sfre, tfre, sampling);
	for (ss = 0; ss <= NSTEPS; ss++)
	{
		pHybrid->set_Timestep(ss);
		pHybrid->Calculate_Incident_Fiend();

>>>>>>>> Hybrid_2024.0:General_FDTD_solver/FDTD_3D_Hybrid.cpp
		PML1->Fx_update();
		PML1->Fy_update();
		PML1->Fz_update();
		
		pHybrid->Apply_sourceE();
		//pHybrid->Set_In_zeroE();
		pGeneric->J_field_update();
		pGeneric->E_field_update();
		pHybrid->Set_In_zeroE();

		pGeneric->H_field_update();
		pHybrid->Apply_sourceH();
		pHybrid->Set_In_zeroH();

		PML1->Gx_update();
		PML1->Gy_update();
		PML1->Gz_update();

		pSnap->Set_Current_iter(ss);
		pSnap->Field_Obs_Points();
<<<<<<<< HEAD:General_FDTD_solver/FDTD_Terminal.cpp

		// pSnap->Save_Composition();
		// pSnap->store_snapshot_Time(ss);
		// pSnap->store_current_Time(ss);
		//pOptical->save_optical_chirality_Time(ss);
========
		if (nosnap > 0 && nosnap <= NSTEPS)
			if (((ss % snapshot_iter) == 0) || ss == NSTEPS)
			{
				pSnap->Field_Obs_Plane(snap_index);
				snap_index++;
			}
>>>>>>>> Hybrid_2024.0:General_FDTD_solver/FDTD_3D_Hybrid.cpp

		if ((ss % 100) == 0)
		{
			printf("%06.2lf[%%]\t%2.2le\t%2.2le\t%2.2le\t%2.2le\t%2.2le\t%2.2le\t%2.2le\n", ((double)ss / (double)NSTEPS) * 100, 0,
				Ex[110][60][60], Ey[110][60][60], Ez[110][60][60],
				Ex[110][60][60], Ey[110][60][60], Ez[110][60][60]);
		}
	} // *** Leap Frog Time Marching END ***//
	//pOptical->Calculate_optical_chirality_Frequency();
	//pOptical->save_optical_chirality_Frequency();

<<<<<<<< HEAD:General_FDTD_solver/FDTD_Terminal.cpp
========
	delete pHybrid;
>>>>>>>> Hybrid_2024.0:General_FDTD_solver/FDTD_3D_Hybrid.cpp
	delete PML1;
	delete pGeneric;
	delete pSnap;
}

<<<<<<<< HEAD:General_FDTD_solver/FDTD_Terminal.cpp
void Tmain::Run_Initiate()
{
	// Nstep,
	// PML, PBC
	// 컨텐츠 판단 (plate가 포함된 metasurface, cavity filter, )
	// mesh 정보1. (tx, ty, tz, dx, dy, dz, dt)
	// mesh 정보2. mediaEx, mediaEy, mediaEz
	// 
	// Boundary condition 판단
	// (PBC, screwed PBC, 
	
	// NSTEPS -> source class에서 계산후 넘겨 받기.
	NSTEPS = Source1->getNSTEPS();
	printf("NSTEPS = %d\n", NSTEPS);
	

	// pGeneric->set_PBC_Screwed_Parameters();
	media = pGeneric->get_media();
	media_Ex = pGeneric->get_media_Ex();
	media_Ey = pGeneric->get_media_Ey();
	media_Ez = pGeneric->get_media_Ez();

	// PML load
	PML1 = new TFDTD_2PML(Ex, Ey, Ez, Hx, Hy, Hz);
	PML1->PML_Initial_Setting(PMLX, PMLY, PMLZ);
	PML1->Set_Media(media_Ex, media_Ey, media_Ez);
	// pPMLPBC = new TFDTD_2PML_PBC(Ex, Ey, Ez, Hx, Hy, Hz);
	// pPMLPBC->PML_Initial_Setting(PMLX, PMLY, PMLZ);
	// pPMLPBC->Set_Kappa_coefficient(0, 1);
	// pPMLPBC->Set_Media(media_Ex, media_Ey, media_Ez);

	// Snapshot
	// NSTPES를 넘겨주고 snap class로부터 snap저장 주기 계산,
	// snap의 개수는 json에서 결정.
	// 이후 n * (Delta t)의 n값을 넘겨주고 데이터를 저장을 판단.
	nosnap = NSTEPS / 1000;
	pSnap = new TSaveSnapshot_expend(Ex, Ey, Ez, Hx, Hy, Hz);
	pSnap->Read_save_info();
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
	pSnap->set_Media(media);
	pSnap->set_T0(T0);
	pSnap->set_steady_start(0 * T0);
}

void Tmain::Source_infomation_upload()
{
	Source1 = new TSource_Excitation(Ex, Ey, Ez, Hx, Hy, Hz);
	// snap_expand class함수로, snap class에서 처리
	// 또한, f0, BW, sampling 등은 json파일에서 획득
	Source1->Set_dsdt(dx, dy, dz, dt);
========
void Tmain::Source_infomation_upload_Hybrid()
{
	Source1 = new TSource_Excitation(Ex, Ey, Ez, Hx, Hy, Hz);
	pHybrid = new TFDTD_Hybrid_mod_with_MOM(Ex, Ey, Ez, Hx, Hy, Hz);
>>>>>>>> Hybrid_2024.0:General_FDTD_solver/FDTD_3D_Hybrid.cpp
	Source1->Set_Pulse_parameters();
	
	f0 = Source1->get_f0();
	BW = Source1->get_BW();
<<<<<<<< HEAD:General_FDTD_solver/FDTD_Terminal.cpp


	
========
	Source1->Set_dt(dt);
	pHybrid->set_dt(dt);
>>>>>>>> Hybrid_2024.0:General_FDTD_solver/FDTD_3D_Hybrid.cpp
	number_of_source = Source1->Get_Number_of_Source();
	pHsi = Source1->Get_hybrid_data();
	pHybrid->insert_parameters(pHsi);
	pHybrid->set_surface_current_data();
	SourceType = Source1->Get_sourceType();
<<<<<<<< HEAD:General_FDTD_solver/FDTD_Terminal.cpp


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
========
>>>>>>>> Hybrid_2024.0:General_FDTD_solver/FDTD_3D_Hybrid.cpp
}

void Tmain::Assign_Generic_with_Get_mesh_parameters()
{
	pGeneric = new TFDTD_GenericAlgorithm();
	pGeneric->Initial_Setting();
	//pChiral->Generic_dispersive_media_Initiate();

	txsize = pGeneric->get_txsize();
	tysize = pGeneric->get_tysize();
	tzsize = pGeneric->get_tzsize();

	PMLX = pGeneric->get_PMLX();
	PMLY = pGeneric->get_PMLY();
	PMLZ = pGeneric->get_PMLZ();

	dx = pGeneric->get_dx();
	dy = pGeneric->get_dy();
	dz = pGeneric->get_dz();
	dt = pGeneric->get_dt();
	total_pol = pGeneric->get_total_pol();
}
