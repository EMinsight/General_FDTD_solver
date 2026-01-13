#include "FDTD_Terminal.h"

Tmain::Tmain()
{
	__try
	{
		_mkdir("result");
		_mkdir("result\\media");
		_mkdir("result\\optical_chirality");
		_mkdir("result\\optical_chirality_frequency");
		strcpy(JsonFileName, "./IEM_solver_input_format.json");
		eps0 = 1e-9 / 36 / m_pi;
	}
	__except (EXCEPTION_EXECUTE_HANDLER)
	{
	}
	FILE *pf_ompthread;
	pf_ompthread = fopen("Thread_num.txt", "r");
	active_thread_num = 10;// n_Thread * 0.4;
	fscanf(pf_ompthread, "%d", &active_thread_num);
	sFc.num_thread = active_thread_num;
	fclose(pf_ompthread);

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

	if (availPhyMem < txsize * (tysize + 1) * (tzsize + 1) * sizeof(double))
	{
		printf("Ex field 메모리 할당 실패.\n");
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
		printf("Ey field 메모리 할당 실패.\n");
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
		printf("Ez field 메모리 할당 실패.\n");
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
		printf("Jx field 메모리 할당 실패.\n");
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
		printf("Jy field 메모리 할당 실패.\n");
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
		printf("Jz field 메모리 할당 실패.\n");
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
		printf("Hx field 메모리 할당 실패.\n");
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
		printf("Hy field 메모리 할당 실패.\n");
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
		printf("Hz field 메모리 할당 실패.\n");
		return;
	}
	else
		Hz = Tptr3<double>(txsize, tysize, tzsize + 1);

	printf("Field memory allocation success.\n");
	sFc.Ex_address = Ex;	sFc.Ey_address = Ey;	sFc.Ez_address = Ez;
	sFc.Jx_address = Jx;	sFc.Jy_address = Jy;	sFc.Jz_address = Jz;
	sFc.Hx_address = Hx;	sFc.Hy_address = Hy;	sFc.Hz_address = Hz;
}

void Tmain::TFDTD_main_routine()
{
	// Nstep,
	// PML, PBC
	// ������ �Ǵ� (plate�� ���Ե� metasurface, cavity filter, )
	// mesh ����1. (tx, ty, tz, dx, dy, dz, dt)
	// mesh ����2. mediaEx, mediaEy, mediaEz
	// PML, PBC load (PBC�� PML�� update���� ���� �ٸ�.
	// ���� update ���� �������ֱ� ���� PBC�Ǵܺ��� ����.

	pBC = new FDTD_Boundary_Condition(&sFc);
	PBCmod = pBC->Get_BC_mode();


	if (PBCmod)
		Generic_PBC_FDTD_mode();
	else
		Generic_FDTD_mode();

	// PML �Ǵ�.
	// pBC->Set_EHJ_Fields(&sFc);
	pBC->Initicate();
	 // Assign Ex, Ey, Ez, Jx, Jy, Jz, Hx, Hy, and Hz
	

	Source_infomation_upload();
	Post_processing_info();
	if (PBCmod)
		Time_marching_PBC();
	else
		Time_marching();

	getchar();
	return;
}

void Tmain::Time_marching()
{
	double Pulse1, Pulse2;
	//pbProgress->Max = NSTEPS;
	int nth = (int)log10((double)NSTEPS);
	// int Sx0, Sy0, Sz0, Sx1, Sy1, Sz1;

	if (SourceType == 3)
	{
		// pSnap->get_Obs_point(&Sx0, &Sy0, &Sz0, &Sx1, &Sy1, &Sz1);
	}
	else
	{
		// pSnap->get_Obs_point(&Sx0, &Sy0, &Sz0, &Sx1, &Sy1, &Sz1);
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
	omp_set_dynamic(0);
	// omp_set_nested(1);
	omp_set_max_active_levels(2);
	//cout << "omp_get_max_active_levels : " << omp_get_max_active_levels() << endl;
	// omp_set_dynamic(0);
	// FILE *pf_ompthread;
	// pf_ompthread = fopen("Thread_num.txt", "r");
	// active_thread_num = 10;// n_Thread * 0.4;
	// fscanf(pf_ompthread, "%d", &active_thread_num);
	// sFc.num_thread = active_thread_num;
	// fclose(pf_ompthread);
	omp_set_num_threads(active_thread_num);
	printf("active thread : %d \n", active_thread_num);
	
	int pp = NSTEPS / 100;
	// pSnap->set_Frequency_range(f0, sfre, tfre, sampling);
	// pOptical->Initiate_optical_option_Frequency(f0, sfre, tfre, sampling);
	int snap_start = 0 * NSTEPS;
	int snap_due = NSTEPS - snap_start;
	int save_due = snap_due / nosnap;
	

	double start_time = GetTickCount64();
	for (ss = 0; ss <= NSTEPS; ss++)
	{
		//	pbProgress->Position= ss;
		//	Pulse - point
		//	time1 = clock();
		// Time marching ���ο� source�� ������Ʈ�� �ϱ� ����
		// source�� ���¿� ���� edge, planewave, user define �� ���� �й��ϴ� cpp ������ �߰�.
		
		pBC->Boundary_Condition_PMLF();

		// Apply source E
		pGeneric->J_field_update();
		pGeneric->E_field_update();
		Source1->Inject_SourceE(ss);

		// Set in zeror E
		Source1->Inject_Source_add(1);
		pGeneric->H_field_update();
		Source1->Inject_SourceH(ss);
		// Apply source H
		// Set in zeror H
		Source1->Inject_Source_add(2);
		pBC->Boundary_Condition_PMLG();
		pPost->Save(ss);
		// pSnap->Save_Composition();
		// pSnap->store_snapshot_Time(ss);
		// pSnap->store_current_Time(ss);
		// pOptical->save_optical_chirality_Time(ss);
		// pPost->Save(ss);
		// pSnap->Set_Current_iter(ss);
		// pSnap->Field_Obs_Points();
		// if (ss > snap_start > 0 && ss <= NSTEPS)
		// {
		// 	if ( (int)(ss-snap_start) % (int)(snap_due/nosnap) == 0)
		// 	{
//		// 		if (((ss % snapshot_iter) == 0) || ss == NSTEPS)
//		// 		{
		// 			pSnap->Field_Obs_Plane(snap_index);
		// 			snap_index++;
//		// 		}
		// 	}
		// }
		if ((ss % 100) == 0)
		{
			printf("%06.2lf[%%]\t%2.2le\t%2.2le\t%2.2le\t%2.2le\t%2.2le\t%2.2le\t%2.2le\n", ((double)ss / (double)NSTEPS) * 100, Pulse,
			Ex[50][50][50], Ey[50][50][50], Ez[50][50][50],
			Ex[50][50][50], Ey[50][50][50], Ez[50][50][50]);
		}
		
		//   pbProgress->Position = ss;
	} // *** Leap Frog Time Marching END ***//
	double finish_time = GetTickCount64();
	FILE* pf_time;
	char savecount[256];
	sprintf(savecount, "../count_%d.txt", active_thread_num);
	pf_time = fopen(savecount, "w");
	fprintf(pf_time,  "%f", finish_time - start_time);
	fclose(pf_time);
	printf("finished\n");
	//pOptical->Calculate_optical_chirality_Frequency();
	//pOptical->save_optical_chirality_Frequency();
	if (PBCmod)
		delete pGeneric_PBC;
	else
		delete pGeneric;
	delete pPost;
}

void Tmain::Time_marching_PBC()
{
	double Pulse1, Pulse2;
	//pbProgress->Max = NSTEPS;
	int nth = (int)log10((double)NSTEPS);
	int Sx0, Sy0, Sz0, Sx1, Sy1, Sz1;

	// if (SourceType == 3)
	// {
	// 	pSnap->get_Obs_point(&Sx0, &Sy0, &Sz0, &Sx1, &Sy1, &Sz1);
	// }
	// else
	// {
	// 	pSnap->get_Obs_point(&Sx0, &Sy0, &Sz0, &Sx1, &Sy1, &Sz1);
	// 	//Sx0 = sst[0].Sx1;
	// 	//Sy0 = sst[0].Sy1; 
	// 	//Sz0 = sst[0].Sz1; 
	// }
	// 
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
	omp_set_num_threads(round(n_Thread * 0.4));
	
	                  // nested parallelism Ȱ��ȭ

	int pp = NSTEPS / 100;
	// pSnap->set_Frequency_range(f0, sfre, tfre, sampling);
	// pOptical->Initiate_optical_option_Frequency(f0, sfre, tfre, sampling);
	int snap_start = 0 * NSTEPS;
	int snap_due = NSTEPS - snap_start;
	int save_due = snap_due / nosnap;

	for (ss = 0; ss <= NSTEPS; ss++)
	{
		//	pbProgress->Position= ss;
		//	Pulse - point
		//	time1 = clock();
		// Time marching ���ο� source�� ������Ʈ�� �ϱ� ����
		// source�� ���¿� ���� edge, planewave, user define �� ���� �й��ϴ� cpp ������ �߰�.

		pBC->Boundary_Condition_PMLF();
		// Apply source E
		Source1->Inject_SourceE(ss);
		pGeneric_PBC->J_field_update();
		pGeneric_PBC->PBC_Save_cal();
		pBC->Boundary_Condition_PBC();
		pGeneric_PBC->E_field_update();


		Source1->Inject_SourceH(ss);
		// Set in zeror E
		Source1->Inject_Source_add(1);
		pGeneric_PBC->H_field_update();
		// Apply source H
		// Set in zeror H
		Source1->Inject_Source_add(2);
		pBC->Boundary_Condition_PMLG();

		pPost->Save(ss);
		// pSnap->Save_Composition();
		// pSnap->store_snapshot_Time(ss);
		// pSnap->store_current_Time(ss);
		//pOptical->save_optical_chirality_Time(ss);

		// pSnap->Set_Current_iter(ss);
		// pSnap->Field_Obs_Points();
		// if (ss > snap_start > 0 && ss <= NSTEPS)
		// {
		// 	if ((int)(ss - snap_start) % (int)(snap_due / nosnap) == 0)
		// 	{
		// 		//				if (((ss % snapshot_iter) == 0) || ss == NSTEPS)
		// 		//				{
		// 		pSnap->Field_Obs_Plane(snap_index);
		// 		snap_index++;
		// 		//				}
		// 	}
		// }
		if ((ss % 100) == 0)
		{
			printf("%06.2lf[%%]\t%2.2le\t%2.2le\t%2.2le\t%2.2le\t%2.2le\t%2.2le\t%2.2le\n", ((double)ss / (double)NSTEPS) * 100, Pulse,
				Ex[50][50][50], Ey[50][50][50], Ez[50][50][50],
				Ex[50][50][50], Ey[50][50][50], Ez[50][50][50]);
		}
		//   pbProgress->Position = ss;
	} // *** Leap Frog Time Marching END ***//
	//pOptical->Calculate_optical_chirality_Frequency();
	//pOptical->save_optical_chirality_Frequency();
	if (PBCmod)
		delete pGeneric_PBC;
	else
		delete pGeneric;
	delete pPost;
}

void Tmain::Post_processing_info()
{
	// NSTEPS -> source class���� ����� �Ѱ� �ޱ�.
	NSTEPS = Source1->getNSTEPS();
	sFc.NSTEPS = NSTEPS;
	printf("NSTEPS = %d\n", NSTEPS);

	// Snapshot
	// NSTPES�� �Ѱ��ְ� snap class�κ��� snap���� �ֱ� ���,
	// snap�� ������ json���� ����.
	// ���� n * (Delta t)�� n���� �Ѱ��ְ� �����͸� ������ �Ǵ�.
	nosnap = 100;//NSTEPS / 30;
	// pSnap = new TSaveSnapshot_expend(Ex, Ey, Ez, Hx, Hy, Hz);
	
	pPost = new TPostprocessing(&sFc, pMp);
	pPost->Prepare2save();
	// pSnap = new TSaveSnapshot_expend(&sFc, &sMp);
	// pSnap->Read_save_info();
	// pSnap->set_TimeDiff_NstepMax(dt, NSTEPS);
	// pSnap->set_spatial_and_dt(min1D(dx), min1D(dy), min1D(dz), dt);
	// pSnap->Set_max_iter();
	// 
	// if (nosnap > 0 && nosnap <= NSTEPS)
	// {
	// 	pSnap->get_mesh_length(txsize, tysize, tzsize);
	// 	snapshot_iter = (int)((double)NSTEPS / (double)nosnap + 0.5);	// snapshot no. iteration
	// 	snap_index = 0;
	// }
	// double T0 = 1 / f0;
	// pSnap->set_Media(media);
	// //pSnap->set_Media(media_Ex, media_Ey, media_Ez);
	// pSnap->set_T0(T0);
	// pSnap->set_steady_start(0 * T0);
}

void Tmain::Generic_FDTD_mode()
{
	pGeneric = new TFDTD_GenericAlgorithm();
	// pGeneric->set_PMLX(PMLX);
	// pGeneric->set_PMLY(PMLY);
	// pGeneric->set_PMLZ(PMLZ);
	pGeneric->Initial_Setting(&sFc);
	//pChiral->Generic_dispersive_media_Initiate();
	txsize = sFc.txsize;
	tysize = sFc.tysize;
	tzsize = sFc.tzsize;
	total_pol = sFc.total_pol;
	Initiallized_Fields();

	pGeneric->Field_Initiate(&sFc);
	sFc = pGeneric->get_Field_address(&sFc);
	pMp = pGeneric->get_material_properties();
	dx = sFc.dx;
	dy = sFc.dy;
	dz = sFc.dz;
	dt = sFc.dt;
	// log ���� �ʿ�.
	
	media = sFc.media;
	media_Ex = sFc.Media_Ex;
	media_Ey = sFc.Media_Ey;
	media_Ez = sFc.Media_Ez;
	cout << total_pol << endl;
}

void Tmain::Generic_PBC_FDTD_mode()
{
	pGeneric_PBC = new TFDTD_Generic_PBC();
	// pGeneric_PBC->set_PMLX(PMLX);
	// pGeneric_PBC->set_PMLY(PMLY);
	// pGeneric_PBC->set_PMLZ(PMLZ);
	pGeneric_PBC->Initial_Setting(&sFc);
	txsize = sFc.txsize;
	tysize = sFc.tysize;
	tzsize = sFc.tzsize;
	total_pol = sFc.total_pol;
	Initiallized_Fields();

	pGeneric_PBC->Field_Initiate(&sFc);
	sFc = pGeneric_PBC->get_Field_address(&sFc);
	pMp = pGeneric_PBC->get_material_properties();

	dx = sFc.dx;
	dy = sFc.dy;
	dz = sFc.dz;
	dt = sFc.dt;
	media = sFc.media;
	media_Ex = sFc.Media_Ex;
	media_Ey = sFc.Media_Ey;
	media_Ez = sFc.Media_Ez;
}

void Tmain::Source_infomation_upload()
{
	Source1 = new TSource_Excitation(&sFc);
	// snap_expand class�Լ���, snap class���� ó��
	// ����, f0, BW, sampling ���� json���Ͽ��� ȹ��
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
