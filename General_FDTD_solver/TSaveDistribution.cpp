#include "TSaveDistribution.h"

TSaveDistribution::TSaveDistribution()
{
	_mkdir("result\\distribution");
	snap_count = 0;
}

void TSaveDistribution::Prepare2save()
{
	// harmonic과 Pulse저장 방법 다르게.
	Read_save_info();
	Read_save_mod();
}

void TSaveDistribution::Save(int Nststep_iter)
{
	// Gaussian modulated sinewave -> full save mod
	// Harmonic -> steady-state save mod
	
	if (save_sta <= Nststep_iter && save_end >= Nststep_iter && strcmp(Domain, "Time") == 0)
		Save_POI(Nststep_iter);
	else if ((save_sta == Nststep_iter || save_end == Nststep_iter) && strcmp(Domain, "Frequency") == 0)
		Save_POI_frequency(Nststep_iter);
}

void TSaveDistribution::Save_POI(int Nststep_iter)
{
	if (Nststep_iter == save_sta + snap_count * snap_iter)
	{
		for (int pp = 0; pp < total_no_planes; pp++)
		{
			switch (sopd[pp].s_plane)
			{
			case 1: //  xy 평면
				// Ex, Ey, Ez norm
				XY_plane_plot(pp, Plane_save_name[pp]);
				// Ex, Ey, Ez each
				Field_XY_2D_Each(pp, snap_count, Plane_save_name[pp]);
				break;

			case 2: //  yz 평면
				// Ex, Ey, Ez norm
				YZ_plane_plot(pp, Plane_save_name[pp]);
				// Ex, Ey, Ez each
				Field_YZ_2D_Each(pp, snap_count, Plane_save_name[pp]);
				break;

			case 3: //  xz 평면
				// Ex, Ey, Ez norm
				XZ_plane_plot(pp, Plane_save_name[pp]);
				// Ex, Ey, Ez each
				Field_XZ_2D_Each(pp, snap_count, Plane_save_name[pp]);
				break;
				// default: // 굴곡.
				// 	sprintf(Plane_save_name[pp], "./result/distribution/XY_Plane[z=%2.2e]/Field_distribution_%04d.dat", sopd[ii].r_th_real, snap_count);

				break;
			}
		}
		snap_count++;
	}
}

void TSaveDistribution::Read_save_mod()
{
	// steady-sate 는 0.85 Steps ~ 0.95 Steps, time으로 저장, (frequency 저장은 추후 업데이트)
	// steady-state 구간은 추후 json에서 읽어들어올 예정.
	ifstream ifs(JsonFileName);
	IStreamWrapper isw(ifs);
	Jdocument.ParseStream(isw);
	double steady_points = (0.01 / pFc->f0 / pFc->dt)/(double)(pFc->NSTEPS);
	no_snap = round(1 / pFc->f0 / pFc->dt); //PulseName, Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["PulseName"].GetInt()); // json에서 개수 입력 필요.

	// saveMod Harmonic, Frequency 선택하도록 json 업데이트
	if (strcmp(PulseName, "HARMONIC") == 0)
	{
		// snapshot 개수를 받아 특정 주기마다 snap저장하는 알고리즘 필요.
		// set_steady_start(0.85);
		// set_steady_end(0.85 + steady_points);
		// strcpy(Domain, "Frequency");
		set_steady_start(0);
		set_steady_end(1);
		strcpy(Domain, "Time");
		sampling = 1;
		pFc->f0;
	}
	else if (strcmp(PulseName, "GAUSSIAN MODULATED SINEWAVE") == 0)
	{
		set_steady_start(0);
		set_steady_end(1);
		strcpy(Domain, "Time");
	}
	if ((save_end - save_sta) < (no_snap - 1))
		snap_iter = 1;
	else
		snap_iter = (save_end - save_sta) / (no_snap - 1);
	ifs.close();
}

void TSaveDistribution::Read_save_info()
{
	ifstream ifs(JsonFileName);
	IStreamWrapper isw(ifs);
	Jdocument.ParseStream(isw);

	// 관찰평면의 갯수 및 관찰평면 지정
	//fscanf(pFdtd, "%d\n", &no_planes);
	no_planes = Jdocument["SIMULATION_ENVIORONMENT"]["Distribution"].Size();
	strcpy(PulseName, Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["PulseName"].GetString());

	double* temp_coordx, * temp_coordy, * temp_coordz;
	double maxx, maxy, maxz;

	maxx = txsize + 1;
	maxy = tysize + 1;
	maxz = tzsize + 1;

	if (strcmp(PulseName, "HARMONIC") == 0)
	{
		int noFlatnum = 0;
		if( Jdocument["SIMULATION_ENVIORONMENT"].HasMember("Distribution_notFlat"))
			noFlatnum = Jdocument["SIMULATION_ENVIORONMENT"]["Distribution_notFlat"].Size();
		int*** pEnvelop_Material = Tptr3<int>(2, txsize, tysize);

		ActualNum = 0;
		psEnvelop_Material = Tptr1<struct_no_planes_envelop>(1);
		struct_no_planes_envelop* psEnvelop_Material_temp;
		for (int ii = 0; ii < noFlatnum; ii++)
		{
			strcpy(matName, Jdocument["SIMULATION_ENVIORONMENT"]["Distribution_notFlat"][ii]["Materials"].GetString());
			// ex.)matName == Gold, PDMS, Copper...
			flag = TEnvelop_Extract::Chiral_Envelop_Extract(pFc, pMp, pEnvelop_Material, matName);
			// step1. 물질이 존재하면, 
			if (flag)
			{
				// step2. psEnvelop_Material의 크기로 psEnvelop_Material_temp생성
				// step3. 현재 psEnvelop_Material를 psEnvelop_Material_temp에 복제.
					// step 3-1. psEnvelop_Material 할당 해제
					// step 3-2. psEnvelop_Materialtemp -> psEnvelop_Material
				// step4. psEnvelop_Material 초기화.
				// step5. Chiral_Envelop_Extract 결과를 psEnvelop_Material에 입력.
				// step6. psEnvelop_Material_temp 를 할당 해제 (**중요** 내부 동적할당 객체들 전부  해제 해야함.)

				ActualNum++;
				// step2
				psEnvelop_Material_temp = Tptr1<struct_no_planes_envelop>(ActualNum);
				for (int jj = 0; jj < ActualNum; jj++)
					psEnvelop_Material_temp[jj].pEnvelop = Tptr3<int>(2, txsize, tysize);
				// step3, 현재 psEnvelop_Material를 psEnvelop_Material_temp에 복제. (처음은 ActualNum-1 ==0 으로 실행 x)
				for (int jj = 0; jj < ActualNum - 1; jj++)
				{
					strcpy(psEnvelop_Material_temp[jj].Name, psEnvelop_Material[jj].Name);
					for (int nn = 0; nn < 2; nn++)
						for (int xx = 0; xx < txsize; xx++)
							memcpy(psEnvelop_Material_temp[jj].pEnvelop[nn][xx], psEnvelop_Material[jj].pEnvelop[nn][xx], tysize * sizeof(int));
					// step 3-1 할당 해제
					memfree(&psEnvelop_Material[jj].pEnvelop, 2, txsize);
				}
				memfree(&psEnvelop_Material);
				// step 4. psEnvelop_Material 초기화.
				psEnvelop_Material = Tptr1<struct_no_planes_envelop>(ActualNum);
				for (int jj = 0; jj < ActualNum; jj++)
					psEnvelop_Material[jj].pEnvelop = Tptr3<int>(2, txsize, tysize);
				// step 4-1. 이전 temp 데이터 가져오기.
				for (int jj = 0; jj < ActualNum - 1; jj++)
				{
					strcpy(psEnvelop_Material[jj].Name, psEnvelop_Material_temp[jj].Name);
					for (int nn = 0; nn < 2; nn++)
						for (int xx = 0; xx < txsize; xx++)
							memcpy(psEnvelop_Material[jj].pEnvelop[nn][xx], psEnvelop_Material_temp[jj].pEnvelop[nn][xx], tysize * sizeof(int));
				}

				// step5. Chiral_Envelop_Extract 결과를 psEnvelop_Material에 입력.
				strcpy(psEnvelop_Material[ActualNum - 1].Name, matName);
				for (int nn = 0; nn < 2; nn++)
					for (int xx = 0; xx < txsize; xx++)
						// indexing은 0 ~ A-1
						memcpy(psEnvelop_Material[ActualNum - 1].pEnvelop[nn][xx], pEnvelop_Material[nn][xx], tysize * sizeof(int));

				// step6. psEnvelop_Material_temp 를 할당 해제 (**중요** 내부 동적할당 객체들 전부  해제 해야함.)
				for (int jj = 0; jj < ActualNum; jj++)
					memfree(&psEnvelop_Material_temp[jj].pEnvelop, 2, txsize);
				memfree(&psEnvelop_Material_temp);
			}
			else
			{
				// step1.1. 존재하지 않으면 continue.
				ActualNum = ActualNum;
				continue;
			}
		}
		if (noFlatnum == 0)
		{
			memfree(&pEnvelop_Material, 2, txsize);
			memfree(&psEnvelop_Material);
		}

		// matName에 저장된 만큼  
		if (ActualNum == 0)
		{
			total_no_planes = no_planes;
			sopd = Tptr1<struct_Obs_Plane>(no_planes);
			Plane_save_name = Tptr2<char>(no_planes, 256);
		}
		else
		{
			total_no_planes = no_planes + 2 * ActualNum;
			sopd = Tptr1<struct_Obs_Plane>(total_no_planes);
			Plane_save_name = Tptr2<char>(total_no_planes, 256);
			for (int ii = 0; ii < ActualNum; ii++)
			{
				sopd[no_planes + 2 * ii].s_plane = 1;
				sopd[no_planes + 2 * ii + 1].s_plane = 1;
				sopd[no_planes + 2 * ii].is_sopdf = true;
				sopd[no_planes + 2 * ii + 1].is_sopdf = true;
			}
		}
	}
	else // Gaussian
	{
		total_no_planes = no_planes;
		sopd = Tptr1<struct_Obs_Plane>(total_no_planes);
		Plane_save_name = Tptr2<char>(total_no_planes, 256);
	}
	// Obs_save_name = Tptr1<string>(no_planes);
	temp_coordx = Tptr1<double>(maxx);
	temp_coordy = Tptr1<double>(maxy);
	temp_coordz = Tptr1<double>(maxz);
	for (int ii = 0; ii < maxx; ii++)
		temp_coordx[ii] = Jdocument["MESH"]["CellLocation"]["xcoordinate"][ii].GetDouble();

	for (int ii = 0; ii < maxy; ii++)
		temp_coordy[ii] = Jdocument["MESH"]["CellLocation"]["ycoordinate"][ii].GetDouble();

	for (int ii = 0; ii < maxz; ii++)
		temp_coordz[ii] = Jdocument["MESH"]["CellLocation"]["zcoordinate"][ii].GetDouble();

	char plane[3];
	int splane = 0;
	for (int nn = 0; nn < no_planes; nn++)
	{
		//fscanf(pFdtd, "%d,%d\n", &sopd[ii].s_plane, &sopd[ii].r_th);
		strcpy(plane, Jdocument["SIMULATION_ENVIORONMENT"]["Distribution"][nn]["Plane"].GetString());
		sopd[nn].r_th_real = Jdocument["SIMULATION_ENVIORONMENT"]["Distribution"][nn]["Location"].GetDouble();
		if (strcmp(plane, "XY") == 0)
		{
			sopd[nn].s_plane = 1;
			for (int ii = 0; ii < maxz; ii++)
				if (temp_coordz[ii] == sopd[nn].r_th_real)
					sopd[nn].r_th = ii;

			sprintf(Plane_save_name[nn], "./result/distribution/XY_Plane[z=%2.2e]", sopd[nn].r_th_real);
			_mkdir(Plane_save_name[nn]);
		}
		if (strcmp(plane, "YZ") == 0)
		{
			sopd[nn].s_plane = 2;
			for (int ii = 0; ii < maxx; ii++)
				if (temp_coordx[ii] == sopd[nn].r_th_real)
					sopd[nn].r_th = ii;

			sprintf(Plane_save_name[nn], "./result/distribution/YZ_Plane[x=%2.2e]", sopd[nn].r_th_real);
			_mkdir(Plane_save_name[nn]);
		}
		if (strcmp(plane, "XZ") == 0)
		{
			sopd[nn].s_plane = 3;
			for (int ii = 0; ii < maxy; ii++)
				if (temp_coordy[ii] == sopd[nn].r_th_real)
					sopd[nn].r_th = ii;

			sprintf(Plane_save_name[nn], "./result/distribution/XZ_Plane[y=%2.2e]", sopd[nn].r_th_real);
			_mkdir(Plane_save_name[nn]);
		}
	}

	int nn = 0;
	for (int nn1 = 0; nn1 < ActualNum; nn1++)
	{
		for (int nn2 = 0; nn2 < 2; nn2++)
		{
			nn = nn1 * 2 + nn2;
			if (nn2 == 0)
				sprintf(Plane_save_name[no_planes + nn], "./result/distribution/XY_Plane_Envelop[%s]_up", psEnvelop_Material[nn1].Name);
			else if (nn2 == 1)
				sprintf(Plane_save_name[no_planes + nn], "./result/distribution/XY_Plane_Envelop[%s]_down", psEnvelop_Material[nn1].Name);
			_mkdir(Plane_save_name[no_planes + nn]);
		}
	}
	if (strcmp(PulseName, "HARMONIC") == 0)
		MemoryAssign_for_FieldDistribution();
	for (int ii = 0; ii < total_no_planes; ii++)
		cout << Plane_save_name[ii] << endl;

	memfree(&temp_coordx);
	memfree(&temp_coordy);
	memfree(&temp_coordz);
	ifs.close();
}

// time domain
// Ex, Ey, Ez norm
void TSaveDistribution::XY_plane_plot(int kk, char *Fname)
{
	sprintf(pFileName, "/Field_distribution_%04d.dat", snap_count);
	strcpy(pFileName_tmp, Fname);
	strcat(pFileName_tmp, pFileName);
	pf_Plane = fopen(pFileName_tmp, "wb");   // append, write, field point.

	E_value_2D = Tptr2<double>(txsize, tysize);
	for (int ii = 0; ii < txsize; ii++)
	{
		for (int jj = 0; jj < tysize; jj++)
		{
			ExV = (EOx[ii][jj][sopd[kk].r_th] +
				EOx[ii][jj + 1][sopd[kk].r_th] +
				EOx[ii][jj][sopd[kk].r_th + 1] +
				EOx[ii][jj + 1][sopd[kk].r_th + 1]) / 4;
			EyV = (EOy[ii][jj][sopd[kk].r_th] +
				EOy[ii + 1][jj][sopd[kk].r_th] +
				EOy[ii][jj][sopd[kk].r_th + 1] +
				EOy[ii + 1][jj][sopd[kk].r_th + 1]) / 4;
			EzV = (EOz[ii][jj][sopd[kk].r_th] +
				EOz[ii + 1][jj][sopd[kk].r_th] +
				EOz[ii][jj + 1][sopd[kk].r_th] +
				EOz[ii + 1][jj + 1][sopd[kk].r_th]) / 4;

			HxV = (HOx[ii][jj][sopd[kk].r_th] +
				HOx[ii + 1][jj][sopd[kk].r_th]) / 2;
			HyV = (HOy[ii][jj][sopd[kk].r_th] +
				HOy[ii][jj + 1][sopd[kk].r_th]) / 2;
			HzV = (HOz[ii][jj][sopd[kk].r_th] +
				HOz[ii][jj][sopd[kk].r_th + 1]) / 2;
			E_value_2D[ii][jj] = sqrt(ExV * ExV + EyV * EyV + EzV * EzV);

		}
		fwrite(E_value_2D[ii], sizeof(double), tysize, pf_Plane);

	}
	memfree(&E_value_2D, txsize);
	fclose(pf_Plane);
}
void TSaveDistribution::YZ_plane_plot(int ii, char* Fname)
{
	sprintf(pFileName, "/Field_distribution_%04d.dat", snap_count);
	strcpy(pFileName_tmp, Fname);
	strcat(pFileName_tmp, pFileName);
	pf_Plane = fopen(pFileName_tmp, "wb");   // append, write, field point.
	//  yz 평면
	E_value_2D = Tptr2<double>(tysize, tzsize);
	for (int jj = 0; jj < tysize; jj++)
	{
		for (int kk = 0; kk < tzsize; kk++)
		{

			ExV = (EOx[sopd[ii].r_th][jj][kk] +
				EOx[sopd[ii].r_th][jj + 1][kk] +
				EOx[sopd[ii].r_th][jj][kk + 1] +
				EOx[sopd[ii].r_th][jj + 1][kk + 1]) / 4;
			EyV = (EOy[sopd[ii].r_th][jj][kk] +
				EOy[sopd[ii].r_th + 1][jj][kk] +
				EOy[sopd[ii].r_th][jj][kk + 1] +
				EOy[sopd[ii].r_th + 1][jj][kk + 1]) / 4;
			EzV = (EOz[sopd[ii].r_th][jj][kk] +
				EOz[sopd[ii].r_th + 1][jj][kk] +
				EOz[sopd[ii].r_th][jj + 1][kk] +
				EOz[sopd[ii].r_th + 1][jj + 1][kk]) / 4;

			HxV = (HOx[sopd[ii].r_th][jj][kk] +
				HOx[sopd[ii].r_th + 1][jj][kk]) / 2;
			HyV = (HOy[sopd[ii].r_th][jj][kk] +
				HOy[sopd[ii].r_th][jj + 1][kk]) / 2;
			HzV = (HOz[sopd[ii].r_th][jj][kk] +
				HOz[sopd[ii].r_th][jj][kk + 1]) / 2;

			E_value_2D[jj][kk] = sqrt(ExV * ExV + EyV * EyV + EzV * EzV);

		}
		fwrite(E_value_2D[jj], sizeof(double), tzsize, pf_Plane);

	}
	memfree(&E_value_2D, tysize);
	fclose(pf_Plane);
}
void TSaveDistribution::XZ_plane_plot(int jj, char* Fname)
{	
	sprintf(pFileName, "/Field_distribution_%04d.dat", snap_count);
	strcpy(pFileName_tmp, Fname);
	strcat(pFileName_tmp, pFileName);
	pf_Plane = fopen(pFileName_tmp, "wb");   // append, write, field point.
	//  zx 평면
	E_value_2D = Tptr2<double>(txsize, tzsize);
	for (int ii = 0; ii < txsize; ii++)
	{
		for (int kk = 0; kk < tzsize; kk++)
		{
			ExV = (EOx[ii][sopd[jj].r_th][kk] +
				EOx[ii][sopd[jj].r_th + 1][kk] +
				EOx[ii][sopd[jj].r_th][kk + 1] +
				EOx[ii][sopd[jj].r_th + 1][kk + 1]) / 4;
			EyV = (EOy[ii][sopd[jj].r_th][kk] +
				EOy[ii + 1][sopd[jj].r_th][kk] +
				EOy[ii][sopd[jj].r_th][kk + 1] +
				EOy[ii + 1][sopd[jj].r_th][kk + 1]) / 4;
			EzV = (EOz[ii][sopd[jj].r_th][kk] +
				EOz[ii + 1][sopd[jj].r_th][kk] +
				EOz[ii][sopd[jj].r_th + 1][kk] +
				EOz[ii + 1][sopd[jj].r_th + 1][kk]) / 4;

			HxV = (HOx[ii][sopd[jj].r_th][kk] +
				HOx[ii + 1][sopd[jj].r_th][kk]) / 2;
			HyV = (HOy[ii][sopd[jj].r_th][kk] +
				HOy[ii][sopd[jj].r_th + 1][kk]) / 2;
			HzV = (HOz[ii][sopd[jj].r_th][kk] +
				HOz[ii][sopd[jj].r_th][kk + 1]) / 2;

			E_value_2D[ii][kk] = sqrt(ExV * ExV + EyV * EyV + EzV * EzV);

		}
		fwrite(E_value_2D[ii], sizeof(double), tzsize, pf_Plane);

	}
	memfree(&E_value_2D, txsize);
	fclose(pf_Plane);
}

// Ex, Ey, Ez each
// XY plane
void TSaveDistribution::Field_XY_2D_Each(int kk, int ss, char* Fname)
{
	double** E_value_2D_Ex = Tptr2<double>(txsize, tysize);
	double** E_value_2D_Ey = Tptr2<double>(txsize, tysize);
	double** E_value_2D_Ez = Tptr2<double>(txsize, tysize);

	double** H_value_2D_Hx = Tptr2<double>(txsize, tysize);
	double** H_value_2D_Hy = Tptr2<double>(txsize, tysize);
	double** H_value_2D_Hz = Tptr2<double>(txsize, tysize);

	FILE* pf_ex, * pf_ey, * pf_ez;
	FILE* pf_hx, * pf_hy, * pf_hz;

	char Plane_save_name_Ex[256], Plane_save_name_Hx[256];
	char Plane_save_name_Ey[256], Plane_save_name_Hy[256];
	char Plane_save_name_Ez[256], Plane_save_name_Hz[256];

	strcpy(pFileName_tmp, Fname);
	sprintf(pFileName, "/Field_distribution_%04d_Ex.dat",  ss);
	strcat(pFileName_tmp, pFileName);
	strcpy(Plane_save_name_Ex, pFileName_tmp);

	strcpy(pFileName_tmp, Fname);
	sprintf(pFileName, "/Field_distribution_%04d_Ey.dat",  ss);
	strcat(pFileName_tmp, pFileName);
	strcpy(Plane_save_name_Ey, pFileName_tmp);

	strcpy(pFileName_tmp, Fname);
	sprintf(pFileName, "/Field_distribution_%04d_Ez.dat",  ss);
	strcat(pFileName_tmp, pFileName);
	strcpy(Plane_save_name_Ez, pFileName_tmp);

	strcpy(pFileName_tmp, Fname);
	sprintf(pFileName, "/Field_distribution_%04d_Hx.dat",  ss);
	strcat(pFileName_tmp, pFileName);
	strcpy(Plane_save_name_Hx, pFileName_tmp);

	strcpy(pFileName_tmp, Fname);
	sprintf(pFileName, "/Field_distribution_%04d_Hy.dat",  ss);
	strcat(pFileName_tmp, pFileName);
	strcpy(Plane_save_name_Hy, pFileName_tmp);

	strcpy(pFileName_tmp, Fname);
	sprintf(pFileName, "/Field_distribution_%04d_Hz.dat", ss);
	strcat(pFileName_tmp, pFileName);
	strcpy(Plane_save_name_Hz, pFileName_tmp);

	pf_ex = fopen(Plane_save_name_Ex, "wb"); pf_hx = fopen(Plane_save_name_Hx, "wb");
	pf_ey = fopen(Plane_save_name_Ey, "wb"); pf_hy = fopen(Plane_save_name_Hy, "wb");
	pf_ez = fopen(Plane_save_name_Ez, "wb"); pf_hz = fopen(Plane_save_name_Hz, "wb");
	
	findex = (kk - no_planes) / 2;
	for (int ii = 0; ii < txsize; ii++)
	{
		for (int jj = 0; jj < tysize; jj++)
		{
			if (sopd[kk].is_sopdf)
				if ((kk - no_planes) % 2 == 0) // 0, up
					pos = psEnvelop_Material[findex].pEnvelop[0][ii][jj];
				else // 1, down
					pos = psEnvelop_Material[findex].pEnvelop[1][ii][jj];
			else
				pos = sopd[kk].r_th;

			ExV = (EOx[ii][jj][pos] + EOx[ii][jj + 1][pos] +
				EOx[ii][jj][pos + 1] + EOx[ii][jj + 1][pos + 1]) / 4;
			EyV = (EOy[ii][jj][pos] + EOy[ii + 1][jj][pos] +
				EOy[ii][jj][pos + 1] + EOy[ii + 1][jj][pos + 1]) / 4;
			EzV = (EOz[ii][jj][pos] + EOz[ii + 1][jj][pos] +
				EOz[ii][jj + 1][pos] + EOz[ii + 1][jj + 1][pos]) / 4;

			HxV = (HOx[ii][jj][pos] + HOx[ii + 1][jj][pos]) / 2;
			HyV = (HOy[ii][jj][pos] + HOy[ii][jj + 1][pos]) / 2;
			HzV = (HOz[ii][jj][pos] + HOz[ii][jj][pos + 1]) / 2;

			E_value_2D_Ex[ii][jj] = ExV; H_value_2D_Hx[ii][jj] = HxV;
			E_value_2D_Ey[ii][jj] = EyV; H_value_2D_Hy[ii][jj] = HyV;
			E_value_2D_Ez[ii][jj] = EzV; H_value_2D_Hz[ii][jj] = HzV;
		}

		fwrite(E_value_2D_Ex[ii], sizeof(double), tysize, pf_ex);
		fwrite(E_value_2D_Ey[ii], sizeof(double), tysize, pf_ey);
		fwrite(E_value_2D_Ez[ii], sizeof(double), tysize, pf_ez);

		fwrite(H_value_2D_Hx[ii], sizeof(double), tysize, pf_hx);
		fwrite(H_value_2D_Hy[ii], sizeof(double), tysize, pf_hy);
		fwrite(H_value_2D_Hz[ii], sizeof(double), tysize, pf_hz);
	}

	memfree(&E_value_2D_Ex, txsize); memfree(&H_value_2D_Hx, txsize);
	memfree(&E_value_2D_Ey, txsize); memfree(&H_value_2D_Hy, txsize);
	memfree(&E_value_2D_Ez, txsize); memfree(&H_value_2D_Hz, txsize);

	fclose(pf_ex); fclose(pf_hx);
	fclose(pf_ey); fclose(pf_hy);
	fclose(pf_ez); fclose(pf_hz);
}
// YZ plane
void TSaveDistribution::Field_YZ_2D_Each(int ii, int ss, char* Fname)
{

	double** E_value_2D_Ex = Tptr2<double>(tysize, tzsize);
	double** E_value_2D_Ey = Tptr2<double>(tysize, tzsize);
	double** E_value_2D_Ez = Tptr2<double>(tysize, tzsize);

	double** H_value_2D_Hx = Tptr2<double>(tysize, tzsize);
	double** H_value_2D_Hy = Tptr2<double>(tysize, tzsize);
	double** H_value_2D_Hz = Tptr2<double>(tysize, tzsize);

	FILE* pf_ex, * pf_ey, * pf_ez;
	FILE* pf_hx, * pf_hy, * pf_hz;

	char Plane_save_name_Ex[256], Plane_save_name_Hx[256];
	char Plane_save_name_Ey[256], Plane_save_name_Hy[256];
	char Plane_save_name_Ez[256], Plane_save_name_Hz[256];

	strcpy(pFileName_tmp, Fname);
	sprintf(pFileName, "/Field_distribution_%04d_Ex.dat", ss);
	strcat(pFileName_tmp, pFileName);
	strcpy(Plane_save_name_Ex, pFileName_tmp);

	strcpy(pFileName_tmp, Fname);
	sprintf(pFileName, "/Field_distribution_%04d_Ey.dat", ss);
	strcat(pFileName_tmp, pFileName);
	strcpy(Plane_save_name_Ey, pFileName_tmp);

	strcpy(pFileName_tmp, Fname);
	sprintf(pFileName, "/Field_distribution_%04d_Ez.dat", ss);
	strcat(pFileName_tmp, pFileName);
	strcpy(Plane_save_name_Ez, pFileName_tmp);

	strcpy(pFileName_tmp, Fname);
	sprintf(pFileName, "/Field_distribution_%04d_Hx.dat", ss);
	strcat(pFileName_tmp, pFileName);
	strcpy(Plane_save_name_Hx, pFileName_tmp);

	strcpy(pFileName_tmp, Fname);
	sprintf(pFileName, "/Field_distribution_%04d_Hy.dat", ss);
	strcat(pFileName_tmp, pFileName);
	strcpy(Plane_save_name_Hy, pFileName_tmp);

	strcpy(pFileName_tmp, Fname);
	sprintf(pFileName, "/Field_distribution_%04d_Hz.dat", ss);
	strcat(pFileName_tmp, pFileName);
	strcpy(Plane_save_name_Hz, pFileName_tmp);

	pf_ex = fopen(Plane_save_name_Ex, "wb"); pf_hx = fopen(Plane_save_name_Hx, "wb");
	pf_ey = fopen(Plane_save_name_Ey, "wb"); pf_hy = fopen(Plane_save_name_Hy, "wb");
	pf_ez = fopen(Plane_save_name_Ez, "wb"); pf_hz = fopen(Plane_save_name_Hz, "wb");

	for (int jj = 0; jj < tysize; jj++)
	{
		for (int kk = 0; kk < tzsize; kk++)
		{
			ExV = (EOx[sopd[ii].r_th][jj][kk] +
				EOx[sopd[ii].r_th][jj + 1][kk] +
				EOx[sopd[ii].r_th][jj][kk + 1] +
				EOx[sopd[ii].r_th][jj + 1][kk + 1]) / 4;
			EyV = (EOy[sopd[ii].r_th][jj][kk] +
				EOy[sopd[ii].r_th + 1][jj][kk] +
				EOy[sopd[ii].r_th][jj][kk + 1] +
				EOy[sopd[ii].r_th + 1][jj][kk + 1]) / 4;
			EzV = (EOz[sopd[ii].r_th][jj][kk] +
				EOz[sopd[ii].r_th + 1][jj][kk] +
				EOz[sopd[ii].r_th][jj + 1][kk] +
				EOz[sopd[ii].r_th + 1][jj + 1][kk]) / 4;

			HxV = (HOx[sopd[ii].r_th][jj][kk] +
				HOx[sopd[ii].r_th + 1][jj][kk]) / 2;
			HyV = (HOy[sopd[ii].r_th][jj][kk] +
				HOy[sopd[ii].r_th][jj + 1][kk]) / 2;
			HzV = (HOz[sopd[ii].r_th][jj][kk] +
				HOz[sopd[ii].r_th][jj][kk + 1]) / 2;

			E_value_2D_Ex[jj][kk] = ExV;  H_value_2D_Hx[jj][kk] = HxV;
			E_value_2D_Ey[jj][kk] = EyV;	 H_value_2D_Hy[jj][kk] = HyV;
			E_value_2D_Ez[jj][kk] = EzV;	 H_value_2D_Hz[jj][kk] = HzV;
		}

		fwrite(E_value_2D_Ex[jj], sizeof(double), tzsize, pf_ex);
		fwrite(E_value_2D_Ey[jj], sizeof(double), tzsize, pf_ey);
		fwrite(E_value_2D_Ez[jj], sizeof(double), tzsize, pf_ez);
		fwrite(H_value_2D_Hx[jj], sizeof(double), tzsize, pf_hx);
		fwrite(H_value_2D_Hy[jj], sizeof(double), tzsize, pf_hy);
		fwrite(H_value_2D_Hz[jj], sizeof(double), tzsize, pf_hz);
	}

	memfree(&E_value_2D_Ex, tysize); memfree(&H_value_2D_Hx, tysize);
	memfree(&E_value_2D_Ey, tysize); memfree(&H_value_2D_Hy, tysize);
	memfree(&E_value_2D_Ez, tysize); memfree(&H_value_2D_Hz, tysize);

	fclose(pf_ex); fclose(pf_hx);
	fclose(pf_ey); fclose(pf_hy);
	fclose(pf_ez); fclose(pf_hz);
}

void TSaveDistribution::Field_XZ_2D_Each(int jj, int ss, char* Fname)
{
	double** E_value_2D_Ex = Tptr2<double>(txsize, tzsize);
	double** E_value_2D_Ey = Tptr2<double>(txsize, tzsize);
	double** E_value_2D_Ez = Tptr2<double>(txsize, tzsize);

	double** H_value_2D_Hx = Tptr2<double>(txsize, tzsize);
	double** H_value_2D_Hy = Tptr2<double>(txsize, tzsize);
	double** H_value_2D_Hz = Tptr2<double>(txsize, tzsize);

	FILE* pf_ex, * pf_ey, * pf_ez;
	FILE* pf_hx, * pf_hy, * pf_hz;

	char Plane_save_name_Ex[256], Plane_save_name_Hx[256];
	char Plane_save_name_Ey[256], Plane_save_name_Hy[256];
	char Plane_save_name_Ez[256], Plane_save_name_Hz[256];

	strcpy(pFileName_tmp, Fname);
	sprintf(pFileName, "/Field_distribution_%04d_Ex.dat",  ss);
	strcat(pFileName_tmp, pFileName);
	strcpy(Plane_save_name_Ex, pFileName_tmp);

	strcpy(pFileName_tmp, Fname);
	sprintf(pFileName, "/Field_distribution_%04d_Ey.dat",  ss);
	strcat(pFileName_tmp, pFileName);
	strcpy(Plane_save_name_Ey, pFileName_tmp);

	strcpy(pFileName_tmp, Fname);
	sprintf(pFileName, "/Field_distribution_%04d_Ez.dat",  ss);
	strcat(pFileName_tmp, pFileName);
	strcpy(Plane_save_name_Ez, pFileName_tmp);

	strcpy(pFileName_tmp, Fname);
	sprintf(pFileName, "/Field_distribution_%04d_Hx.dat",  ss);
	strcat(pFileName_tmp, pFileName);
	strcpy(Plane_save_name_Hx, pFileName_tmp);

	strcpy(pFileName_tmp, Fname);
	sprintf(pFileName, "/Field_distribution_%04d_Hy.dat",  ss);
	strcat(pFileName_tmp, pFileName);
	strcpy(Plane_save_name_Hy, pFileName_tmp);

	strcpy(pFileName_tmp, Fname);
	sprintf(pFileName, "/Field_distribution_%04d_Hz.dat",  ss);
	strcat(pFileName_tmp, pFileName);
	strcpy(Plane_save_name_Hz, pFileName_tmp);

	pf_ex = fopen(Plane_save_name_Ex, "wb"); pf_hx = fopen(Plane_save_name_Hx, "wb");
	pf_ey = fopen(Plane_save_name_Ey, "wb"); pf_hy = fopen(Plane_save_name_Hy, "wb");
	pf_ez = fopen(Plane_save_name_Ez, "wb"); pf_hz = fopen(Plane_save_name_Hz, "wb");

	for (int ii = 0; ii < txsize; ii++)
	{
		for (int kk = 0; kk < tzsize; kk++)
		{
			ExV = (EOx[ii][sopd[jj].r_th][kk] +
				EOx[ii][sopd[jj].r_th + 1][kk] +
				EOx[ii][sopd[jj].r_th][kk + 1] +
				EOx[ii][sopd[jj].r_th + 1][kk + 1]) / 4;
			EyV = (EOy[ii][sopd[jj].r_th][kk] +
				EOy[ii + 1][sopd[jj].r_th][kk] +
				EOy[ii][sopd[jj].r_th][kk + 1] +
				EOy[ii + 1][sopd[jj].r_th][kk + 1]) / 4;
			EzV = (EOz[ii][sopd[jj].r_th][kk] +
				EOz[ii + 1][sopd[jj].r_th][kk] +
				EOz[ii][sopd[jj].r_th + 1][kk] +
				EOz[ii + 1][sopd[jj].r_th + 1][kk]) / 4;

			HxV = (HOx[ii][sopd[jj].r_th][kk] +
				HOx[ii + 1][sopd[jj].r_th][kk]) / 2;
			HyV = (HOy[ii][sopd[jj].r_th][kk] +
				HOy[ii][sopd[jj].r_th + 1][kk]) / 2;
			HzV = (HOz[ii][sopd[jj].r_th][kk] +
				HOz[ii][sopd[jj].r_th][kk + 1]) / 2;

			E_value_2D_Ex[ii][kk] = ExV; H_value_2D_Hx[ii][kk] = HxV;
			E_value_2D_Ey[ii][kk] = EyV;	H_value_2D_Hy[ii][kk] = HyV;
			E_value_2D_Ez[ii][kk] = EzV;	H_value_2D_Hz[ii][kk] = HzV;
		}

		fwrite(E_value_2D_Ex[ii], sizeof(double), tzsize, pf_ex);
		fwrite(E_value_2D_Ey[ii], sizeof(double), tzsize, pf_ey);
		fwrite(E_value_2D_Ez[ii], sizeof(double), tzsize, pf_ez);

		fwrite(H_value_2D_Hx[ii], sizeof(double), tzsize, pf_hx);
		fwrite(H_value_2D_Hy[ii], sizeof(double), tzsize, pf_hy);
		fwrite(H_value_2D_Hz[ii], sizeof(double), tzsize, pf_hz);
	}
	memfree(&E_value_2D_Ex, txsize); memfree(&H_value_2D_Hx, txsize);
	memfree(&E_value_2D_Ey, txsize); memfree(&H_value_2D_Hy, txsize);
	memfree(&E_value_2D_Ez, txsize); memfree(&H_value_2D_Hz, txsize);

	fclose(pf_ex); fclose(pf_hx);
	fclose(pf_ey); fclose(pf_hy);
	fclose(pf_ez); fclose(pf_hz);
}

// Frequency domain
void TSaveDistribution::save_FieldDistribution_assign_file_pointer(int nn)
{
	strcpy(Fname_save_field_Pathcpy_amp, Plane_save_name[nn]);
	strcpy(Fname_save_field_Pathcpy_phase, Plane_save_name[nn]);
	sprintf(Fname_save_field_Name_temp, "/Ex_Amp[%1.2le].dat", pFc->f0);
	strcpy(Fname_save_field_Amp, strcat(Fname_save_field_Pathcpy_amp, Fname_save_field_Name_temp));
	sprintf(Fname_save_field_Name_temp, "/Ex_Phase[%1.2le].dat", pFc->f0);
	strcpy(Fname_save_field_Phase, strcat(Fname_save_field_Pathcpy_phase, Fname_save_field_Name_temp));
	pFEx_Amp = fopen(Fname_save_field_Amp, "wb");
	pFEx_Phase = fopen(Fname_save_field_Phase, "wb");

	strcpy(Fname_save_field_Pathcpy_amp, Plane_save_name[nn]);
	strcpy(Fname_save_field_Pathcpy_phase, Plane_save_name[nn]);
	sprintf(Fname_save_field_Name_temp, "/Ey_Amp[%1.2le].dat", pFc->f0);
	strcpy(Fname_save_field_Amp, strcat(Fname_save_field_Pathcpy_amp, Fname_save_field_Name_temp));
	sprintf(Fname_save_field_Name_temp, "/Ey_Phase[%1.2le].dat", pFc->f0);
	strcpy(Fname_save_field_Phase, strcat(Fname_save_field_Pathcpy_phase, Fname_save_field_Name_temp));
	pFEy_Amp = fopen(Fname_save_field_Amp, "wb");
	pFEy_Phase = fopen(Fname_save_field_Phase, "wb");

	strcpy(Fname_save_field_Pathcpy_amp, Plane_save_name[nn]);
	strcpy(Fname_save_field_Pathcpy_phase, Plane_save_name[nn]);
	sprintf(Fname_save_field_Name_temp, "/Ez_Amp[%1.2le].dat", pFc->f0);
	strcpy(Fname_save_field_Amp, strcat(Fname_save_field_Pathcpy_amp, Fname_save_field_Name_temp));
	sprintf(Fname_save_field_Name_temp, "/Ez_Phase[%1.2le].dat", pFc->f0);
	strcpy(Fname_save_field_Phase, strcat(Fname_save_field_Pathcpy_phase, Fname_save_field_Name_temp));
	pFEz_Amp = fopen(Fname_save_field_Amp, "wb");
	pFEz_Phase = fopen(Fname_save_field_Phase, "wb");
}

void TSaveDistribution::Save_POI_frequency(int Nststep_iter)
{
	int steady_state_offset = 0;
	complex <double> temp;
	double f1, f2;
	if (Nststep_iter == save_sta)
	{
		tt1 = Nststep_iter;
		for (int pp = 0; pp < total_no_planes; pp++)
		{
			switch (sopd[pp].s_plane)
			{
			case 1:
				SaveFieldXY(pp, 1);
				break;
			case 2:
				SaveFieldYZ(pp, 1);
				break;
			case 3:
				SaveFieldXZ(pp, 1);
				break;
			}
		}
		// save 3D

	}
	if (Nststep_iter == save_end)
	{
		// optical chirality 계산 가능.
		tt2 = Nststep_iter;
		for (int pp = 0; pp < total_no_planes; pp++)
		{
			switch (sopd[pp].s_plane)
			{
			case 1:
				SaveFieldXY(pp, 2);
				break;
			case 2:
				SaveFieldYZ(pp, 2);
				break;
			case 3:
				SaveFieldXZ(pp, 2);
				break;
			}
		}

		Calcuate_Amp_Phase();
		save_FieldDistribution();

		// cal_amp_phae 3D
		// save 3D
	}
}
void TSaveDistribution::save_FieldDistribution()
{
	for (int pp = 0; pp < total_no_planes; pp++)
	{
		//for (int ff = 0; ff < sampling; ff++)
		//{}
		save_FieldDistribution_assign_file_pointer(pp);

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
					pSave_ExAmp[jj] = sopds[pp].AmpEx[ii][jj];
					pSave_ExPhase[jj] = sopds[pp].PhaseEx[ii][jj];
					pSave_EyAmp[jj] = sopds[pp].AmpEy[ii][jj];
					pSave_EyPhase[jj] = sopds[pp].PhaseEy[ii][jj];
					pSave_EzAmp[jj] = sopds[pp].AmpEz[ii][jj];
					pSave_EzPhase[jj] = sopds[pp].PhaseEz[ii][jj];
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
					pSave_ExAmp[kk] = sopds[pp].AmpEx[jj][kk];
					pSave_ExPhase[kk] = sopds[pp].PhaseEx[jj][kk];
					pSave_EyAmp[kk] = sopds[pp].AmpEy[jj][kk];
					pSave_EyPhase[kk] = sopds[pp].PhaseEy[jj][kk];
					pSave_EzAmp[kk] = sopds[pp].AmpEz[jj][kk];
					pSave_EzPhase[kk] = sopds[pp].PhaseEz[jj][kk];
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
					pSave_ExAmp[kk] = sopds[pp].AmpEx[ii][kk];
					pSave_ExPhase[kk] = sopds[pp].PhaseEx[ii][kk];
					pSave_EyAmp[kk] = sopds[pp].AmpEy[ii][kk];
					pSave_EyPhase[kk] = sopds[pp].PhaseEy[ii][kk];
					pSave_EzAmp[kk] = sopds[pp].AmpEz[ii][kk];
					pSave_EzPhase[kk] = sopds[pp].PhaseEz[ii][kk];
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

void TSaveDistribution::SaveFieldXY(int pp, int t1t2)
{
	if (t1t2 == 1)
	{
		findex = (pp - no_planes) / 2;
		for (int ii = 0; ii < txsize; ii++)
			for (int jj = 0; jj < tysize; jj++)
			{
				if (sopd[pp].is_sopdf)
					if ((pp - no_planes) % 2 == 0) // 0, up
						pos = psEnvelop_Material[findex].pEnvelop[0][ii][jj];
					else // 1, down
						pos = psEnvelop_Material[findex].pEnvelop[1][ii][jj];
				else
					pos = sopd[pp].r_th;

				sopds[pp].Ex_save_plane1[ii][jj] = (EOx[ii][jj][pos] + EOx[ii][jj + 1][pos] +
					EOx[ii][jj][pos + 1] + EOx[ii][jj + 1][pos + 1]) / 4.0;
				sopds[pp].Ey_save_plane1[ii][jj] = (EOy[ii][jj][pos] + EOy[ii + 1][jj][pos] +
					EOy[ii][jj][pos + 1] + EOy[ii + 1][jj][pos + 1]) / 4.0;
				sopds[pp].Ez_save_plane1[ii][jj] = (EOz[ii][jj][pos] + EOz[ii + 1][jj][pos] +
					EOz[ii][jj + 1][pos + 1] + EOz[ii + 1][jj + 1][pos]) / 4.0;
			}
	}
	if (t1t2 == 2)
	{
		findex = (pp - no_planes) / 2;
		for (int ii = 0; ii < txsize; ii++)
			for (int jj = 0; jj < tysize; jj++)
			{
				if (sopd[pp].is_sopdf)
					if ((pp - no_planes) % 2 == 0) // 0
						pos = psEnvelop_Material[findex].pEnvelop[0][ii][jj];
					else // 1
						pos = psEnvelop_Material[findex].pEnvelop[1][ii][jj];
				else
					pos = sopd[pp].r_th;

				sopds[pp].Ex_save_plane2[ii][jj] = (EOx[ii][jj][pos] + EOx[ii][jj + 1][pos] +
					EOx[ii][jj][pos + 1] + EOx[ii][jj + 1][pos + 1]) / 4.0;
				sopds[pp].Ey_save_plane2[ii][jj] = (EOy[ii][jj][pos] + EOy[ii + 1][jj][pos] +
					EOy[ii][jj][pos + 1] + EOy[ii + 1][jj][pos + 1]) / 4.0;
				sopds[pp].Ez_save_plane2[ii][jj] = (EOz[ii][jj][pos] + EOz[ii + 1][jj][pos] +
					EOz[ii][jj + 1][pos + 1] + EOz[ii + 1][jj + 1][pos]) / 4.0;
			}
	}

}
void TSaveDistribution::SaveFieldYZ(int pp, int t1t2)
{
	if (t1t2 == 1)
	{
		findex = (pp - no_planes) / 2;
		for (int jj = 0; jj < tysize; jj++)
			for (int kk = 0; kk < tzsize; kk++)
			{
				if (sopd[pp].is_sopdf)
					if ((pp - no_planes) % 2 == 0) // 0
						pos = psEnvelop_Material[findex].pEnvelop[0][jj][kk];
					else // 1
						pos = psEnvelop_Material[findex].pEnvelop[1][jj][kk];
				else
					pos = sopd[pp].r_th;

				sopds[pp].Ex_save_plane1[jj][kk] = (EOx[pos][jj][kk] + EOx[pos][jj + 1][kk] +
					EOx[pos][jj][kk + 1] + EOx[pos][jj + 1][kk + 1]) / 4.0;
				sopds[pp].Ey_save_plane1[jj][kk] = (EOy[pos][jj][kk] + EOy[pos + 1][jj][kk] +
					EOy[pos][jj][kk + 1] + EOy[pos + 1][jj][kk + 1]) / 4.0;
				sopds[pp].Ez_save_plane1[jj][kk] = (EOz[pos][jj][kk] + EOz[pos + 1][jj][kk] +
					EOz[pos][jj][kk + 1] + EOz[pos + 1][jj + 1][kk]) / 4.0;
			}
	}
	if (t1t2 == 2)
	{
		findex = (pp - no_planes) / 2;
		for (int jj = 0; jj < tysize; jj++)
			for (int kk = 0; kk < tzsize; kk++)
			{
				if (sopd[pp].is_sopdf)
					if ((pp - no_planes) % 2 == 0) // 0
						pos = psEnvelop_Material[pp - no_planes - 1].pEnvelop[0][jj][kk];
					else // 1
						pos = psEnvelop_Material[pp - no_planes - 1].pEnvelop[1][jj][kk];
				else
					pos = sopd[pp].r_th;

				sopds[pp].Ex_save_plane2[jj][kk] = (EOx[pos][jj][kk] + EOx[pos][jj + 1][kk] +
					EOx[pos][jj][kk + 1] + EOx[pos][jj + 1][kk + 1]) / 4.0;
				sopds[pp].Ey_save_plane2[jj][kk] = (EOy[pos][jj][kk] + EOy[pos + 1][jj][kk] +
					EOy[pos][jj][kk + 1] + EOy[pos + 1][jj][kk + 1]) / 4.0;
				sopds[pp].Ez_save_plane2[jj][kk] = (EOz[pos][jj][kk] + EOz[pos + 1][jj][kk] +
					EOz[pos][jj][kk + 1] + EOz[pos + 1][jj + 1][kk]) / 4.0;
			}
	}
}
void TSaveDistribution::SaveFieldXZ(int pp, int t1t2)
{
	if (t1t2 == 1)
	{
		findex = (pp - no_planes) / 2;
		for (int ii = 0; ii < txsize; ii++)
			for (int kk = 0; kk < tzsize; kk++)
			{
				if (sopd[pp].is_sopdf)
					if ((pp - no_planes) % 2 == 0) // 0
						pos = psEnvelop_Material[pp - no_planes - 1].pEnvelop[0][ii][kk];
					else // 1
						pos = psEnvelop_Material[pp - no_planes - 1].pEnvelop[1][ii][kk];
				else
					pos = sopd[pp].r_th;

				sopds[pp].Ex_save_plane1[ii][kk] = (EOx[ii][pos][kk] + EOx[ii][pos + 1][kk] +
					EOx[ii][pos][kk + 1] + EOx[ii][pos + 1][kk + 1]) / 4.0;
				sopds[pp].Ey_save_plane1[ii][kk] = (EOy[ii][pos][kk] + EOy[ii + 1][pos][kk] +
					EOy[ii][pos][kk + 1] + EOy[ii + 1][pos][kk + 1]) / 4.0;
				sopds[pp].Ez_save_plane1[ii][kk] = (EOz[ii][pos][kk] + EOz[ii + 1][pos][kk] +
					EOz[ii][pos + 1][kk + 1] + EOz[ii + 1][pos + 1][kk]) / 4.0;
			}
	}
	if (t1t2 == 2)
	{
		findex = (pp - no_planes) / 2;
		for (int ii = 0; ii < txsize; ii++)
			for (int kk = 0; kk < tzsize; kk++)
			{
				if (sopd[pp].is_sopdf)
					if ((pp - no_planes) % 2 == 0) // 0
						pos = psEnvelop_Material[pp - no_planes - 1].pEnvelop[0][ii][kk];
					else // 1
						pos = psEnvelop_Material[pp - no_planes - 1].pEnvelop[1][ii][kk];
				else
					pos = sopd[pp].r_th;

				sopds[pp].Ex_save_plane2[ii][kk] = (EOx[ii][pos][kk] + EOx[ii][pos + 1][kk] +
					EOx[ii][pos][kk + 1] + EOx[ii][pos + 1][kk + 1]) / 4.0;
				sopds[pp].Ey_save_plane2[ii][kk] = (EOy[ii][pos][kk] + EOy[ii + 1][pos][kk] +
					EOy[ii][pos][kk + 1] + EOy[ii + 1][pos][kk + 1]) / 4.0;
				sopds[pp].Ez_save_plane2[ii][kk] = (EOz[ii][pos][kk] + EOz[ii + 1][pos][kk] +
					EOz[ii][pos + 1][kk + 1] + EOz[ii + 1][pos + 1][kk]) / 4.0;
			}
	}
}


void TSaveDistribution::Calcuate_Amp_PhaseXY(int pp)
{
	double f1, f2;
	string Save_Distribution_debug = "Test_Amp_Phase.txt";
	ofstream pf(Save_Distribution_debug.data());

	for (int ii = 0; ii < txsize; ii++)
		for (int jj = 0; jj < tysize; jj++)
		{
			f1 = sopds[pp].Ex_save_plane1[ii][jj];
			f2 = sopds[pp].Ex_save_plane2[ii][jj];

			//Phase = atan2(cos(2 * pi*f0*dt*tt2) - f2 / f1, sin(2 * pi*f0*tt2));
			Phase = atan2((f2 * sin(2 * pi * f0 * dt * tt1) - f1 * sin(2 * pi * f0 * dt * tt2)),
				(f1 * cos(2 * pi * f0 * dt * tt2) - f2 * cos(2 * pi * f0 * dt * tt1)));
			Amp = f1 / sin(2 * pi * f0 * dt * tt1 + Phase);

			//pf << "f1 = " << f1 << " " << "f2 =" << f2 << " " << "tt1 = " << tt1 << " " << "tt2 = " << tt2 << " "
			//	<< "Amp = " << Amp << " " << "Phase = " << Phase << endl;
			Amp_Phase_Tunning();
			sopds[pp].AmpEx[ii][jj] = Amp;
			sopds[pp].PhaseEx[ii][jj] = Phase;


			f1 = sopds[pp].Ey_save_plane1[ii][jj];
			f2 = sopds[pp].Ey_save_plane2[ii][jj];
			//Phase = atan2(cos(2 * pi*f0*dt*tt2) - f2 / f1, sin(2 * pi*f0*tt2));
			Phase = atan2((f2 * sin(2 * pi * f0 * dt * tt1) - f1 * sin(2 * pi * f0 * dt * tt2)),
				(f1 * cos(2 * pi * f0 * dt * tt2) - f2 * cos(2 * pi * f0 * dt * tt1)));
			Amp = f1 / sin(2 * pi * f0 * dt * tt1 + Phase);
			Amp_Phase_Tunning();
			sopds[pp].AmpEy[ii][jj] = Amp;
			sopds[pp].PhaseEy[ii][jj] = Phase;


			f1 = sopds[pp].Ez_save_plane1[ii][jj];
			f2 = sopds[pp].Ez_save_plane2[ii][jj];
			//Phase = atan2(cos(2 * pi*f0*dt*tt2) - f2 / f1, sin(2 * pi*f0*tt2));
			Phase = atan2((f2 * sin(2 * pi * f0 * dt * tt1) - f1 * sin(2 * pi * f0 * dt * tt2)),
				(f1 * cos(2 * pi * f0 * dt * tt2) - f2 * cos(2 * pi * f0 * dt * tt1)));
			Amp = f1 / sin(2 * pi * f0 * dt * tt1 + Phase);
			Amp_Phase_Tunning();
			sopds[pp].AmpEz[ii][jj] = Amp;
			sopds[pp].PhaseEz[ii][jj] = Phase;
		}

	pf.close();
}
void TSaveDistribution::Calcuate_Amp_PhaseYZ(int pp)
{
	double f1, f2;
	for (int jj = 0; jj < tysize; jj++)
		for (int kk = 0; kk < tzsize; kk++)
		{
			f1 = sopds[pp].Ex_save_plane1[jj][kk];
			f2 = sopds[pp].Ex_save_plane2[jj][kk];
			//Phase = atan2(cos(2 * pi*f0*dt*tt2) - f2 / f1, sin(2 * pi*f0*tt2));
			Phase = atan2((f2 * sin(2 * pi * f0 * dt * tt1) - f1 * sin(2 * pi * f0 * dt * tt2)),
				(f1 * cos(2 * pi * f0 * dt * tt2) - f2 * cos(2 * pi * f0 * dt * tt1)));
			Amp = f1 / sin(2 * pi * f0 * dt * tt1 + Phase);
			Amp_Phase_Tunning();
			sopds[pp].AmpEx[jj][kk] = Amp;
			sopds[pp].PhaseEx[jj][kk] = Phase;


			f1 = sopds[pp].Ey_save_plane1[jj][kk];
			f2 = sopds[pp].Ey_save_plane2[jj][kk];
			//Phase = atan2(cos(2 * pi*f0*dt*tt2) - f2 / f1, sin(2 * pi*f0*tt2));
			Phase = atan2((f2 * sin(2 * pi * f0 * dt * tt1) - f1 * sin(2 * pi * f0 * dt * tt2)),
				(f1 * cos(2 * pi * f0 * dt * tt2) - f2 * cos(2 * pi * f0 * dt * tt1)));
			Amp = f1 / sin(2 * pi * f0 * dt * tt1 + Phase);
			Amp_Phase_Tunning();
			sopds[pp].AmpEy[jj][kk] = Amp;
			sopds[pp].PhaseEy[jj][kk] = Phase;


			f1 = sopds[pp].Ez_save_plane1[jj][kk];
			f2 = sopds[pp].Ez_save_plane2[jj][kk];
			//Phase = atan2(cos(2 * pi*f0*dt*tt2) - f2 / f1, sin(2 * pi*f0*tt2));
			Phase = atan2((f2 * sin(2 * pi * f0 * dt * tt1) - f1 * sin(2 * pi * f0 * dt * tt2)),
				(f1 * cos(2 * pi * f0 * dt * tt2) - f2 * cos(2 * pi * f0 * dt * tt1)));
			Amp = f1 / sin(2 * pi * f0 * dt * tt1 + Phase);
			Amp_Phase_Tunning();
			sopds[pp].AmpEz[jj][kk] = Amp;
			sopds[pp].PhaseEz[jj][kk] = Phase;
		}
}

void TSaveDistribution::Calcuate_Amp_PhaseXZ(int pp)
{
	double f1, f2;
	for (int ii = 0; ii < txsize; ii++)
		for (int kk = 0; kk < tzsize; kk++)
		{
			f1 = sopds[pp].Ex_save_plane1[ii][kk];
			f2 = sopds[pp].Ex_save_plane2[ii][kk];
			//Phase = atan2(cos(2 * pi*f0*dt*tt2) - f2 / f1, sin(2 * pi*f0*tt2));
			Phase = atan2((f2 * sin(2 * pi * f0 * dt * tt1) - f1 * sin(2 * pi * f0 * dt * tt2)),
				(f1 * cos(2 * pi * f0 * dt * tt2) - f2 * cos(2 * pi * f0 * dt * tt1)));
			Amp = f1 / sin(2 * pi * f0 * dt * tt1 + Phase);
			Amp_Phase_Tunning();
			sopds[pp].AmpEx[ii][kk] = Amp;
			sopds[pp].PhaseEx[ii][kk] = Phase;


			f1 = sopds[pp].Ey_save_plane1[ii][kk];
			f2 = sopds[pp].Ey_save_plane2[ii][kk];
			//Phase = atan2(cos(2 * pi*f0*dt*tt2) - f2 / f1, sin(2 * pi*f0*tt2));
			Phase = atan2((f2 * sin(2 * pi * f0 * dt * tt1) - f1 * sin(2 * pi * f0 * dt * tt2)),
				(f1 * cos(2 * pi * f0 * dt * tt2) - f2 * cos(2 * pi * f0 * dt * tt1)));
			Amp = f1 / sin(2 * pi * f0 * dt * tt1 + Phase);
			Amp_Phase_Tunning();
			sopds[pp].AmpEy[ii][kk] = Amp;
			sopds[pp].PhaseEy[ii][kk] = Phase;


			f1 = sopds[pp].Ez_save_plane1[ii][kk];
			f2 = sopds[pp].Ez_save_plane2[ii][kk];
			//Phase = atan2(cos(2 * pi*f0*dt*tt2) - f2 / f1, sin(2 * pi*f0*tt2));
			Phase = atan2((f2 * sin(2 * pi * f0 * dt * tt1) - f1 * sin(2 * pi * f0 * dt * tt2)),
				(f1 * cos(2 * pi * f0 * dt * tt2) - f2 * cos(2 * pi * f0 * dt * tt1)));
			Amp = f1 / sin(2 * pi * f0 * dt * tt1 + Phase);
			Amp_Phase_Tunning();
			sopds[pp].AmpEz[ii][kk] = Amp;
			sopds[pp].PhaseEz[ii][kk] = Phase;
		}
}

void TSaveDistribution::Calcuate_Amp_Phase()
{
	for (int pp = 0; pp < total_no_planes; pp++)
	{
		switch (sopd[pp].s_plane)
		{
		case 1:
			Calcuate_Amp_PhaseXY(pp);
			break;
		case 2:
			Calcuate_Amp_PhaseYZ(pp);
			break;
		case 3:
			Calcuate_Amp_PhaseXZ(pp);
			break;
		}
	}
}


void TSaveDistribution::MemoryAssign_for_FieldDistribution()
{
	// 주파수 별로 sopds 개수 할당,
	// 각 주파수의 평면값 할당.
	// 공통 : 주파수, 평면.
	sopds = Tptr1<struct_Obs_Plane_frequency>(total_no_planes);
	for (int pp = 0; pp < total_no_planes; pp++)
	{
		switch (sopd[pp].s_plane)
		{
		case 1:
			sopds[pp].Ex_save_plane1 = Tptr2<double>(txsize, tysize);
			sopds[pp].Ey_save_plane1 = Tptr2<double>(txsize, tysize);
			sopds[pp].Ez_save_plane1 = Tptr2<double>(txsize, tysize);
			sopds[pp].Ex_save_plane2 = Tptr2<double>(txsize, tysize);
			sopds[pp].Ey_save_plane2 = Tptr2<double>(txsize, tysize);
			sopds[pp].Ez_save_plane2 = Tptr2<double>(txsize, tysize);

			sopds[pp].AmpEx = Tptr2<double>(txsize, tysize);
			sopds[pp].AmpEy = Tptr2<double>(txsize, tysize);
			sopds[pp].AmpEz = Tptr2<double>(txsize, tysize);
			sopds[pp].PhaseEx = Tptr2<double>(txsize, tysize);
			sopds[pp].PhaseEy = Tptr2<double>(txsize, tysize);
			sopds[pp].PhaseEz = Tptr2<double>(txsize, tysize);
			break;
		case 2:
			sopds[pp].Ex_save_plane1 = Tptr2<double>(tysize, tzsize);
			sopds[pp].Ey_save_plane1 = Tptr2<double>(tysize, tzsize);
			sopds[pp].Ez_save_plane1 = Tptr2<double>(tysize, tzsize);
			sopds[pp].Ex_save_plane2 = Tptr2<double>(tysize, tzsize);
			sopds[pp].Ey_save_plane2 = Tptr2<double>(tysize, tzsize);
			sopds[pp].Ez_save_plane2 = Tptr2<double>(tysize, tzsize);
			sopds[pp].AmpEx = Tptr2<double>(tysize, tzsize);
			sopds[pp].AmpEy = Tptr2<double>(tysize, tzsize);
			sopds[pp].AmpEz = Tptr2<double>(tysize, tzsize);
			sopds[pp].PhaseEx = Tptr2<double>(tysize, tzsize);
			sopds[pp].PhaseEy = Tptr2<double>(tysize, tzsize);
			sopds[pp].PhaseEz = Tptr2<double>(tysize, tzsize);
			break;
		case 3:
			sopds[pp].Ex_save_plane1 = Tptr2<double>(txsize, tzsize);
			sopds[pp].Ey_save_plane1 = Tptr2<double>(txsize, tzsize);
			sopds[pp].Ez_save_plane1 = Tptr2<double>(txsize, tzsize);
			sopds[pp].Ex_save_plane2 = Tptr2<double>(txsize, tzsize);
			sopds[pp].Ey_save_plane2 = Tptr2<double>(txsize, tzsize);
			sopds[pp].Ez_save_plane2 = Tptr2<double>(txsize, tzsize);
			sopds[pp].AmpEx = Tptr2<double>(txsize, tzsize);
			sopds[pp].AmpEy = Tptr2<double>(txsize, tzsize);
			sopds[pp].AmpEz = Tptr2<double>(txsize, tzsize);
			sopds[pp].PhaseEx = Tptr2<double>(txsize, tzsize);
			sopds[pp].PhaseEy = Tptr2<double>(txsize, tzsize);
			sopds[pp].PhaseEz = Tptr2<double>(txsize, tzsize);
			break;
		}
	}
}
