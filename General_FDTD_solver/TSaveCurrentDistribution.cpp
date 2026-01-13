#include "TSaveCurrentDistribution.h"
TSaveCurrentDistribution::TSaveCurrentDistribution()
{
	_mkdir("result\\Current");
}

/*
void TSaveCurrentDistribution::save_Current_density()
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
				strcpy(Fname_Current_density_Amp, strcat(Fname_Current_density_Path, Fname_Current_density_Name_temp));
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
void TSaveCurrentDistribution::save_nXh()
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
void TSaveCurrentDistribution::Read_save_info()
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
		int noFlatnum = Jdocument["SIMULATION_ENVIORONMENT"]["Distribution_notFlat"].Size();
		int*** pEnvelop_Material = Tptr3<int>(2, txsize, tysize);

		ActualNum = 0;
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
	for (int nn = 0; nn < ActualNum; nn++)
	{
		sprintf(Plane_save_name[no_planes + 2 * nn], "./result/distribution/XY_Plane_Envelop[%s]_up", psEnvelop_Material[nn].Name);
		sprintf(Plane_save_name[no_planes + 2 * nn + 1], "./result/distribution/XY_Plane_Envelop[%s]_down", psEnvelop_Material[nn].Name);
		_mkdir(Plane_save_name[no_planes + 2 * nn]);
		_mkdir(Plane_save_name[no_planes + 2 * nn + 1]);
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
*/