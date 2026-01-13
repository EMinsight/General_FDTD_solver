//---------------------------------------------------------------------------
#pragma hdrstop
#include "TSaveSnapshot.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)

__fastcall TSaveSnapshot::TSaveSnapshot(double ***Ex, double ***Ey, double ***Ez, double ***Hx, double ***Hy, double ***Hz)
{
	this->EOx = Ex;
	this->EOy = Ey;
	this->EOz = Ez;
	this->HOx = Hx;
	this->HOy = Hy;
	this->HOz = Hz;
	strcpy(JsonFileName, "./IEM_solver_input_format.json");
	isChild = 0;
}
__fastcall TSaveSnapshot::~TSaveSnapshot()
{
	memfree(&sop); 
	memfree(&sopd);
}

void TSaveSnapshot::get_Obs_point(int *xx1, int *yy1, int *zz1, int *xx2, int *yy2, int *zz2)
{
    if(no_points<1)
    {
        return;
    }
    else if(no_points<2)
    {
        *xx1 = sop[0].Ox;
        *yy1 = sop[0].Oy;
        *zz1 = sop[0].Oz;
    }
    else if(no_points>=2)
    {
        *xx1 = sop[0].Ox;
        *yy1 = sop[0].Oy;
        *zz1 = sop[0].Oz;
        *xx2 = sop[1].Ox;
        *yy2 = sop[1].Oy;
        *zz2 = sop[1].Oz;
    }
}

void TSaveSnapshot::set_TimeDiff_NstepMax(double dt, int Nststep_max)
{
	this->dt = dt;
	if (this->Nststep_Max != Nststep_max)
		this->Nststep_Max = Nststep_max;
	Nststep_save = Nststep_max*0.01;
}

void TSaveSnapshot::set_spatial_and_dt(double dx, double dy, double dz, double dt)
{
	this->dx = dx;
	this->dy = dy;
	this->dz = dz;
	this->dt = dt;
}

void TSaveSnapshot::Set_max_iter()
{
	time_obs = Tptr1<double>(Nststep_save);
	time_comp = Tptr1<double>(Nststep_save);
	for (int ii = 0; ii < no_points; ii++)
	{
		sop[ii].Ex_save_array = Tptr1<double>(Nststep_save);
		sop[ii].Ey_save_array = Tptr1<double>(Nststep_save);
		sop[ii].Ez_save_array = Tptr1<double>(Nststep_save);

		sop[ii].Hx_save_array = Tptr1<double>(Nststep_save);
		sop[ii].Hy_save_array = Tptr1<double>(Nststep_save);
		sop[ii].Hz_save_array = Tptr1<double>(Nststep_save);
	}
	for (int ii = 0; ii < no_composit; ii++)
	{
		pSOC[ii].SumEx = Tptr1<double>(Nststep_save);
		pSOC[ii].SumEy = Tptr1<double>(Nststep_save);
	}
}

void TSaveSnapshot::get_mesh_length(int txsize, int tysize, int tzsize)
{
	this->txsize = txsize;
	this->tysize = tysize;
	this->tzsize = tzsize;
}

void TSaveSnapshot::Read_save_info()
{
	int dummy;
	Document Jdocument;
	ifstream ifs(JsonFileName);
	IStreamWrapper isw(ifs);
	Jdocument.ParseStream(isw);

	no_points = Jdocument["SIMULATION_ENVIORONMENT"]["Observation"].Size();
	// 관찰점의 갯수 및 관찰점 지정
	sop = Tptr1<struct_Obs_Point>(no_points);
	Obs_save_name = Tptr2<char>(no_points, 256);
	double* temp_coordx, * temp_coordy, * temp_coordz;
	double maxx, maxy, maxz;

	txsize = Jdocument["MESH"]["Voxelsize"]["txsize"].GetDouble();
	tysize = Jdocument["MESH"]["Voxelsize"]["tysize"].GetDouble();
	tzsize = Jdocument["MESH"]["Voxelsize"]["tzsize"].GetDouble();

	maxx = txsize+1;
	maxy = tysize+1;
	maxz = tzsize+1;
	temp_coordx = Tptr1<double>(maxx);
	temp_coordy = Tptr1<double>(maxy);
	temp_coordz = Tptr1<double>(maxz);
	for (int ii = 0; ii < maxx; ii++)
		temp_coordx[ii] = Jdocument["MESH"]["CellLocation"]["xcoordinate"][ii].GetDouble();

	for (int ii = 0; ii < maxy; ii++)
		temp_coordy[ii] = Jdocument["MESH"]["CellLocation"]["ycoordinate"][ii].GetDouble();

	for (int ii = 0; ii < maxz; ii++)
		temp_coordz[ii] = Jdocument["MESH"]["CellLocation"]["zcoordinate"][ii].GetDouble();


	for (int nn = 0; nn < no_points; nn++)
	{
		sop[nn].Rox = Jdocument["SIMULATION_ENVIORONMENT"]["Observation"][nn]["Location"][0].GetDouble();
		sop[nn].Roy = Jdocument["SIMULATION_ENVIORONMENT"]["Observation"][nn]["Location"][1].GetDouble();
		sop[nn].Roz = Jdocument["SIMULATION_ENVIORONMENT"]["Observation"][nn]["Location"][2].GetDouble();
		
		for (int ii = 0; ii < maxx; ii++)
			if (temp_coordx[ii] == sop[nn].Rox)
				sop[nn].Ox = ii;
			else
				continue;
		for (int ii = 0; ii < maxy; ii++)
			if (temp_coordy[ii] == sop[nn].Roy)
				sop[nn].Oy = ii;
			else
				continue;

		for (int ii = 0; ii < maxz; ii++)
			if (temp_coordz[ii] == sop[nn].Roz)
				sop[nn].Oz = ii;
			else
				continue;

		sprintf(Obs_save_name[nn], "./result/observation/Observation_[x=%2.2e,y=%2.2e,z=%2.2e]_field_plot.dat", sop[nn].Rox, sop[nn].Roy, sop[nn].Roz);
		pf_Obs = fopen(Obs_save_name[nn], "w");
		fclose(pf_Obs);
	}


	// 관찰평면의 갯수 및 관찰평면 지정
	//fscanf(pFdtd, "%d\n", &no_planes);
	no_planes = Jdocument["SIMULATION_ENVIORONMENT"]["Distribution"].Size();
	sopd      = Tptr1<struct_Obs_Plane>(no_planes);
	Plane_save_name = Tptr2<char>(no_planes, 256);
	char plane[3];
	char Save_Distribution_Path[1024];
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

			sprintf(Save_Distribution_Path, "./result/distribution/XY_Plane[z=%2.2e]", sopd[nn].r_th_real);
			_mkdir(Save_Distribution_Path);
		}
		if (strcmp(plane, "YZ") == 0)
		{
			sopd[nn].s_plane = 2;
			for (int ii = 0; ii < maxx; ii++)
				if (temp_coordx[ii] == sopd[nn].r_th_real)
					sopd[nn].r_th = ii;

			sprintf(Save_Distribution_Path, "./result/distribution/YZ_Plane[x=%2.2e]", sopd[nn].r_th_real);
			_mkdir(Save_Distribution_Path);
		}
		if (strcmp(plane, "XZ") == 0)
		{
			sopd[nn].s_plane = 3;
			for (int ii = 0; ii < maxy; ii++)
				if (temp_coordy[ii] == sopd[nn].r_th_real)
					sopd[nn].r_th = ii;

			sprintf(Save_Distribution_Path, "./result/distribution/XZ_Plane[y=%2.2e]", sopd[nn].r_th_real);
			_mkdir(Save_Distribution_Path);
		}
	}

	memfree(&temp_coordx);
	memfree(&temp_coordy);
	memfree(&temp_coordz);
	ifs.close();
}

void TSaveSnapshot::Field_Obs_Points()
{
	ss_loc = Nststep_iter % Nststep_save; // 0 ~ Nststep_save-1
	time_obs[ss_loc] = dt*Nststep_iter;
	for (int ii = 0; ii < no_points; ii++)
	{
		// need to define Filename_pfField
		ExV = (EOx[sop[ii].Ox    ][sop[ii].Oy    ][sop[ii].Oz    ] +
		       EOx[sop[ii].Ox - 1][sop[ii].Oy    ][sop[ii].Oz    ]) / 2;
		EyV = (EOy[sop[ii].Ox    ][sop[ii].Oy    ][sop[ii].Oz    ] +
		   	   EOy[sop[ii].Ox    ][sop[ii].Oy - 1][sop[ii].Oz    ]) / 2;
		EzV = (EOz[sop[ii].Ox    ][sop[ii].Oy    ][sop[ii].Oz    ] +
		   	   EOz[sop[ii].Ox    ][sop[ii].Oy    ][sop[ii].Oz - 1]) / 2;
		HxV = (HOx[sop[ii].Ox    ][sop[ii].Oy    ][sop[ii].Oz    ] +
		   	   HOx[sop[ii].Ox    ][sop[ii].Oy + 1][sop[ii].Oz    ] +
		   	   HOx[sop[ii].Ox    ][sop[ii].Oy    ][sop[ii].Oz + 1] +
		   	   HOx[sop[ii].Ox    ][sop[ii].Oy + 1][sop[ii].Oz + 1]) / 4;
		HyV = (HOy[sop[ii].Ox    ][sop[ii].Oy    ][sop[ii].Oz    ] +
		   	   HOy[sop[ii].Ox + 1][sop[ii].Oy    ][sop[ii].Oz    ] +
		   	   HOy[sop[ii].Ox    ][sop[ii].Oy    ][sop[ii].Oz + 1] +
		   	   HOy[sop[ii].Ox + 1][sop[ii].Oy    ][sop[ii].Oz + 1]) / 4;
		HzV = (HOz[sop[ii].Ox    ][sop[ii].Oy    ][sop[ii].Oz    ] +
		   	   HOz[sop[ii].Ox + 1][sop[ii].Oy    ][sop[ii].Oz    ] +
		   	   HOz[sop[ii].Ox    ][sop[ii].Oy + 1][sop[ii].Oz    ] +
		   	   HOz[sop[ii].Ox + 1][sop[ii].Oy + 1][sop[ii].Oz    ]) / 4;

		sop[ii].Ex_save_array[ss_loc] = ExV;
		sop[ii].Ey_save_array[ss_loc] = EyV;
		sop[ii].Ez_save_array[ss_loc] = EzV;			   
		sop[ii].Hx_save_array[ss_loc] = HxV;
		sop[ii].Hy_save_array[ss_loc] = HyV;
		sop[ii].Hz_save_array[ss_loc] = HzV;
	}

	if (Nststep_iter % Nststep_save == Nststep_save-1)
	{
		for (int ii = 0; ii < no_points; ii++)
		{
			FILE *pf_Obs;
			pf_Obs = fopen(Obs_save_name[ii], "a"); // append, write, field point.
			for (int jj = 0; jj < Nststep_save; jj++)
			{
				fprintf(pf_Obs, "%e\t%e\t%e\t%e\t%e\t%e\t%e\r\n", time_obs[jj], sop[ii].Ex_save_array[jj], sop[ii].Ey_save_array[jj], sop[ii].Ez_save_array[jj],
																		sop[ii].Hx_save_array[jj], sop[ii].Hy_save_array[jj], sop[ii].Hz_save_array[jj] );
			}
			fclose(pf_Obs);

			// 데이터 초기화
			
			memset(sop[ii].Ex_save_array, 0, sizeof(double)*Nststep_save);
			memset(sop[ii].Ey_save_array, 0, sizeof(double)*Nststep_save);
			memset(sop[ii].Ez_save_array, 0, sizeof(double)*Nststep_save);
			memset(sop[ii].Hx_save_array, 0, sizeof(double)*Nststep_save);
			memset(sop[ii].Hy_save_array, 0, sizeof(double)*Nststep_save);

		}
		memset(time_obs, 0, sizeof(double) * Nststep_save);
	}
	else if (Nststep_iter == Nststep_Max - 1)
	{
		for (int ii = 0; ii < no_points; ii++)
		{
			FILE *pf_Obs;
			pf_Obs = fopen(Obs_save_name[ii], "a"); // append, write, field point.
			for (int jj = 0; jj < ss_loc; jj++)
			{
				fprintf(pf_Obs, "%e\t%e\t%e\t%e\t%e\t%e\t%e\r\n", time_obs[jj], sop[ii].Ex_save_array[jj], sop[ii].Ey_save_array[jj], sop[ii].Ez_save_array[jj],
					sop[ii].Hx_save_array[jj], sop[ii].Hy_save_array[jj], sop[ii].Hz_save_array[jj]);
			}
			fclose(pf_Obs);

			// 데이터 초기화
			memset(sop[ii].Ex_save_array, 0, sizeof(double)*Nststep_save);
			memset(sop[ii].Ey_save_array, 0, sizeof(double)*Nststep_save);
			memset(sop[ii].Ez_save_array, 0, sizeof(double)*Nststep_save);
			memset(sop[ii].Hx_save_array, 0, sizeof(double)*Nststep_save);
			memset(sop[ii].Hy_save_array, 0, sizeof(double)*Nststep_save);
		}
		memset(time_obs, 0, sizeof(double) * Nststep_save);
	}
}

void TSaveSnapshot::Save_Composition()
{
	ss_loc = Nststep_iter % Nststep_save; // 0 ~ Nststep_save-1
	time_comp[ss_loc] = dt * Nststep_iter;
	double comxy = txsize * tysize;
	for (int ii = 0; ii < no_composit; ii++)
	{
		sumEx = 0; 	sumEy = 0;
		for (int xx = 0; xx < txsize  ; xx++)
		for (int yy = 0; yy < tysize  ; yy++)
		{
			ExV = (EOx[xx][yy][pSOC[ii].Oz] +
				   EOx[xx][yy + 1][pSOC[ii].Oz]) / 2;
			sumEx += ExV;
		}
		for (int xx = 0; xx < txsize    ; xx++)
		for (int yy = 0; yy < tysize    ; yy++)
		{
			EyV = (EOy[xx][yy][pSOC[ii].Oz] +
				   EOy[xx + 1][yy][pSOC[ii].Oz]) / 2;
			sumEy += EyV;
		}

		pSOC[ii].SumEx[ss_loc] = sumEx/ comxy;
		pSOC[ii].SumEy[ss_loc] = sumEy/ comxy;
	}

	if (Nststep_iter % Nststep_save == Nststep_save - 1)
	{
		for (int ii = 0; ii < no_composit; ii++)
		{
			FILE *pf_Obs;
			pf_Obs = fopen(Comp_save_name[ii], "a"); // append, write, field point.
			for (int jj = 0; jj < Nststep_save; jj++)
			{
				fprintf(pf_Obs, "%le\t%le\t%le\r\n", time_comp[jj], pSOC[ii].SumEx[jj], pSOC[ii].SumEy[jj]);
			}
			fclose(pf_Obs);


			memset(pSOC[ii].SumEx, 0, sizeof(double)*Nststep_save);
			memset(pSOC[ii].SumEy, 0, sizeof(double)*Nststep_save);
		}
		memset(time_comp, 0, sizeof(double) * Nststep_save);
	}
	else if (Nststep_iter == Nststep_Max-1)
	{
		for (int ii = 0; ii < no_composit; ii++)
		{
			FILE *pf_Obs;
			pf_Obs = fopen(Comp_save_name[ii], "a"); // append, write, field point.
			for (int jj = 0; jj < ss_loc; jj++)
			{
				fprintf(pf_Obs, "%le\t%le\t%le\r\n", time_comp[jj], pSOC[ii].SumEx[jj], pSOC[ii].SumEy[jj]);
			}
			fclose(pf_Obs);
			memset(pSOC[ii].SumEx, 0, sizeof(double)*Nststep_save);
			memset(pSOC[ii].SumEy, 0, sizeof(double)*Nststep_save);
		}
	}
}

void TSaveSnapshot::XY_plane_plot(int kk)
{
	E_value_2D = Tptr2<double>(txsize, tysize);
	for (int ii = 0; ii < txsize; ii++)
	{
		for (int jj = 0; jj < tysize; jj++)
		{
			ExV = (EOx[ii    ][jj    ][sopd[kk].r_th    ] + 
			       EOx[ii    ][jj + 1][sopd[kk].r_th    ] +
			       EOx[ii    ][jj    ][sopd[kk].r_th + 1] +
			       EOx[ii    ][jj + 1][sopd[kk].r_th + 1]) / 4;
			EyV = (EOy[ii    ][jj    ][sopd[kk].r_th    ] +
				   EOy[ii + 1][jj    ][sopd[kk].r_th    ] +
				   EOy[ii    ][jj    ][sopd[kk].r_th + 1] +
				   EOy[ii + 1][jj    ][sopd[kk].r_th + 1]) / 4;
			EzV = (EOz[ii    ][jj    ][sopd[kk].r_th    ] +
				   EOz[ii + 1][jj    ][sopd[kk].r_th    ] +
				   EOz[ii    ][jj + 1][sopd[kk].r_th    ] +
				   EOz[ii + 1][jj + 1][sopd[kk].r_th    ]) / 4;

			HxV = (HOx[ii    ][jj    ][sopd[kk].r_th    ] +
				   HOx[ii + 1][jj    ][sopd[kk].r_th    ]) / 2;
			HyV = (HOy[ii    ][jj    ][sopd[kk].r_th    ] +
				   HOy[ii    ][jj + 1][sopd[kk].r_th    ]) / 2;
			HzV = (HOz[ii    ][jj    ][sopd[kk].r_th    ] +
				   HOz[ii    ][jj    ][sopd[kk].r_th + 1]) / 2;
		   E_value_2D[ii][jj] = sqrt(ExV*ExV + EyV*EyV + EzV*EzV);

		}
		fwrite(E_value_2D[ii], sizeof(double), tysize, pf_Plane);

	}
	memfree(&E_value_2D, txsize);
	fclose(pf_Plane);
}

void TSaveSnapshot::YZ_plane_plot(int ii)
{
	//  yz 평면
	E_value_2D = Tptr2<double>(tysize, tzsize);
	for (int jj = 0; jj < tysize; jj++)
	{
		for (int kk = 0; kk < tzsize; kk++)
		{

			ExV = (EOx[sopd[ii].r_th    ][jj    ][kk    ] +
				   EOx[sopd[ii].r_th    ][jj + 1][kk    ] +
				   EOx[sopd[ii].r_th    ][jj    ][kk + 1] +
				   EOx[sopd[ii].r_th    ][jj + 1][kk + 1]) / 4;
			EyV = (EOy[sopd[ii].r_th    ][jj    ][kk    ] +
				   EOy[sopd[ii].r_th + 1][jj    ][kk    ] +
				   EOy[sopd[ii].r_th    ][jj    ][kk + 1] +
				   EOy[sopd[ii].r_th + 1][jj    ][kk + 1]) / 4;
			EzV = (EOz[sopd[ii].r_th    ][jj    ][kk    ] +
				   EOz[sopd[ii].r_th + 1][jj    ][kk    ] +
				   EOz[sopd[ii].r_th    ][jj + 1][kk    ] +
				   EOz[sopd[ii].r_th + 1][jj + 1][kk    ]) / 4;

			HxV = (HOx[sopd[ii].r_th    ][jj    ][kk    ] +
				   HOx[sopd[ii].r_th + 1][jj    ][kk    ]) / 2;
			HyV = (HOy[sopd[ii].r_th    ][jj    ][kk    ] +
				   HOy[sopd[ii].r_th    ][jj + 1][kk    ]) / 2;
			HzV = (HOz[sopd[ii].r_th    ][jj    ][kk    ] +
				   HOz[sopd[ii].r_th    ][jj    ][kk + 1]) / 2;

		   E_value_2D[jj][kk] = sqrt(ExV*ExV + EyV*EyV + EzV*EzV);

		}
		fwrite(E_value_2D[jj], sizeof(double), tzsize, pf_Plane);

	}
	memfree(&E_value_2D, tysize);
	fclose(pf_Plane);
}

void TSaveSnapshot::XZ_plane_plot(int jj)
{
	//  zx 평면
	E_value_2D = Tptr2<double>(txsize, tzsize);
	for (int ii = 0; ii < txsize; ii++)
	{
		for (int kk = 0; kk < tzsize; kk++)
		{
			ExV = (EOx[ii    ][sopd[jj].r_th    ][kk    ] +
				   EOx[ii    ][sopd[jj].r_th + 1][kk    ] +
				   EOx[ii    ][sopd[jj].r_th    ][kk + 1] +
				   EOx[ii    ][sopd[jj].r_th + 1][kk + 1]) / 4;
			EyV = (EOy[ii    ][sopd[jj].r_th    ][kk    ] +
				   EOy[ii + 1][sopd[jj].r_th    ][kk    ] +
				   EOy[ii    ][sopd[jj].r_th    ][kk + 1] +
				   EOy[ii + 1][sopd[jj].r_th    ][kk + 1]) / 4;
			EzV = (EOz[ii    ][sopd[jj].r_th    ][kk    ] +
				   EOz[ii + 1][sopd[jj].r_th    ][kk    ] +
				   EOz[ii    ][sopd[jj].r_th + 1][kk    ] +
				   EOz[ii + 1][sopd[jj].r_th + 1][kk    ]) / 4;

			HxV = (HOx[ii    ][sopd[jj].r_th    ][kk    ] +
				   HOx[ii + 1][sopd[jj].r_th    ][kk    ]) / 2;
			HyV = (HOy[ii    ][sopd[jj].r_th    ][kk    ] +
				   HOy[ii    ][sopd[jj].r_th + 1][kk    ]) / 2;
			HzV = (HOz[ii    ][sopd[jj].r_th    ][kk    ] +
				   HOz[ii    ][sopd[jj].r_th    ][kk + 1]) / 2;

			   E_value_2D[ii][kk] = sqrt(ExV*ExV + EyV*EyV + EzV*EzV);

		}
		fwrite(E_value_2D[ii], sizeof(double), tzsize, pf_Plane);

	}
	memfree(&E_value_2D, txsize);
	fclose(pf_Plane);
}

void TSaveSnapshot::Field_Obs_Plane(int ss)
{
	for (int ii = 0; ii < no_planes; ii++)
	{
		switch (sopd[ii].s_plane)
		{
			//  xy 평면
		case 1:
			sprintf(Plane_save_name[ii], "./result/distribution/XY_Plane[z=%2.2e]/Field_distribution_%04d.dat", sopd[ii].r_th_real, ss);
			pf_Plane = fopen(Plane_save_name[ii], "wb");   // append, write, field point.
			XY_plane_plot(ii);
			Field_XY_2D_Each(ii, ss);
			break;
			//  yz 평면
		case 2:
			sprintf(Plane_save_name[ii], "./result/distribution/YZ_Plane[x=%2.2e]/Field_distribution_%04d.dat", sopd[ii].r_th_real, ss);
			pf_Plane = fopen(Plane_save_name[ii], "wb");   // append, write, field point.
			YZ_plane_plot(ii);
			Field_YZ_2D_Each(ii, ss);
			break;
			//  xz 평면
		case 3:
			sprintf(Plane_save_name[ii], "./result/distribution/XZ_Plane[y=%2.2e]/Field_distribution_%04d.dat", sopd[ii].r_th_real, ss);
			pf_Plane = fopen(Plane_save_name[ii], "wb");   // append, write, field point.
			XZ_plane_plot(ii);
			Field_XZ_2D_Each(ii, ss);
			break;
		}
	}
}

// XY plane
void TSaveSnapshot::Field_XY_2D_Each(int kk, int ss)
{
	double **E_value_2D_Ex = Tptr2<double>(txsize, tysize);
	double **E_value_2D_Ey = Tptr2<double>(txsize, tysize);
	double **E_value_2D_Ez = Tptr2<double>(txsize, tysize);

	double **H_value_2D_Hx = Tptr2<double>(txsize, tysize);
	double **H_value_2D_Hy = Tptr2<double>(txsize, tysize);
	double **H_value_2D_Hz = Tptr2<double>(txsize, tysize);

	FILE *pf_ex, *pf_ey, *pf_ez;
	FILE *pf_hx, *pf_hy, *pf_hz;

	char Plane_save_name_Ex[256], Plane_save_name_Hx[256];
	char Plane_save_name_Ey[256], Plane_save_name_Hy[256];
	char Plane_save_name_Ez[256], Plane_save_name_Hz[256];

	sprintf(Plane_save_name_Ex, "./result/distribution/XY_Plane[z=%2.2e]/Field_distribution_%04d_Ex.dat", sopd[kk].r_th_real, ss);
	sprintf(Plane_save_name_Ey, "./result/distribution/XY_Plane[z=%2.2e]/Field_distribution_%04d_Ey.dat", sopd[kk].r_th_real, ss);
	sprintf(Plane_save_name_Ez, "./result/distribution/XY_Plane[z=%2.2e]/Field_distribution_%04d_Ez.dat", sopd[kk].r_th_real, ss);																		
	sprintf(Plane_save_name_Hx, "./result/distribution/XY_Plane[z=%2.2e]/Field_distribution_%04d_Hx.dat", sopd[kk].r_th_real, ss);
	sprintf(Plane_save_name_Hy, "./result/distribution/XY_Plane[z=%2.2e]/Field_distribution_%04d_Hy.dat", sopd[kk].r_th_real, ss);
	sprintf(Plane_save_name_Hz, "./result/distribution/XY_Plane[z=%2.2e]/Field_distribution_%04d_Hz.dat", sopd[kk].r_th_real, ss);


	pf_ex = fopen(Plane_save_name_Ex, "wb"); pf_hx = fopen(Plane_save_name_Hx, "wb");
	pf_ey = fopen(Plane_save_name_Ey, "wb"); pf_hy = fopen(Plane_save_name_Hy, "wb");
	pf_ez = fopen(Plane_save_name_Ez, "wb"); pf_hz = fopen(Plane_save_name_Hz, "wb");

	for (int ii = 0; ii < txsize; ii++)
	{
		for (int jj = 0; jj < tysize; jj++)
		{
			ExV = (EOx[ii    ][jj    ][sopd[kk].r_th    ] +
				   EOx[ii    ][jj + 1][sopd[kk].r_th    ] +
				   EOx[ii    ][jj    ][sopd[kk].r_th + 1] +
				   EOx[ii    ][jj + 1][sopd[kk].r_th + 1]) / 4;
			EyV = (EOy[ii    ][jj    ][sopd[kk].r_th    ] +
				   EOy[ii + 1][jj    ][sopd[kk].r_th    ] +
				   EOy[ii    ][jj    ][sopd[kk].r_th + 1] +
				   EOy[ii + 1][jj    ][sopd[kk].r_th + 1]) / 4;
			EzV = (EOz[ii    ][jj    ][sopd[kk].r_th    ] +
				   EOz[ii + 1][jj    ][sopd[kk].r_th    ] +
				   EOz[ii    ][jj + 1][sopd[kk].r_th    ] +
				   EOz[ii + 1][jj + 1][sopd[kk].r_th    ]) / 4;

			HxV = (HOx[ii    ][jj    ][sopd[kk].r_th    ] +
				   HOx[ii + 1][jj    ][sopd[kk].r_th    ]) / 2;
			HyV = (HOy[ii    ][jj    ][sopd[kk].r_th    ] +
				   HOy[ii    ][jj + 1][sopd[kk].r_th    ]) / 2;
			HzV = (HOz[ii    ][jj    ][sopd[kk].r_th    ] +
				   HOz[ii    ][jj    ][sopd[kk].r_th + 1]) / 2;

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
void TSaveSnapshot::Field_YZ_2D_Each(int ii, int ss)
{

	double **E_value_2D_Ex = Tptr2<double>(tysize, tzsize);
	double **E_value_2D_Ey = Tptr2<double>(tysize, tzsize);
	double **E_value_2D_Ez = Tptr2<double>(tysize, tzsize);
	
	double **H_value_2D_Hx = Tptr2<double>(tysize, tzsize);
	double **H_value_2D_Hy = Tptr2<double>(tysize, tzsize);
	double **H_value_2D_Hz = Tptr2<double>(tysize, tzsize);

	FILE *pf_ex, *pf_ey, *pf_ez;
	FILE *pf_hx, *pf_hy, *pf_hz;

	char Plane_save_name_Ex[256], Plane_save_name_Hx[256];
	char Plane_save_name_Ey[256], Plane_save_name_Hy[256];
	char Plane_save_name_Ez[256], Plane_save_name_Hz[256];

	sprintf(Plane_save_name_Ex, "./result/distribution/YZ_Plane[x=%2.2e]/Field_distribution_%04d_Ex.dat", sopd[ii].r_th_real, ss);
	sprintf(Plane_save_name_Ey, "./result/distribution/YZ_Plane[x=%2.2e]/Field_distribution_%04d_Ey.dat", sopd[ii].r_th_real, ss);
	sprintf(Plane_save_name_Ez, "./result/distribution/YZ_Plane[x=%2.2e]/Field_distribution_%04d_Ez.dat", sopd[ii].r_th_real, ss);																											  
	sprintf(Plane_save_name_Hx, "./result/distribution/YZ_Plane[x=%2.2e]/Field_distribution_%04d_Hx.dat", sopd[ii].r_th_real, ss);
	sprintf(Plane_save_name_Hy, "./result/distribution/YZ_Plane[x=%2.2e]/Field_distribution_%04d_Hy.dat", sopd[ii].r_th_real, ss);
	sprintf(Plane_save_name_Hz, "./result/distribution/YZ_Plane[x=%2.2e]/Field_distribution_%04d_Hz.dat", sopd[ii].r_th_real, ss);

	pf_ex = fopen(Plane_save_name_Ex, "wb"); pf_hx = fopen(Plane_save_name_Hx, "wb");
	pf_ey = fopen(Plane_save_name_Ey, "wb"); pf_hy = fopen(Plane_save_name_Hy, "wb");
	pf_ez = fopen(Plane_save_name_Ez, "wb"); pf_hz = fopen(Plane_save_name_Hz, "wb");

	for (int jj = 0; jj < tysize; jj++)
	{
		for (int kk = 0; kk < tzsize; kk++)
		{
			ExV = (EOx[sopd[ii].r_th    ][jj    ][kk    ] +
				   EOx[sopd[ii].r_th    ][jj + 1][kk    ] +
				   EOx[sopd[ii].r_th    ][jj    ][kk + 1] + 
				   EOx[sopd[ii].r_th    ][jj + 1][kk + 1]) / 4;
			EyV = (EOy[sopd[ii].r_th    ][jj    ][kk    ] +
				   EOy[sopd[ii].r_th + 1][jj    ][kk    ] +
				   EOy[sopd[ii].r_th    ][jj    ][kk + 1] +
				   EOy[sopd[ii].r_th + 1][jj    ][kk + 1]) / 4;
			EzV = (EOz[sopd[ii].r_th    ][jj    ][kk    ] +
				   EOz[sopd[ii].r_th + 1][jj    ][kk    ] +
				   EOz[sopd[ii].r_th    ][jj + 1][kk    ] +
				   EOz[sopd[ii].r_th + 1][jj + 1][kk    ]) / 4;

			HxV = (HOx[sopd[ii].r_th    ][jj    ][kk    ] +
				   HOx[sopd[ii].r_th + 1][jj    ][kk    ]) / 2;
			HyV = (HOy[sopd[ii].r_th    ][jj    ][kk    ] +
				   HOy[sopd[ii].r_th    ][jj + 1][kk    ]) / 2;
			HzV = (HOz[sopd[ii].r_th    ][jj    ][kk    ] +
				   HOz[sopd[ii].r_th    ][jj    ][kk + 1]) / 2;

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

void TSaveSnapshot::Field_XZ_2D_Each(int jj, int ss)
{
	double **E_value_2D_Ex = Tptr2<double>(txsize, tzsize);
	double **E_value_2D_Ey = Tptr2<double>(txsize, tzsize);
	double **E_value_2D_Ez = Tptr2<double>(txsize, tzsize);
													
	double **H_value_2D_Hx = Tptr2<double>(txsize, tzsize);
	double **H_value_2D_Hy = Tptr2<double>(txsize, tzsize);
	double **H_value_2D_Hz = Tptr2<double>(txsize, tzsize);

	FILE *pf_ex, *pf_ey, *pf_ez;
	FILE *pf_hx, *pf_hy, *pf_hz;

	char Plane_save_name_Ex[256], Plane_save_name_Hx[256];
	char Plane_save_name_Ey[256], Plane_save_name_Hy[256];
	char Plane_save_name_Ez[256], Plane_save_name_Hz[256];

	sprintf(Plane_save_name_Ex, "./result/distribution/XZ_Plane[y=%2.2e]/Field_distribution_%04d_Ex.dat", sopd[jj].r_th_real, ss);
	sprintf(Plane_save_name_Ey, "./result/distribution/XZ_Plane[y=%2.2e]/Field_distribution_%04d_Ey.dat", sopd[jj].r_th_real, ss);
	sprintf(Plane_save_name_Ez, "./result/distribution/XZ_Plane[y=%2.2e]/Field_distribution_%04d_Ez.dat", sopd[jj].r_th_real, ss);
	sprintf(Plane_save_name_Hx, "./result/distribution/XZ_Plane[y=%2.2e]/Field_distribution_%04d_Hx.dat", sopd[jj].r_th_real, ss);
	sprintf(Plane_save_name_Hy, "./result/distribution/XZ_Plane[y=%2.2e]/Field_distribution_%04d_Hy.dat", sopd[jj].r_th_real, ss);
	sprintf(Plane_save_name_Hz, "./result/distribution/XZ_Plane[y=%2.2e]/Field_distribution_%04d_Hz.dat", sopd[jj].r_th_real, ss);

	pf_ex = fopen(Plane_save_name_Ex, "wb"); pf_hx = fopen(Plane_save_name_Hx, "wb");
	pf_ey = fopen(Plane_save_name_Ey, "wb"); pf_hy = fopen(Plane_save_name_Hy, "wb");
	pf_ez = fopen(Plane_save_name_Ez, "wb"); pf_hz = fopen(Plane_save_name_Hz, "wb");

	for (int ii = 0; ii < txsize; ii++)
	{
		for (int kk = 0; kk < tzsize; kk++)
		{
			ExV = (EOx[ii    ][sopd[jj].r_th    ][kk    ] +
				   EOx[ii    ][sopd[jj].r_th + 1][kk    ] +
				   EOx[ii    ][sopd[jj].r_th    ][kk + 1] +
				   EOx[ii    ][sopd[jj].r_th + 1][kk + 1]) / 4;
			EyV = (EOy[ii    ][sopd[jj].r_th    ][kk    ] +
				   EOy[ii + 1][sopd[jj].r_th    ][kk    ] +
				   EOy[ii    ][sopd[jj].r_th    ][kk + 1] +
				   EOy[ii + 1][sopd[jj].r_th    ][kk + 1]) / 4;
			EzV = (EOz[ii    ][sopd[jj].r_th    ][kk    ] +
				   EOz[ii + 1][sopd[jj].r_th    ][kk    ] +
				   EOz[ii    ][sopd[jj].r_th + 1][kk    ] +
				   EOz[ii + 1][sopd[jj].r_th + 1][kk    ]) / 4;

			HxV = (HOx[ii    ][sopd[jj].r_th    ][kk    ] +
				   HOx[ii + 1][sopd[jj].r_th    ][kk    ]) / 2;
			HyV = (HOy[ii    ][sopd[jj].r_th    ][kk    ] +
				   HOy[ii    ][sopd[jj].r_th + 1][kk    ]) / 2;
			HzV = (HOz[ii    ][sopd[jj].r_th    ][kk    ] +
				   HOz[ii    ][sopd[jj].r_th    ][kk + 1]) / 2;

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

void TSaveSnapshot::Field_3D_All()
{
	for (int ii = 0; ii < txsize; ii++)
		for (int jj = 0; jj < tysize; jj++)
		{
			for (int kk = 0; kk < tzsize; kk++)
			{
				ExV = (EOx[ii    ][jj    ][kk    ] + 
					   EOx[ii    ][jj + 1][kk    ] +
					   EOx[ii    ][jj    ][kk + 1] +
					   EOx[ii    ][jj + 1][kk + 1]) / 4;
				EyV = (EOy[ii    ][jj    ][kk    ] +
					   EOy[ii + 1][jj    ][kk    ] +
					   EOy[ii    ][jj    ][kk + 1] +
					   EOy[ii + 1][jj    ][kk + 1]) / 4;
				EzV = (EOz[ii    ][jj    ][kk    ] +
					   EOz[ii + 1][jj    ][kk    ] +
					   EOz[ii    ][jj + 1][kk    ] +
					   EOz[ii + 1][jj + 1][kk    ]) / 4;

				HxV = (HOx[ii    ][jj    ][kk    ] +
					   HOx[ii + 1][jj    ][kk    ]) / 2;
				HyV = (HOy[ii    ][jj    ][kk    ] +
					   HOy[ii    ][jj + 1][kk    ]) / 2;
				HzV = (HOz[ii    ][jj    ][kk    ] +
					   HOz[ii    ][jj    ][kk + 1]) / 2;

				E_value[ii][jj][kk] = sqrt(ExV*ExV + EyV*EyV + EzV*EzV);
			}
		}
}

// void TSaveSnapshot::Initiate_optical_option()
// {
// 	TExctract_chiral_surface();
// 	Assign_optical_chirality_plane();
// }

void TSaveSnapshot::set_Media(int ***Media)
{
	this->Media = Media;
}

////
////void TSaveSnapshot::TExctract_chiral_surface()
////{
////	int Air = 0;
////	int PDMS = 2;
////	int Au = 3;
////	Extract_Optical_Chirality_plane = Tptr2<int>(txsize, tysize);
////	// 2차원 변환을 할 평면 일기 (= composit과 동일 위치)
////	// Method 1.0 : PDMS와 Air, Au의 접합면 검출. (z 방향으로만)
////	for (int ii = 0; ii < txsize; ii++)
////	for (int jj = 0; jj < tysize; jj++)
////	for (int kk = 0; kk < tzsize-1; kk++)
////	{
////		if ( Media[ii][jj][kk] == PDMS)
////		{
////			// 카이랄 경계를 envelop.
////			if ( Media[ii][jj][kk+1] == Au || Media[ii][jj][kk + 1] == Air)
////				Extract_Optical_Chirality_plane[ii][jj] = kk + 1;
////		}
////		else
////		{
////			continue;
////		}
////	}
////}
////
////
////void TSaveSnapshot::Assign_optical_chirality_plane()
////{
////	// 5개의 면,
////	// 1. 카이랄 평면 2등분해서 처음 가운데 끝
////	// 2. 투과면 특정지점
////	// 3. 카이랄 경계면
////
////	// 주파수 샘플링: 0.1 THz - 3 THz 101개의 샘플링.
////
////	// 각 지점에서 Ex, Ey, Ez, Hx, Hy, Hz 추출. : 약 1.5 GB 
////	// 공간 index +1 지점까지 최대 3 GB 메모리 공간이 필요.
////
////	pOptical_plane = Tptr3<struct_optical_chirality>(txsize, tysize, no_composit + 1);
////	sampling = 101;
////	for (int ii = 0; ii < txsize; ii++)
////	for (int jj = 0; jj < tysize; jj++)
////	for (int kk = 0; kk < no_composit+1; kk++)
////	{
////		pOptical_plane[ii][jj][kk].Ex_save_array1_re = Tptr1<double>(sampling); pOptical_plane[ii][jj][kk].Ex_save_array1_im = Tptr1<double>(sampling);
////		pOptical_plane[ii][jj][kk].Ey_save_array1_re = Tptr1<double>(sampling);	pOptical_plane[ii][jj][kk].Ey_save_array1_im = Tptr1<double>(sampling);
////		pOptical_plane[ii][jj][kk].Ez_save_array1_re = Tptr1<double>(sampling);	pOptical_plane[ii][jj][kk].Ez_save_array1_im = Tptr1<double>(sampling);
////		pOptical_plane[ii][jj][kk].Hx_save_array1_re = Tptr1<double>(sampling);	pOptical_plane[ii][jj][kk].Hx_save_array1_im = Tptr1<double>(sampling);
////		pOptical_plane[ii][jj][kk].Hy_save_array1_re = Tptr1<double>(sampling);	pOptical_plane[ii][jj][kk].Hy_save_array1_im = Tptr1<double>(sampling);
////		pOptical_plane[ii][jj][kk].Hz_save_array1_re = Tptr1<double>(sampling);	pOptical_plane[ii][jj][kk].Hz_save_array1_im = Tptr1<double>(sampling);
////
////		pOptical_plane[ii][jj][kk].Ex_save_array2_re = Tptr1<double>(sampling); pOptical_plane[ii][jj][kk].Ex_save_array2_im = Tptr1<double>(sampling);
////		pOptical_plane[ii][jj][kk].Ey_save_array2_re = Tptr1<double>(sampling);	pOptical_plane[ii][jj][kk].Ey_save_array2_im = Tptr1<double>(sampling);
////		pOptical_plane[ii][jj][kk].Ez_save_array2_re = Tptr1<double>(sampling);	pOptical_plane[ii][jj][kk].Ez_save_array2_im = Tptr1<double>(sampling);
////		pOptical_plane[ii][jj][kk].Hx_save_array2_re = Tptr1<double>(sampling);	pOptical_plane[ii][jj][kk].Hx_save_array2_im = Tptr1<double>(sampling);
////		pOptical_plane[ii][jj][kk].Hy_save_array2_re = Tptr1<double>(sampling);	pOptical_plane[ii][jj][kk].Hy_save_array2_im = Tptr1<double>(sampling);
////		pOptical_plane[ii][jj][kk].Hz_save_array2_re = Tptr1<double>(sampling);	pOptical_plane[ii][jj][kk].Hz_save_array2_im = Tptr1<double>(sampling);
////	}
////	
////	pSave_Temp_re = Tptr1<double>(101);
////	pSave_Temp_im = Tptr1<double>(101);
////	pFourier = new TFourierTransform();
////	double sFre = 0.1e12;
////	double tFre = 3e12;
////	double df = (tFre - sFre) / sampling;
////	pFourier->setParameters(sFre, tFre, df, dt);
////}
////
////void TSaveSnapshot::Transform_xy_plane_dft(int xx, int yy)
////{
////#pragma omp parallel 
////	{
////		for (int nn = 0; nn < no_composit; nn++)
////		{
////			int zz = pSOC[nn].Oz;
////
////			Ex1_of = EOx[xx][yy][zz];
////			Ey1_of = EOy[xx][yy][zz];
////			Ez1_of = EOz[xx][yy][zz];
////			Hx1_of = HOx[xx][yy][zz];
////			Hy1_of = HOy[xx][yy][zz];
////			Hz1_of = HOz[xx][yy][zz];
////			Ex2_of = EOx[xx][yy][zz + 1];
////			Ey2_of = EOy[xx][yy][zz + 1];
////			Ez2_of = EOz[xx][yy][zz + 1];
////			Hx2_of = HOx[xx][yy][zz + 1];
////			Hy2_of = HOy[xx][yy][zz + 1];
////			Hz2_of = HOz[xx][yy][zz + 1];
////
////			pFourier->onthefly_dft(Ex1_of, Nststep_iter);
////			pSave_Temp_re = pFourier->Get_Frequency_Domain_Response_Real();
////			pSave_Temp_im = pFourier->Get_Frequency_Domain_Response_Imag();
////#pragma omp for schedule(dynamic)
////			for (int ff = 0; ff < sampling; ff++)
////			{
////				pOptical_plane[xx][yy][nn].Ex_save_array1_re[ff] = pOptical_plane[xx][yy][nn].Ex_save_array1_re[ff] + pSave_Temp_re[ff];
////				pOptical_plane[xx][yy][nn].Ex_save_array1_im[ff] = pOptical_plane[xx][yy][nn].Ex_save_array1_im[ff] + pSave_Temp_im[ff];
////			}
////			pFourier->onthefly_dft(Ey1_of, Nststep_iter);
////			pSave_Temp_re = pFourier->Get_Frequency_Domain_Response_Real();
////			pSave_Temp_im = pFourier->Get_Frequency_Domain_Response_Imag();
////#pragma omp for schedule(dynamic)
////			for (int ff = 0; ff < sampling; ff++)
////			{
////				pOptical_plane[xx][yy][nn].Ey_save_array1_re[ff] = pOptical_plane[xx][yy][nn].Ey_save_array1_re[ff] + pSave_Temp_re[ff];
////				pOptical_plane[xx][yy][nn].Ey_save_array1_im[ff] = pOptical_plane[xx][yy][nn].Ey_save_array1_im[ff] + pSave_Temp_im[ff];
////			}
////
////			pFourier->onthefly_dft(Ez1_of, Nststep_iter);
////			pSave_Temp_re = pFourier->Get_Frequency_Domain_Response_Real();
////			pSave_Temp_im = pFourier->Get_Frequency_Domain_Response_Imag();
////#pragma omp for schedule(dynamic)
////			for (int ff = 0; ff < sampling; ff++)
////			{
////				pOptical_plane[xx][yy][nn].Ez_save_array1_re[ff] = pOptical_plane[xx][yy][nn].Ez_save_array1_re[ff] + pSave_Temp_re[ff];
////				pOptical_plane[xx][yy][nn].Ez_save_array1_im[ff] = pOptical_plane[xx][yy][nn].Ez_save_array1_im[ff] + pSave_Temp_im[ff];
////			}
////
////			pFourier->onthefly_dft(Ex2_of, Nststep_iter);
////			pSave_Temp_re = pFourier->Get_Frequency_Domain_Response_Real();
////			pSave_Temp_im = pFourier->Get_Frequency_Domain_Response_Imag();
////#pragma omp for schedule(dynamic)
////			for (int ff = 0; ff < sampling; ff++)
////			{
////				pOptical_plane[xx][yy][nn].Ex_save_array2_re[ff] = pOptical_plane[xx][yy][nn].Ex_save_array2_re[ff] + pSave_Temp_re[ff];
////				pOptical_plane[xx][yy][nn].Ex_save_array2_im[ff] = pOptical_plane[xx][yy][nn].Ex_save_array2_im[ff] + pSave_Temp_im[ff];
////			}
////
////			pFourier->onthefly_dft(Ey2_of, Nststep_iter);
////			pSave_Temp_re = pFourier->Get_Frequency_Domain_Response_Real();
////			pSave_Temp_im = pFourier->Get_Frequency_Domain_Response_Imag();
////#pragma omp for schedule(dynamic)
////			for (int ff = 0; ff < sampling; ff++)
////			{
////				pOptical_plane[xx][yy][nn].Ey_save_array2_re[ff] = pOptical_plane[xx][yy][nn].Ey_save_array2_re[ff] + pSave_Temp_re[ff];
////				pOptical_plane[xx][yy][nn].Ey_save_array2_im[ff] = pOptical_plane[xx][yy][nn].Ey_save_array2_im[ff] + pSave_Temp_im[ff];
////			}
////
////			pFourier->onthefly_dft(Ez2_of, Nststep_iter);
////			pSave_Temp_re = pFourier->Get_Frequency_Domain_Response_Real();
////			pSave_Temp_im = pFourier->Get_Frequency_Domain_Response_Imag();
////#pragma omp for schedule(dynamic)
////			for (int ff = 0; ff < sampling; ff++)
////			{
////				pOptical_plane[xx][yy][nn].Ez_save_array2_re[ff] = pOptical_plane[xx][yy][nn].Ez_save_array2_re[ff] + pSave_Temp_re[ff];
////				pOptical_plane[xx][yy][nn].Ez_save_array2_im[ff] = pOptical_plane[xx][yy][nn].Ez_save_array2_im[ff] + pSave_Temp_im[ff];
////			}
////			//
////			pFourier->onthefly_dft(Hx1_of, Nststep_iter);
////			pSave_Temp_re = pFourier->Get_Frequency_Domain_Response_Real();
////			pSave_Temp_im = pFourier->Get_Frequency_Domain_Response_Imag();
////#pragma omp for schedule(dynamic)
////			for (int ff = 0; ff < sampling; ff++)
////			{
////				pOptical_plane[xx][yy][nn].Hx_save_array1_re[ff] = pOptical_plane[xx][yy][nn].Hx_save_array1_re[ff] + pSave_Temp_re[ff];
////				pOptical_plane[xx][yy][nn].Hx_save_array1_im[ff] = pOptical_plane[xx][yy][nn].Hx_save_array1_im[ff] + pSave_Temp_im[ff];
////			}
////
////			pFourier->onthefly_dft(Hy1_of, Nststep_iter);
////			pSave_Temp_re = pFourier->Get_Frequency_Domain_Response_Real();
////			pSave_Temp_im = pFourier->Get_Frequency_Domain_Response_Imag();
////#pragma omp for schedule(dynamic)
////			for (int ff = 0; ff < sampling; ff++)
////			{
////				pOptical_plane[xx][yy][nn].Hy_save_array1_re[ff] = pOptical_plane[xx][yy][nn].Hy_save_array1_re[ff] + pSave_Temp_re[ff];
////				pOptical_plane[xx][yy][nn].Hy_save_array1_im[ff] = pOptical_plane[xx][yy][nn].Hy_save_array1_im[ff] + pSave_Temp_im[ff];
////			}
////
////			pFourier->onthefly_dft(Hz1_of, Nststep_iter);
////			pSave_Temp_re = pFourier->Get_Frequency_Domain_Response_Real();
////			pSave_Temp_im = pFourier->Get_Frequency_Domain_Response_Imag();
////#pragma omp for schedule(dynamic)
////			for (int ff = 0; ff < sampling; ff++)
////			{
////				pOptical_plane[xx][yy][nn].Hz_save_array1_re[ff] = pOptical_plane[xx][yy][nn].Hz_save_array1_re[ff] + pSave_Temp_re[ff];
////				pOptical_plane[xx][yy][nn].Hz_save_array1_im[ff] = pOptical_plane[xx][yy][nn].Hz_save_array1_im[ff] + pSave_Temp_im[ff];
////			}
////
////			pFourier->onthefly_dft(Hx2_of, Nststep_iter);
////			pSave_Temp_re = pFourier->Get_Frequency_Domain_Response_Real();
////			pSave_Temp_im = pFourier->Get_Frequency_Domain_Response_Imag();
////#pragma omp for schedule(dynamic)
////			for (int ff = 0; ff < sampling; ff++)
////			{
////				pOptical_plane[xx][yy][nn].Hx_save_array2_re[ff] = pOptical_plane[xx][yy][nn].Hx_save_array2_re[ff] + pSave_Temp_re[ff];
////				pOptical_plane[xx][yy][nn].Hx_save_array2_im[ff] = pOptical_plane[xx][yy][nn].Hx_save_array2_im[ff] + pSave_Temp_im[ff];
////			}
////
////			pFourier->onthefly_dft(Hy2_of, Nststep_iter);
////			pSave_Temp_re = pFourier->Get_Frequency_Domain_Response_Real();
////			pSave_Temp_im = pFourier->Get_Frequency_Domain_Response_Imag();
////#pragma omp for schedule(dynamic)
////			for (int ff = 0; ff < sampling; ff++)
////			{
////				pOptical_plane[xx][yy][nn].Hy_save_array2_re[ff] = pOptical_plane[xx][yy][nn].Hy_save_array2_re[ff] + pSave_Temp_re[ff];
////				pOptical_plane[xx][yy][nn].Hy_save_array2_im[ff] = pOptical_plane[xx][yy][nn].Hy_save_array2_im[ff] + pSave_Temp_im[ff];
////			}
////
////			pFourier->onthefly_dft(Hz2_of, Nststep_iter);
////			pSave_Temp_re = pFourier->Get_Frequency_Domain_Response_Real();
////			pSave_Temp_im = pFourier->Get_Frequency_Domain_Response_Imag();
////#pragma omp for schedule(dynamic)
////			for (int ff = 0; ff < sampling; ff++)
////			{
////				pOptical_plane[xx][yy][nn].Hz_save_array2_re[ff] = pOptical_plane[xx][yy][nn].Hz_save_array2_re[ff] + pSave_Temp_re[ff];
////				pOptical_plane[xx][yy][nn].Hz_save_array2_im[ff] = pOptical_plane[xx][yy][nn].Hz_save_array2_im[ff] + pSave_Temp_im[ff];
////			}
////		}
////	}
////}
////
////void TSaveSnapshot::Transform_chiral_surface_dft(int xx, int yy)
////{
////#pragma omp parallel 
////	{
////		int zz = Extract_Optical_Chirality_plane[xx][yy];
////		int nn = no_composit;
////
////		Ex1_of = EOx[xx][yy][zz];  
////		Ey1_of = EOy[xx][yy][zz];  
////		Ez1_of = EOz[xx][yy][zz];  
////		Hx1_of = HOx[xx][yy][zz];  
////		Hy1_of = HOy[xx][yy][zz];  
////		Hz1_of = HOz[xx][yy][zz];  
////		Ex2_of = EOx[xx][yy][zz+1];
////		Ey2_of = EOy[xx][yy][zz+1];
////		Ez2_of = EOz[xx][yy][zz+1];
////		Hx2_of = HOx[xx][yy][zz+1];
////		Hy2_of = HOy[xx][yy][zz+1];
////		Hz2_of = HOz[xx][yy][zz+1];
////
////		pFourier->onthefly_dft(Ex1_of, Nststep_iter);
////		pSave_Temp_re = pFourier->Get_Frequency_Domain_Response_Real();
////		pSave_Temp_im = pFourier->Get_Frequency_Domain_Response_Imag();
////#pragma omp for schedule(dynamic)
////		for (int ff = 0; ff < sampling; ff++)
////		{
////			pOptical_plane[xx][yy][nn].Ex_save_array1_re[ff] = pOptical_plane[xx][yy][nn].Ex_save_array1_re[ff] + pSave_Temp_re[ff];
////			pOptical_plane[xx][yy][nn].Ex_save_array1_im[ff] = pOptical_plane[xx][yy][nn].Ex_save_array1_im[ff] + pSave_Temp_im[ff];
////		}
////		pFourier->onthefly_dft(Ey1_of, Nststep_iter);
////		pSave_Temp_re = pFourier->Get_Frequency_Domain_Response_Real();
////		pSave_Temp_im = pFourier->Get_Frequency_Domain_Response_Imag();
////#pragma omp for schedule(dynamic)
////		for (int ff = 0; ff < sampling; ff++)
////		{
////			pOptical_plane[xx][yy][nn].Ey_save_array1_re[ff] = pOptical_plane[xx][yy][nn].Ey_save_array1_re[ff] + pSave_Temp_re[ff];
////			pOptical_plane[xx][yy][nn].Ey_save_array1_im[ff] = pOptical_plane[xx][yy][nn].Ey_save_array1_im[ff] + pSave_Temp_im[ff];
////		}
////
////		pFourier->onthefly_dft(Ez1_of, Nststep_iter);
////		pSave_Temp_re = pFourier->Get_Frequency_Domain_Response_Real();
////		pSave_Temp_im = pFourier->Get_Frequency_Domain_Response_Imag();
////#pragma omp for schedule(dynamic)
////		for (int ff = 0; ff < sampling; ff++)
////		{
////			pOptical_plane[xx][yy][nn].Ez_save_array1_re[ff] = pOptical_plane[xx][yy][nn].Ez_save_array1_re[ff] + pSave_Temp_re[ff];
////			pOptical_plane[xx][yy][nn].Ez_save_array1_im[ff] = pOptical_plane[xx][yy][nn].Ez_save_array1_im[ff] + pSave_Temp_im[ff];
////		}
////
////		pFourier->onthefly_dft(Ex2_of, Nststep_iter);
////		pSave_Temp_re = pFourier->Get_Frequency_Domain_Response_Real();
////		pSave_Temp_im = pFourier->Get_Frequency_Domain_Response_Imag();
////#pragma omp for schedule(dynamic)
////		for (int ff = 0; ff < sampling; ff++)
////		{
////			pOptical_plane[xx][yy][nn].Ex_save_array2_re[ff] = pOptical_plane[xx][yy][nn].Ex_save_array2_re[ff] + pSave_Temp_re[ff];
////			pOptical_plane[xx][yy][nn].Ex_save_array2_im[ff] = pOptical_plane[xx][yy][nn].Ex_save_array2_im[ff] + pSave_Temp_im[ff];
////		}
////
////		pFourier->onthefly_dft(Ey2_of, Nststep_iter);
////		pSave_Temp_re = pFourier->Get_Frequency_Domain_Response_Real();
////		pSave_Temp_im = pFourier->Get_Frequency_Domain_Response_Imag();
////#pragma omp for schedule(dynamic)
////		for (int ff = 0; ff < sampling; ff++)
////		{
////			pOptical_plane[xx][yy][nn].Ey_save_array2_re[ff] = pOptical_plane[xx][yy][nn].Ey_save_array2_re[ff] + pSave_Temp_re[ff];
////			pOptical_plane[xx][yy][nn].Ey_save_array2_im[ff] = pOptical_plane[xx][yy][nn].Ey_save_array2_im[ff] + pSave_Temp_im[ff];
////		}
////
////		pFourier->onthefly_dft(Ez2_of, Nststep_iter);
////		pSave_Temp_re = pFourier->Get_Frequency_Domain_Response_Real();
////		pSave_Temp_im = pFourier->Get_Frequency_Domain_Response_Imag();
////#pragma omp for schedule(dynamic)
////		for (int ff = 0; ff < sampling; ff++)
////		{
////			pOptical_plane[xx][yy][nn].Ez_save_array2_re[ff] = pOptical_plane[xx][yy][nn].Ez_save_array2_re[ff] + pSave_Temp_re[ff];
////			pOptical_plane[xx][yy][nn].Ez_save_array2_im[ff] = pOptical_plane[xx][yy][nn].Ez_save_array2_im[ff] + pSave_Temp_im[ff];
////		}
////		//
////		pFourier->onthefly_dft(Hx1_of, Nststep_iter);
////		pSave_Temp_re = pFourier->Get_Frequency_Domain_Response_Real();
////		pSave_Temp_im = pFourier->Get_Frequency_Domain_Response_Imag();
////#pragma omp for schedule(dynamic)
////		for (int ff = 0; ff < sampling; ff++)
////		{
////			pOptical_plane[xx][yy][nn].Hx_save_array1_re[ff] = pOptical_plane[xx][yy][nn].Hx_save_array1_re[ff] + pSave_Temp_re[ff];
////			pOptical_plane[xx][yy][nn].Hx_save_array1_im[ff] = pOptical_plane[xx][yy][nn].Hx_save_array1_im[ff] + pSave_Temp_im[ff];
////		}
////
////		pFourier->onthefly_dft(Hy1_of, Nststep_iter);
////		pSave_Temp_re = pFourier->Get_Frequency_Domain_Response_Real();
////		pSave_Temp_im = pFourier->Get_Frequency_Domain_Response_Imag();
////#pragma omp for schedule(dynamic)
////		for (int ff = 0; ff < sampling; ff++)
////		{
////			pOptical_plane[xx][yy][nn].Hy_save_array1_re[ff] = pOptical_plane[xx][yy][nn].Hy_save_array1_re[ff] + pSave_Temp_re[ff];
////			pOptical_plane[xx][yy][nn].Hy_save_array1_im[ff] = pOptical_plane[xx][yy][nn].Hy_save_array1_im[ff] + pSave_Temp_im[ff];
////		}
////
////		pFourier->onthefly_dft(Hz1_of, Nststep_iter);
////		pSave_Temp_re = pFourier->Get_Frequency_Domain_Response_Real();
////		pSave_Temp_im = pFourier->Get_Frequency_Domain_Response_Imag();
////#pragma omp for schedule(dynamic)
////		for (int ff = 0; ff < sampling; ff++)
////		{
////			pOptical_plane[xx][yy][nn].Hz_save_array1_re[ff] = pOptical_plane[xx][yy][nn].Hz_save_array1_re[ff] + pSave_Temp_re[ff];
////			pOptical_plane[xx][yy][nn].Hz_save_array1_im[ff] = pOptical_plane[xx][yy][nn].Hz_save_array1_im[ff] + pSave_Temp_im[ff];
////		}
////
////		pFourier->onthefly_dft(Hx2_of, Nststep_iter);
////		pSave_Temp_re = pFourier->Get_Frequency_Domain_Response_Real();
////		pSave_Temp_im = pFourier->Get_Frequency_Domain_Response_Imag();
////#pragma omp for schedule(dynamic)
////		for (int ff = 0; ff < sampling; ff++)
////		{
////			pOptical_plane[xx][yy][nn].Hx_save_array2_re[ff] = pOptical_plane[xx][yy][nn].Hx_save_array2_re[ff] + pSave_Temp_re[ff];
////			pOptical_plane[xx][yy][nn].Hx_save_array2_im[ff] = pOptical_plane[xx][yy][nn].Hx_save_array2_im[ff] + pSave_Temp_im[ff];
////		}
////
////		pFourier->onthefly_dft(Hy2_of, Nststep_iter);
////		pSave_Temp_re = pFourier->Get_Frequency_Domain_Response_Real();
////		pSave_Temp_im = pFourier->Get_Frequency_Domain_Response_Imag();
////#pragma omp for schedule(dynamic)
////		for (int ff = 0; ff < sampling; ff++)
////		{
////			pOptical_plane[xx][yy][nn].Hy_save_array2_re[ff] = pOptical_plane[xx][yy][nn].Hy_save_array2_re[ff] + pSave_Temp_re[ff];
////			pOptical_plane[xx][yy][nn].Hy_save_array2_im[ff] = pOptical_plane[xx][yy][nn].Hy_save_array2_im[ff] + pSave_Temp_im[ff];
////		}
////
////		pFourier->onthefly_dft(Hz2_of, Nststep_iter);
////		pSave_Temp_re = pFourier->Get_Frequency_Domain_Response_Real();
////		pSave_Temp_im = pFourier->Get_Frequency_Domain_Response_Imag();
////#pragma omp for schedule(dynamic)
////		for (int ff = 0; ff < sampling; ff++)
////		{
////			pOptical_plane[xx][yy][nn].Hz_save_array2_re[ff] = pOptical_plane[xx][yy][nn].Hz_save_array2_re[ff] + pSave_Temp_re[ff];
////			pOptical_plane[xx][yy][nn].Hz_save_array2_im[ff] = pOptical_plane[xx][yy][nn].Hz_save_array2_im[ff] + pSave_Temp_im[ff];
////		}
////	}
////}
////
////void TSaveSnapshot::store_2D_field_frequency()
////{
////	for (int ii = 0; ii < txsize ; ii++)
////	for (int jj = 0; jj < tysize ; jj++)
////	{
////		Transform_xy_plane_dft(ii, jj);
////		Transform_chiral_surface_dft(ii, jj);
////	}
////}
////
////void TSaveSnapshot::save_optical_chirality(int ss)
////{
////	if (ss % 10000 == 0)
////	{
////		for (int kk = 0; kk < no_composit + 1; kk++)
////		{
////			if (kk < no_composit)
////			{
////				sprintf(Fname_Ex1_re, "./result/optical_chirality/Field_frequency_%d_Ex1_real.dat", pSOC[kk].Oz);
////				sprintf(Fname_Ey1_re, "./result/optical_chirality/Field_frequency_%d_Ey1_real.dat", pSOC[kk].Oz);
////				sprintf(Fname_Ez1_re, "./result/optical_chirality/Field_frequency_%d_Ez1_real.dat", pSOC[kk].Oz);
////				sprintf(Fname_Ex2_re, "./result/optical_chirality/Field_frequency_%d_Ex2_real.dat", pSOC[kk].Oz);
////				sprintf(Fname_Ey2_re, "./result/optical_chirality/Field_frequency_%d_Ey2_real.dat", pSOC[kk].Oz);
////				sprintf(Fname_Ez2_re, "./result/optical_chirality/Field_frequency_%d_Ez2_real.dat", pSOC[kk].Oz);
////				sprintf(Fname_Hx1_re, "./result/optical_chirality/Field_frequency_%d_Hx1_real.dat", pSOC[kk].Oz);
////				sprintf(Fname_Hy1_re, "./result/optical_chirality/Field_frequency_%d_Hy1_real.dat", pSOC[kk].Oz);
////				sprintf(Fname_Hz1_re, "./result/optical_chirality/Field_frequency_%d_Hz1_real.dat", pSOC[kk].Oz);
////				sprintf(Fname_Hx2_re, "./result/optical_chirality/Field_frequency_%d_Hx2_real.dat", pSOC[kk].Oz);
////				sprintf(Fname_Hy2_re, "./result/optical_chirality/Field_frequency_%d_Hy2_real.dat", pSOC[kk].Oz);
////				sprintf(Fname_Hz2_re, "./result/optical_chirality/Field_frequency_%d_Hz2_real.dat", pSOC[kk].Oz);
////
////				sprintf(Fname_Ex1_im, "./result/optical_chirality/Field_frequency_%d_Ex1_imag.dat", pSOC[kk].Oz);
////				sprintf(Fname_Ey1_im, "./result/optical_chirality/Field_frequency_%d_Ey1_imag.dat", pSOC[kk].Oz);
////				sprintf(Fname_Ez1_im, "./result/optical_chirality/Field_frequency_%d_Ez1_imag.dat", pSOC[kk].Oz);
////				sprintf(Fname_Ex2_im, "./result/optical_chirality/Field_frequency_%d_Ex2_imag.dat", pSOC[kk].Oz);
////				sprintf(Fname_Ey2_im, "./result/optical_chirality/Field_frequency_%d_Ey2_imag.dat", pSOC[kk].Oz);
////				sprintf(Fname_Ez2_im, "./result/optical_chirality/Field_frequency_%d_Ez2_imag.dat", pSOC[kk].Oz);
////				sprintf(Fname_Hx1_im, "./result/optical_chirality/Field_frequency_%d_Hx1_imag.dat", pSOC[kk].Oz);
////				sprintf(Fname_Hy1_im, "./result/optical_chirality/Field_frequency_%d_Hy1_imag.dat", pSOC[kk].Oz);
////				sprintf(Fname_Hz1_im, "./result/optical_chirality/Field_frequency_%d_Hz1_imag.dat", pSOC[kk].Oz);
////				sprintf(Fname_Hx2_im, "./result/optical_chirality/Field_frequency_%d_Hx2_imag.dat", pSOC[kk].Oz);
////				sprintf(Fname_Hy2_im, "./result/optical_chirality/Field_frequency_%d_Hy2_imag.dat", pSOC[kk].Oz);
////				sprintf(Fname_Hz2_im, "./result/optical_chirality/Field_frequency_%d_Hz2_imag.dat", pSOC[kk].Oz);
////
////			}
////			else
////			{
////				sprintf(Fname_Ex1_re, "./result/optical_chirality/Field_frequency_surface_Ex1_real.dat");
////				sprintf(Fname_Ey1_re, "./result/optical_chirality/Field_frequency_surface_Ey1_real.dat");
////				sprintf(Fname_Ez1_re, "./result/optical_chirality/Field_frequency_surface_Ez1_real.dat");
////				sprintf(Fname_Ex2_re, "./result/optical_chirality/Field_frequency_surface_Ex2_real.dat");
////				sprintf(Fname_Ey2_re, "./result/optical_chirality/Field_frequency_surface_Ey2_real.dat");
////				sprintf(Fname_Ez2_re, "./result/optical_chirality/Field_frequency_surface_Ez2_real.dat");
////				sprintf(Fname_Hx1_re, "./result/optical_chirality/Field_frequency_surface_Hx1_real.dat");
////				sprintf(Fname_Hy1_re, "./result/optical_chirality/Field_frequency_surface_Hy1_real.dat");
////				sprintf(Fname_Hz1_re, "./result/optical_chirality/Field_frequency_surface_Hz1_real.dat");
////				sprintf(Fname_Hx2_re, "./result/optical_chirality/Field_frequency_surface_Hx2_real.dat");
////				sprintf(Fname_Hy2_re, "./result/optical_chirality/Field_frequency_surface_Hy2_real.dat");
////				sprintf(Fname_Hz2_re, "./result/optical_chirality/Field_frequency_surface_Hz2_real.dat");
////
////				sprintf(Fname_Ex1_im, "./result/optical_chirality/Field_frequency_surface_Ex1_imag.dat");
////				sprintf(Fname_Ey1_im, "./result/optical_chirality/Field_frequency_surface_Ey1_imag.dat");
////				sprintf(Fname_Ez1_im, "./result/optical_chirality/Field_frequency_surface_Ez1_imag.dat");
////				sprintf(Fname_Ex2_im, "./result/optical_chirality/Field_frequency_surface_Ex2_imag.dat");
////				sprintf(Fname_Ey2_im, "./result/optical_chirality/Field_frequency_surface_Ey2_imag.dat");
////				sprintf(Fname_Ez2_im, "./result/optical_chirality/Field_frequency_surface_Ez2_imag.dat");
////				sprintf(Fname_Hx1_im, "./result/optical_chirality/Field_frequency_surface_Hx1_imag.dat");
////				sprintf(Fname_Hy1_im, "./result/optical_chirality/Field_frequency_surface_Hy1_imag.dat");
////				sprintf(Fname_Hz1_im, "./result/optical_chirality/Field_frequency_surface_Hz1_imag.dat");
////				sprintf(Fname_Hx2_im, "./result/optical_chirality/Field_frequency_surface_Hx2_imag.dat");
////				sprintf(Fname_Hy2_im, "./result/optical_chirality/Field_frequency_surface_Hy2_imag.dat");
////				sprintf(Fname_Hz2_im, "./result/optical_chirality/Field_frequency_surface_Hz2_imag.dat");
////			}
////
////			pF_Ex1_re = fopen(Fname_Ex1_re, "wb"); pF_Ex1_im = fopen(Fname_Ex1_im, "wb");
////			pF_Ey1_re = fopen(Fname_Ey1_re, "wb"); pF_Ey1_im = fopen(Fname_Ey1_im, "wb");
////			pF_Ez1_re = fopen(Fname_Ez1_re, "wb"); pF_Ez1_im = fopen(Fname_Ez1_im, "wb");
////			pF_Ex2_re = fopen(Fname_Ex2_re, "wb"); pF_Ex2_im = fopen(Fname_Ex2_im, "wb");
////			pF_Ey2_re = fopen(Fname_Ey2_re, "wb"); pF_Ey2_im = fopen(Fname_Ey2_im, "wb");
////			pF_Ez2_re = fopen(Fname_Ez2_re, "wb"); pF_Ez2_im = fopen(Fname_Ez2_im, "wb");
////			pF_Hx1_re = fopen(Fname_Hx1_re, "wb"); pF_Hx1_im = fopen(Fname_Hx1_im, "wb");
////			pF_Hy1_re = fopen(Fname_Hy1_re, "wb"); pF_Hy1_im = fopen(Fname_Hy1_im, "wb");
////			pF_Hz1_re = fopen(Fname_Hz1_re, "wb"); pF_Hz1_im = fopen(Fname_Hz1_im, "wb");
////			pF_Hx2_re = fopen(Fname_Hx2_re, "wb"); pF_Hx2_im = fopen(Fname_Hx2_im, "wb");
////			pF_Hy2_re = fopen(Fname_Hy2_re, "wb"); pF_Hy2_im = fopen(Fname_Hy2_im, "wb");
////			pF_Hz2_re = fopen(Fname_Hz2_re, "wb"); pF_Hz2_im = fopen(Fname_Hz2_im, "wb");
////
////			for (int ii = 0; ii < txsize; ii++)
////				for (int jj = 0; jj < tysize; jj++)
////				{
////					fwrite(pOptical_plane[ii][jj][kk].Ex_save_array1_re, sampling, sizeof(double), pF_Ex1_re); fwrite(pOptical_plane[ii][jj][kk].Ex_save_array1_im, sampling, sizeof(double), pF_Ex1_im);
////					fwrite(pOptical_plane[ii][jj][kk].Ey_save_array1_re, sampling, sizeof(double), pF_Ey1_re); fwrite(pOptical_plane[ii][jj][kk].Ey_save_array1_im, sampling, sizeof(double), pF_Ey1_im);
////					fwrite(pOptical_plane[ii][jj][kk].Ez_save_array1_re, sampling, sizeof(double), pF_Ez1_re); fwrite(pOptical_plane[ii][jj][kk].Ez_save_array1_im, sampling, sizeof(double), pF_Ez1_im);
////					fwrite(pOptical_plane[ii][jj][kk].Hx_save_array1_re, sampling, sizeof(double), pF_Hx1_re); fwrite(pOptical_plane[ii][jj][kk].Hx_save_array1_im, sampling, sizeof(double), pF_Hx1_im);
////					fwrite(pOptical_plane[ii][jj][kk].Hy_save_array1_re, sampling, sizeof(double), pF_Hy1_re); fwrite(pOptical_plane[ii][jj][kk].Hy_save_array1_im, sampling, sizeof(double), pF_Hy1_im);
////					fwrite(pOptical_plane[ii][jj][kk].Hz_save_array1_re, sampling, sizeof(double), pF_Hz1_re); fwrite(pOptical_plane[ii][jj][kk].Hz_save_array1_im, sampling, sizeof(double), pF_Hz1_im);
////
////					fwrite(pOptical_plane[ii][jj][kk].Ex_save_array2_re, sampling, sizeof(double), pF_Ex2_re); fwrite(pOptical_plane[ii][jj][kk].Ex_save_array2_im, sampling, sizeof(double), pF_Ex2_im);
////					fwrite(pOptical_plane[ii][jj][kk].Ey_save_array2_re, sampling, sizeof(double), pF_Ey2_re); fwrite(pOptical_plane[ii][jj][kk].Ey_save_array2_im, sampling, sizeof(double), pF_Ey2_im);
////					fwrite(pOptical_plane[ii][jj][kk].Ez_save_array2_re, sampling, sizeof(double), pF_Ez2_re); fwrite(pOptical_plane[ii][jj][kk].Ez_save_array2_im, sampling, sizeof(double), pF_Ez2_im);
////					fwrite(pOptical_plane[ii][jj][kk].Hx_save_array2_re, sampling, sizeof(double), pF_Hx2_re); fwrite(pOptical_plane[ii][jj][kk].Hx_save_array2_im, sampling, sizeof(double), pF_Hx2_im);
////					fwrite(pOptical_plane[ii][jj][kk].Hy_save_array2_re, sampling, sizeof(double), pF_Hy2_re); fwrite(pOptical_plane[ii][jj][kk].Hy_save_array2_im, sampling, sizeof(double), pF_Hy2_im);
////					fwrite(pOptical_plane[ii][jj][kk].Hz_save_array2_re, sampling, sizeof(double), pF_Hz2_re); fwrite(pOptical_plane[ii][jj][kk].Hz_save_array2_im, sampling, sizeof(double), pF_Hz2_im);
////
////				}
////
////			fclose(pF_Ex1_re); fclose(pF_Ex1_im); fclose(pF_Ex2_re); fclose(pF_Ex2_im);
////			fclose(pF_Ey1_re); fclose(pF_Ey1_im); fclose(pF_Ey2_re); fclose(pF_Ey2_im);
////			fclose(pF_Ez1_re); fclose(pF_Ez1_im); fclose(pF_Ez2_re); fclose(pF_Ez2_im);
////			fclose(pF_Hx1_re); fclose(pF_Hx1_im); fclose(pF_Hx2_re); fclose(pF_Hx2_im);
////			fclose(pF_Hy1_re); fclose(pF_Hy1_im); fclose(pF_Hy2_re); fclose(pF_Hy2_im);
////			fclose(pF_Hz1_re); fclose(pF_Hz1_im); fclose(pF_Hz2_re); fclose(pF_Hz2_im);
////
////		}
////	}
////}
