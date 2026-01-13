#include "TSaveObservation.h"

TSaveObservation::TSaveObservation()
{
	_mkdir("result\\observation");
}
void TSaveObservation::Prepare2save()
{
	Read_save_info();
	Assing_save_length();
}

void TSaveObservation::Assing_save_length()
{
	time_obs = Tptr1<double>(Nststep_save);
	for (int ii = 0; ii < no_points; ii++)
	{
		sop[ii].Ex_save_array = Tptr1<double>(Nststep_save);
		sop[ii].Ey_save_array = Tptr1<double>(Nststep_save);
		sop[ii].Ez_save_array = Tptr1<double>(Nststep_save);

		sop[ii].Hx_save_array = Tptr1<double>(Nststep_save);
		sop[ii].Hy_save_array = Tptr1<double>(Nststep_save);
		sop[ii].Hz_save_array = Tptr1<double>(Nststep_save);
	}
}

void TSaveObservation::Read_save_info()
{
	ifstream ifs(JsonFileName);
	IStreamWrapper isw(ifs);
	Jdocument.ParseStream(isw);

	no_points = Jdocument["SIMULATION_ENVIORONMENT"]["Observation"].Size();
	// 관찰점의 갯수 및 관찰점 지정
	sop = Tptr1<struct_Obs_Point>(no_points);
	Obs_save_name = Tptr2<char>(no_points, 256);
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

	memfree(&temp_coordx);
	memfree(&temp_coordy);
	memfree(&temp_coordz);
	ifs.close();
}

void TSaveObservation::Save(int Nststep_iter)
{
	ss_loc = Nststep_iter % Nststep_save; // 0 ~ Nststep_save-1
	time_obs[ss_loc] = dt * Nststep_iter;
	for (int ii = 0; ii < no_points; ii++)
	{
		// need to define Filename_pfField
		ExV = (EOx[sop[ii].Ox][sop[ii].Oy][sop[ii].Oz] +
			EOx[sop[ii].Ox - 1][sop[ii].Oy][sop[ii].Oz]) / 2;
		EyV = (EOy[sop[ii].Ox][sop[ii].Oy][sop[ii].Oz] +
			EOy[sop[ii].Ox][sop[ii].Oy - 1][sop[ii].Oz]) / 2;
		EzV = (EOz[sop[ii].Ox][sop[ii].Oy][sop[ii].Oz] +
			EOz[sop[ii].Ox][sop[ii].Oy][sop[ii].Oz - 1]) / 2;
		HxV = (HOx[sop[ii].Ox][sop[ii].Oy][sop[ii].Oz] +
			HOx[sop[ii].Ox][sop[ii].Oy + 1][sop[ii].Oz] +
			HOx[sop[ii].Ox][sop[ii].Oy][sop[ii].Oz + 1] +
			HOx[sop[ii].Ox][sop[ii].Oy + 1][sop[ii].Oz + 1]) / 4;
		HyV = (HOy[sop[ii].Ox][sop[ii].Oy][sop[ii].Oz] +
			HOy[sop[ii].Ox + 1][sop[ii].Oy][sop[ii].Oz] +
			HOy[sop[ii].Ox][sop[ii].Oy][sop[ii].Oz + 1] +
			HOy[sop[ii].Ox + 1][sop[ii].Oy][sop[ii].Oz + 1]) / 4;
		HzV = (HOz[sop[ii].Ox][sop[ii].Oy][sop[ii].Oz] +
			HOz[sop[ii].Ox + 1][sop[ii].Oy][sop[ii].Oz] +
			HOz[sop[ii].Ox][sop[ii].Oy + 1][sop[ii].Oz] +
			HOz[sop[ii].Ox + 1][sop[ii].Oy + 1][sop[ii].Oz]) / 4;

		sop[ii].Ex_save_array[ss_loc] = ExV;
		sop[ii].Ey_save_array[ss_loc] = EyV;
		sop[ii].Ez_save_array[ss_loc] = EzV;
		sop[ii].Hx_save_array[ss_loc] = HxV;
		sop[ii].Hy_save_array[ss_loc] = HyV;
		sop[ii].Hz_save_array[ss_loc] = HzV;
	}

	if (Nststep_iter % Nststep_save == Nststep_save - 1)
	{
		for (int ii = 0; ii < no_points; ii++)
		{
			FILE* pf_Obs;
			pf_Obs = fopen(Obs_save_name[ii], "a"); // append, write, field point.
			for (int jj = 0; jj < Nststep_save; jj++)
			{
				fprintf(pf_Obs, "%e\t%e\t%e\t%e\t%e\t%e\t%e\r\n", time_obs[jj], sop[ii].Ex_save_array[jj], sop[ii].Ey_save_array[jj], sop[ii].Ez_save_array[jj],
					sop[ii].Hx_save_array[jj], sop[ii].Hy_save_array[jj], sop[ii].Hz_save_array[jj]);
			}
			fclose(pf_Obs);

			// 데이터 초기화

			memset(sop[ii].Ex_save_array, 0, sizeof(double) * Nststep_save);
			memset(sop[ii].Ey_save_array, 0, sizeof(double) * Nststep_save);
			memset(sop[ii].Ez_save_array, 0, sizeof(double) * Nststep_save);
			memset(sop[ii].Hx_save_array, 0, sizeof(double) * Nststep_save);
			memset(sop[ii].Hy_save_array, 0, sizeof(double) * Nststep_save);
			memset(sop[ii].Hz_save_array, 0, sizeof(double) * Nststep_save);

		}
		memset(time_obs, 0, sizeof(double) * Nststep_save);
	}
	else if (Nststep_iter == Nststep_Max - 1)
	{
		for (int ii = 0; ii < no_points; ii++)
		{
			FILE* pf_Obs;
			pf_Obs = fopen(Obs_save_name[ii], "a"); // append, write, field point.
			for (int jj = 0; jj < ss_loc; jj++)
			{
				fprintf(pf_Obs, "%e\t%e\t%e\t%e\t%e\t%e\t%e\r\n", time_obs[jj], sop[ii].Ex_save_array[jj], sop[ii].Ey_save_array[jj], sop[ii].Ez_save_array[jj],
					sop[ii].Hx_save_array[jj], sop[ii].Hy_save_array[jj], sop[ii].Hz_save_array[jj]);
			}
			fclose(pf_Obs);

			// 데이터 초기화
			memset(sop[ii].Ex_save_array, 0, sizeof(double) * Nststep_save);
			memset(sop[ii].Ey_save_array, 0, sizeof(double) * Nststep_save);
			memset(sop[ii].Ez_save_array, 0, sizeof(double) * Nststep_save);
			memset(sop[ii].Hx_save_array, 0, sizeof(double) * Nststep_save);
			memset(sop[ii].Hy_save_array, 0, sizeof(double) * Nststep_save);
		}
		memset(time_obs, 0, sizeof(double) * Nststep_save);
	}
}