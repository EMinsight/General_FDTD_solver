#include "TSaveComposit.h"

TSaveComposit::TSaveComposit()
{
	_mkdir("result\\composit");
}

void TSaveComposit::Prepare2save()
{
	Read_save_info();
	Assing_save_length();
}

void TSaveComposit::Save(int Nststep_iter)
{
	ss_loc = Nststep_iter % Nststep_save; // 0 ~ Nststep_save-1
	time_comp[ss_loc] = dt * Nststep_iter;
	double comxy = txsize * tysize;
	for (int ii = 0; ii < no_composit; ii++)
	{
		sumEx = 0; 	sumEy = 0; sumEz = 0;
		sumHx = 0; 	sumHy = 0; sumHz = 0;
		for (int xx = 0; xx < txsize; xx++)
			for (int yy = 0; yy < tysize; yy++)
			{
				ExV = (EOx[xx][yy][pSOC[ii].Oz] +
					EOx[xx][yy + 1][pSOC[ii].Oz]) / 2;
				sumEx += ExV;

				EyV = (EOy[xx][yy][pSOC[ii].Oz] +
					EOy[xx + 1][yy][pSOC[ii].Oz]) / 2;
				sumEy += EyV;

				EzV = (EOz[xx][yy][pSOC[ii].Oz]);
				sumEz += EzV;

				HxV = (HOx[xx][yy][pSOC[ii].Oz] +
					HOx[xx + 1][yy][pSOC[ii].Oz]) / 2;
				sumHx += HxV;

				HyV = (HOy[xx][yy][pSOC[ii].Oz] +
					HOy[xx][yy + 1][pSOC[ii].Oz]) / 2;
				sumHy += HyV;

				HzV = HOz[xx][yy][pSOC[ii].Oz];
				sumHz += HzV;
			}

		pSOC[ii].SumEx[ss_loc] = (RotationZXZp[0][0] * sumEx + RotationZXZp[0][1] * sumEy + RotationZXZp[0][2] * sumEz) / comxy;
		pSOC[ii].SumEy[ss_loc] = (RotationZXZp[1][0] * sumEx + RotationZXZp[1][1] * sumEy + RotationZXZp[1][2] * sumEz) / comxy;
		pSOC[ii].SumEz[ss_loc] = (RotationZXZp[2][0] * sumEx + RotationZXZp[2][1] * sumEy + RotationZXZp[2][2] * sumEz) / comxy;

		pSOC[ii].SumHx[ss_loc] = (RotationZXZp[0][0] * sumHx + RotationZXZp[0][1] * sumHy + RotationZXZp[0][2] * sumHz) / comxy;
		pSOC[ii].SumHy[ss_loc] = (RotationZXZp[1][0] * sumHx + RotationZXZp[1][1] * sumHy + RotationZXZp[1][2] * sumHz) / comxy;
		pSOC[ii].SumHz[ss_loc] = (RotationZXZp[2][0] * sumHx + RotationZXZp[2][1] * sumHy + RotationZXZp[2][2] * sumHz) / comxy;

		//		pSOC[ii].SumEx[ss_loc] = sumEx / comxy;
		//		pSOC[ii].SumEy[ss_loc] = sumEy / comxy;
	}

	// PRB 논문, abcd mat. ---> e^(-jkl) Ex Ey Hx Hy  --> 저장 
	// Ez 
	if (Nststep_iter % Nststep_save == Nststep_save - 1)
	{
		for (int nn = 0; nn < no_composit; nn++)
		{
			FILE* pf_Obs;
			pf_Obs = fopen(Comp_save_name[nn], "a"); // append, write, field point.
			for (int Nsave = 0; Nsave < Nststep_save; Nsave++)
			{
				fprintf(pf_Obs, "%le\t%le\t%le\t%le\t%le\t%le\t%le\t\r\n",
					time_comp[Nsave], pSOC[nn].SumEx[Nsave], pSOC[nn].SumEy[Nsave], pSOC[nn].SumEz[Nsave], pSOC[nn].SumHx[Nsave], pSOC[nn].SumHy[Nsave], pSOC[nn].SumHz[Nsave]);
			}
			fclose(pf_Obs);
			memset(pSOC[nn].SumEx, 0, sizeof(double) * Nststep_save);
			memset(pSOC[nn].SumEy, 0, sizeof(double) * Nststep_save);
			memset(pSOC[nn].SumEz, 0, sizeof(double) * Nststep_save);
			memset(pSOC[nn].SumHx, 0, sizeof(double) * Nststep_save);
			memset(pSOC[nn].SumHy, 0, sizeof(double) * Nststep_save);
			memset(pSOC[nn].SumHz, 0, sizeof(double) * Nststep_save);
		}
		memset(time_comp, 0, sizeof(double) * Nststep_save);
	}
	else if (Nststep_iter == Nststep_Max - 1)
	{
		for (int nn = 0; nn < no_composit; nn++)
		{
			FILE* pf_Obs;
			pf_Obs = fopen(Comp_save_name[nn], "a"); // append, write, field point.
			for (int Nsave = 0; Nsave < ss_loc; Nsave++)
			{
				fprintf(pf_Obs, "%le\t%le\t%le\t%le\t%le\t%le\t%le\t\r\n",
					time_comp[Nsave], pSOC[nn].SumEx[Nsave], pSOC[nn].SumEy[Nsave], pSOC[nn].SumEz[Nsave], pSOC[nn].SumHx[Nsave], pSOC[nn].SumHy[Nsave], pSOC[nn].SumHz[Nsave]);
			}
			fclose(pf_Obs);
			memset(pSOC[nn].SumEx, 0, sizeof(double) * Nststep_save);
			memset(pSOC[nn].SumEy, 0, sizeof(double) * Nststep_save);
			memset(pSOC[nn].SumEz, 0, sizeof(double) * Nststep_save);
			memset(pSOC[nn].SumHx, 0, sizeof(double) * Nststep_save);
			memset(pSOC[nn].SumHy, 0, sizeof(double) * Nststep_save);
			memset(pSOC[nn].SumHz, 0, sizeof(double) * Nststep_save);
		}
	}
}

void TSaveComposit::Read_save_info()
{
	ifstream ifs(JsonFileName);
	IStreamWrapper isw(ifs);
	Jdocument.ParseStream(isw);
	
	// 추후 Solver 에서 Location 삭제
	// Json에서 Position추가.
	no_composit = Jdocument["SIMULATION_ENVIORONMENT"]["Composition"].Size();

	// 관찰점의 갯수 및 관찰점 지정
	pSOC = Tptr1<struct_Obs_Composit>(no_composit);
	Comp_save_name = Tptr2<char>(no_composit, 256);

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


	no_composit = Jdocument["SIMULATION_ENVIORONMENT"]["Composition"].Size();
	//fscanf(pFdtd, "%d\n", &no_composit);
	pSOC = Tptr1<struct_Obs_Composit>(no_composit);
	Comp_save_name = Tptr2<char>(no_composit, 256);

	double temp = 0;
	
	//if (strcmp("PLANE WAVE", Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["TypeName"].GetString()) == 0)
	//{
	//	double theta = Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["Location"][0]["Theta"].GetDouble();
	//	double fhi = Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["Location"][0]["Fhi"].GetDouble();
	//	double psi = Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["Location"][0]["Polar"].GetDouble();
	//	makeRotation(fhi, theta, psi);
	//}

	for (int nn = 0; nn < no_composit; nn++)
	{
		pSOC[nn].Roz = Jdocument["SIMULATION_ENVIORONMENT"]["Composition"][nn]["Location"].GetDouble();
		for (int ii = 0; ii < maxz; ii++)
			if (temp_coordz[ii] == pSOC[nn].Roz)
				pSOC[nn].Oz = ii;
		sprintf(Comp_save_name[nn], "./result/composit/Observation_[z=%2.2e]_field_plot.dat", pSOC[nn].Roz);
		pf_Obs = fopen(Comp_save_name[nn], "w");
		fclose(pf_Obs);
	}

	memfree(&temp_coordx);
	memfree(&temp_coordy);
	memfree(&temp_coordz);
	ifs.close();
}

void TSaveComposit::Assing_save_length()
{
	time_comp = Tptr1<double>(Nststep_save);
	for (int ii = 0; ii < no_composit; ii++)
	{
		pSOC[ii].SumEx = Tptr1<double>(Nststep_save);
		pSOC[ii].SumEy = Tptr1<double>(Nststep_save);
		pSOC[ii].SumEz = Tptr1<double>(Nststep_save);
		pSOC[ii].SumHx = Tptr1<double>(Nststep_save);
		pSOC[ii].SumHy = Tptr1<double>(Nststep_save);
		pSOC[ii].SumHz = Tptr1<double>(Nststep_save);
	}
}

