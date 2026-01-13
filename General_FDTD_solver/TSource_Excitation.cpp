//---------------------------------------------------------------------------
#pragma hdrstop
#include "TSource_Excitation.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)

//__fastcall TSource_Excitation::TSource_Excitation(double ***Ex, double ***Ey, double ***Ez, double ***Hx, double ***Hy, double ***Hz)
//{
//	this->Ex = Ex;  this->Ey = Ey;  this->Ez = Ez;
//	this->Hx = Hx;  this->Hy = Hy;  this->Hz = Hz;
//	alpha = 0.5;
//	E0 = 1;
//	pi = acos(-1.0);
//	eps0 = 1e-9 / 36 / pi;
//	strcpy(JsonFileName, "./IEM_solver_input_format.json");
//}

TSource_Excitation::TSource_Excitation(struct_Field_carrier *pFc)
{
	this->pFc = pFc;
	this->Ex = pFc->Ex_address;  this->Ey = pFc->Ey_address;  this->Ez = pFc->Ez_address;
	this->Hx = pFc->Hx_address;  this->Hy = pFc->Hy_address;  this->Hz = pFc->Hz_address;
	alpha = 3; //0.5;
	E0 = 1;
	pi = acos(-1.0);
	eps0 = 1e-9 / 36 / pi;
	strcpy(JsonFileName, "./IEM_solver_input_format.json");
}

TSource_Excitation::TSource_Excitation(){}
TSource_Excitation::~TSource_Excitation()
{
	memfree(&sst);
}

int FindIndex(double value, char Coordinate)
{
	char JsonFileName[256];
	strcpy(JsonFileName, "./IEM_solver_input_format.json");
	int index = 0;
	Document Jdocument;
	ifstream ifs(JsonFileName);
	IStreamWrapper isw(ifs);
	Jdocument.ParseStream(isw);

	double maxx, maxy, maxz;
	if (Coordinate == 'X')
	{
		maxx = Jdocument["MESH"]["Voxelsize"]["txsize"].GetDouble() + 1;
		for (int ii = 0; ii < maxx; ii++)
			if (value == Jdocument["MESH"]["CellLocation"]["xcoordinate"][ii].GetDouble())
				index = ii;
	}
	if (Coordinate == 'Y')
	{
		maxy = Jdocument["MESH"]["Voxelsize"]["tysize"].GetDouble() + 1;
		for (int ii = 0; ii < maxy; ii++)
			if (value == Jdocument["MESH"]["CellLocation"]["ycoordinate"][ii].GetDouble())
				index = ii;
	}
	if (Coordinate == 'Z')
	{
		maxz = Jdocument["MESH"]["Voxelsize"]["tzsize"].GetDouble() + 1;
		for (int ii = 0; ii < maxz; ii++)
			if (value == Jdocument["MESH"]["CellLocation"]["zcoordinate"][ii].GetDouble())
				index = ii;
	}

	ifs.close();
	return index;
}



void TSource_Excitation::Set_Pulse_parameters()
{
	Document Jdocument;
	ifstream ifs(JsonFileName);
	IStreamWrapper isw(ifs);
	Jdocument.ParseStream(isw);
	txsize = Jdocument["MESH"]["Voxelsize"]["txsize"].GetDouble();
	tysize = Jdocument["MESH"]["Voxelsize"]["tysize"].GetDouble();
	tzsize = Jdocument["MESH"]["Voxelsize"]["tzsize"].GetDouble();
	// Def Source
	/*x0 = Jdocument["Mesh"]["CellLocation"]["xcoordinate"][0].GetDouble();
	y0 = Jdocument["Mesh"]["CellLocation"]["ycoordinate"][0].GetDouble();
	z0 = Jdocument["Mesh"]["CellLocation"]["zcoordinate"][0].GetDouble();*/

	
	//초기화
	f0 = Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["f0"].GetDouble();
	bw = Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["BW"].GetDouble();
	pFc->f0 = f0;
	pFc->BW = bw;
	Period = Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["Period"].GetDouble();
	NSTEPS = (int)(Period / f0 / dt); // for wideband
	pFc->NSTEPS = NSTEPS;
	char Source_Type[1024];
	char PulseName[1024];
	strcpy(Source_Type, Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["TypeName"].GetString());
	strcpy(PulseName, Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["PulseName"].GetString());
	
	// int opt_source = 2;
	// double f0 = 2.55e12;
	// double fbw = 1.0;
	// Source Type
	if (strcmp(Source_Type, "EDGE SOURCE") == 0)
	{
		// sst[0] ~ sst [??] 까지 소스 이용 가능하도록 코드 작성
		printf("Source Type: EDGE SOURCE\n");
		SourceType = 2;
		no_source = Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["Location"].Size();
		sst = Tptr1<struct_Source_type>(no_source);
		Phase = Tptr1<double>(no_source);

		for (int ii = 0; ii < no_source; ii++)
		{
			sst[ii].Rsx1 = Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["Location"][ii]["Ini"][0].GetDouble();
			sst[ii].Rsx2 = Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["Location"][ii]["Ter"][0].GetDouble();
			sst[ii].Rsy1 = Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["Location"][ii]["Ini"][1].GetDouble();
			sst[ii].Rsy2 = Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["Location"][ii]["Ter"][1].GetDouble();
			sst[ii].Rsz1 = Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["Location"][ii]["Ini"][2].GetDouble();
			sst[ii].Rsz2 = Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["Location"][ii]["Ter"][2].GetDouble();
			// Findindex : 좌표->index값으로 변환
			sst[ii].Sx1 = FindIndex(sst[ii].Rsx1, 'X');
			sst[ii].Sx2 = FindIndex(sst[ii].Rsx2, 'X');
			sst[ii].Sy1 = FindIndex(sst[ii].Rsy1, 'Y');
			sst[ii].Sy2 = FindIndex(sst[ii].Rsy2, 'Y');
			sst[ii].Sz1 = FindIndex(sst[ii].Rsz1, 'Z');
			sst[ii].Sz2 = FindIndex(sst[ii].Rsz2, 'Z');
			Phase[ii] =pi/180*Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["Location"][ii]["Phase"].GetDouble();

			// edge 소스의 방향, 0: point, 1: x dir, 2: y dir, 3: z dir
			sst[ii].point_dir = 0;
			if (sst[ii].Sx1 - sst[ii].Sx2 != 0)
				sst[ii].point_dir = 1;
			if (sst[ii].Sy1 - sst[ii].Sy2 != 0)
				sst[ii].point_dir = 2;
			if (sst[ii].Sz1 - sst[ii].Sz2 != 0)
				sst[ii].point_dir = 3;

		}
	}
	if (strcmp(Source_Type, "LINE SOURCE") == 0) {} // edge의 연장선
	if (strcmp(Source_Type, "CURRENT SOURCE") == 0) { SourceType = 1; }
	if (strcmp(Source_Type, "PLANE WAVE") == 0)
	{
		printf("Source Type: PLANE WAVE TFSF\n");
		SourceType = 3;
		no_source = Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["Location"].Size();
		sst = Tptr1<struct_Source_type>(no_source);

		Phase = Tptr1<double>(no_source);
		sst[0].Rsx1 = Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["Location"][0]["Ini"][0].GetDouble();
		sst[0].Rsx2 = Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["Location"][0]["Ter"][0].GetDouble();
		sst[0].Rsy1 = Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["Location"][0]["Ini"][1].GetDouble();
		sst[0].Rsy2 = Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["Location"][0]["Ter"][1].GetDouble();
		sst[0].Rsz1 = Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["Location"][0]["Ini"][2].GetDouble();
		sst[0].Rsz2 = Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["Location"][0]["Ter"][2].GetDouble();
		
		sst[0].Sx1 = FindIndex(sst[0].Rsx1, 'X');
		sst[0].Sx2 = FindIndex(sst[0].Rsx2, 'X');
		sst[0].Sy1 = FindIndex(sst[0].Rsy1, 'Y');
		sst[0].Sy2 = FindIndex(sst[0].Rsy2, 'Y');
		sst[0].Sz1 = FindIndex(sst[0].Rsz1, 'Z');
		sst[0].Sz2 = FindIndex(sst[0].Rsz2, 'Z');

		printf("TFSF box index: Sx1: %d, Sx2: %d, Sy1: %d, Sy2: %d, Sz1: %d, Sz2: %d\n",
			sst[0].Sx1, sst[0].Sx2, sst[0].Sy1, sst[0].Sy2, sst[0].Sz1, sst[0].Sz2);

		
		PW_Theta = Tptr1<double>(no_source);
		PW_Fhi = Tptr1<double>(no_source);
		PW_Polar = Tptr1<double>(no_source);

		for (int ii = 0; ii < no_source; ii++)
		{
			Phase[ii]    = pi / 180.0 * Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["Location"][ii]["Phase"].GetDouble();
			PW_Theta[ii] = Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["Location"][ii]["Theta"].GetDouble();
			PW_Fhi[ii]   = Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["Location"][ii]["Fhi"].GetDouble();
			// 90 0
			PW_Polar[ii] = Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["Location"][ii]["Polar"].GetDouble();
		}

		printf("Plane wave direction: Theta: %f, Fhi: %f, Polarization: %f\n",PW_Theta[0], PW_Fhi[0], PW_Polar[0]);
		Set_TFSF();
	}
	if (strcmp(Source_Type, "HYBRID MOM") == 0) 
	{
		SourceType = 4;
		pHsi = Tptr1<struct_Hybrid_source_info>(1);
		pHsi[0].f0 = Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["f0"].GetDouble();
		pHsi[0].ds = Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["ds"].GetDouble();
		pHsi[0].Datalen = Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["Datalen"].GetDouble();
		pHsi[0].SPX = Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["SPX"].GetDouble();
		pHsi[0].SPY = Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["SPY"].GetDouble();
		pHsi[0].SPZ = Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["SPZ"].GetDouble();
		strcpy(pHsi[0].Fname_Xminus, Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["Xminus"].GetString());
		strcpy(pHsi[0].Fname_Xplus,  Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["Xplus"].GetString());
		strcpy(pHsi[0].Fname_Yminus, Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["Yminus"].GetString());
		strcpy(pHsi[0].Fname_Yplus,  Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["Yplus"].GetString());
		strcpy(pHsi[0].Fname_Zminus, Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["Zminus"].GetString());
		strcpy(pHsi[0].Fname_Zplus,  Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["Zplus"].GetString());

		Set_Hybrid_19();
	}
	if (strcmp(Source_Type, "HYBRID MOM TFSF") == 0)
	{
		// json에서 Pulsename은 Harmonic / gaussian modulated sinewave 둘 중 하나 지정해야 함.
		// 단일 source면 harmonic, 다중 source면 gaussian modulated sinewave로 지정.
		printf("source type: HYBRID MOM TFSF\n");
		SourceType = 5;

		no_source = Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["sourcefile"].Size();
		printf("no_source=%d\n", no_source);

		// Hybrid TFSF에선 좌표값을 직접 받아와서 인덱스 설정할 때 사용해야 함. source file에 좌표만 들어있기 때문.
		
		int tx = Jdocument["MESH"]["Voxelsize"]["txsize"].GetInt();
		int ty = Jdocument["MESH"]["Voxelsize"]["tysize"].GetInt();
		int tz = Jdocument["MESH"]["Voxelsize"]["tzsize"].GetInt();

		xcoordinate = Tptr1<double>(tx + 1);
		ycoordinate = Tptr1<double>(ty + 1);
		zcoordinate = Tptr1<double>(tz + 1);

		coordinate.xcoordinate = Tptr1<double>(tx + 1);
		coordinate.ycoordinate = Tptr1<double>(ty + 1);
		coordinate.zcoordinate = Tptr1<double>(tz + 1);

		for (int ii = 0; ii <= tx; ii++)
		{
			xcoordinate[ii] = Jdocument["MESH"]["CellLocation"]["xcoordinate"][ii].GetDouble();
			coordinate.xcoordinate[ii] = xcoordinate[ii];
		}

		for (int jj = 0; jj <= ty; jj++)
		{
			ycoordinate[jj] = Jdocument["MESH"]["CellLocation"]["ycoordinate"][jj].GetDouble();
			coordinate.ycoordinate[jj] = ycoordinate[jj];
		}

		for (int kk = 0; kk <= tz; kk++)
		{
			zcoordinate[kk] = Jdocument["MESH"]["CellLocation"]["zcoordinate"][kk].GetDouble();
			coordinate.zcoordinate[kk] = zcoordinate[kk];
		}

		coordinate.xStaggered = Tptr1<double>(tx);
		coordinate.yStaggered = Tptr1<double>(ty);
		coordinate.zStaggered = Tptr1<double>(tz);

		for (int ii = 0; ii < tx; ii++)
			coordinate.xStaggered[ii] = 0.5 * (xcoordinate[ii] + xcoordinate[ii + 1]);

		for (int jj = 0; jj < ty; jj++)
			coordinate.yStaggered[jj] = 0.5 * (ycoordinate[jj] + ycoordinate[jj + 1]);

		for (int kk = 0; kk < tz; kk++)
			coordinate.zStaggered[kk] = 0.5 * (zcoordinate[kk] + zcoordinate[kk + 1]);



		// 이놈들을 Hybrid class로 넘겨주면 된다. 모든 source가 동일한 좌표를 사용하기 때문에 처음 생성할 때 한 번만 넘겨주면 된다.

		pHybTFSF = Tptr1<struct_Hybrid_TFSF_source_info>(no_source);

		// 다중 Hybrid source 입사 구현을 위해 이렇게 수정.
		for (int nn = 0; nn < no_source; nn++)
		{
			double f00 = Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["sourcefile"][nn]["f00"].GetDouble();
			pHybTFSF[nn].f0 = f00;
			strcpy(pHybTFSF[nn].plane1_Ey, Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["sourcefile"][nn]["plane1_Ey"].GetString());
			strcpy(pHybTFSF[nn].plane1_Ez, Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["sourcefile"][nn]["plane1_Ez"].GetString());
			strcpy(pHybTFSF[nn].plane3_Ey, Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["sourcefile"][nn]["plane3_Ey"].GetString());
			strcpy(pHybTFSF[nn].plane3_Ez, Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["sourcefile"][nn]["plane3_Ez"].GetString());
			strcpy(pHybTFSF[nn].plane2_Ex, Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["sourcefile"][nn]["plane2_Ex"].GetString());
			strcpy(pHybTFSF[nn].plane2_Ez, Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["sourcefile"][nn]["plane2_Ez"].GetString());
			strcpy(pHybTFSF[nn].plane4_Ex, Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["sourcefile"][nn]["plane4_Ex"].GetString());
			strcpy(pHybTFSF[nn].plane4_Ez, Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["sourcefile"][nn]["plane4_Ez"].GetString());
			strcpy(pHybTFSF[nn].plane5_Ex, Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["sourcefile"][nn]["plane5_Ex"].GetString());
			strcpy(pHybTFSF[nn].plane5_Ey, Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["sourcefile"][nn]["plane5_Ey"].GetString());
			strcpy(pHybTFSF[nn].plane6_Ex, Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["sourcefile"][nn]["plane6_Ex"].GetString());
			strcpy(pHybTFSF[nn].plane6_Ey, Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["sourcefile"][nn]["plane6_Ey"].GetString());

			// H-fields
			strcpy(pHybTFSF[nn].plane1_Hy, Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["sourcefile"][nn]["plane1_Hy"].GetString());
			strcpy(pHybTFSF[nn].plane1_Hz, Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["sourcefile"][nn]["plane1_Hz"].GetString());
			strcpy(pHybTFSF[nn].plane3_Hy, Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["sourcefile"][nn]["plane3_Hy"].GetString());
			strcpy(pHybTFSF[nn].plane3_Hz, Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["sourcefile"][nn]["plane3_Hz"].GetString());
			strcpy(pHybTFSF[nn].plane2_Hx, Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["sourcefile"][nn]["plane2_Hx"].GetString());
			strcpy(pHybTFSF[nn].plane2_Hz, Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["sourcefile"][nn]["plane2_Hz"].GetString());
			strcpy(pHybTFSF[nn].plane4_Hx, Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["sourcefile"][nn]["plane4_Hx"].GetString());
			strcpy(pHybTFSF[nn].plane4_Hz, Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["sourcefile"][nn]["plane4_Hz"].GetString());
			strcpy(pHybTFSF[nn].plane5_Hx, Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["sourcefile"][nn]["plane5_Hx"].GetString());
			strcpy(pHybTFSF[nn].plane5_Hy, Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["sourcefile"][nn]["plane5_Hy"].GetString());
			strcpy(pHybTFSF[nn].plane6_Hx, Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["sourcefile"][nn]["plane6_Hx"].GetString());
			strcpy(pHybTFSF[nn].plane6_Hy, Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["sourcefile"][nn]["plane6_Hy"].GetString());
		}

		Set_Hybrid_TFSF();
	}

	// Pulse Name
	// 0 : sine
	// 1 : gaussian modulated sine
	// 2 : user define
	if (strcmp(PulseName, "HARMONIC") == 0) { Opt_Pulse = 0; }
	if (strcmp(PulseName, "GAUSSIAN MODULATED SINEWAVE") == 0) { Opt_Pulse = 1; }
	if (strcmp(PulseName, "USER DEFINE") == 0) 
	{
		Opt_Pulse = 2;
		Get_User_define_source();
	}

	if(bw != 0)
	{
		//T = sqrt(2.)*sqrt(log10(2.)) / (m_pi*bw);
		T_0 = 1 / f0;
		T = sqrt(2.)*sqrt(log(20.)) / (pi*bw);
	}
	n_T = 3;
}

void TSource_Excitation::Inject_Edge_Source()
{
	int sx, sy, sz;
	int x_dir, y_dir, z_dir;
	int sx1, sx2, sy1, sy2, sz1, sz2;

	for (int ii = 0; ii < no_source; ii++)
	{
		x_dir = sst[ii].Sx2 - sst[ii].Sx1;
		y_dir = sst[ii].Sy2 - sst[ii].Sy1;
		z_dir = sst[ii].Sz2 - sst[ii].Sz1;

		if (x_dir > 0) // x-dir 정순
		{
			sy = (sst[ii].Sy1 + sst[ii].Sy2) / 2;
			sz = (sst[ii].Sz1 + sst[ii].Sz2) / 2;
			for (int jj = sst[ii].Sx1; jj <= sst[ii].Sx2; jj++)
				Ex[jj][sy][sz] = Ex[jj][sy][sz] - dt / eps0 * Pulse / abs(x_dir);
		}
		else if (x_dir < 0) // x-dir 역순
		{
			sy = (sst[ii].Sy1 + sst[ii].Sy2) / 2;
			sz = (sst[ii].Sz1 + sst[ii].Sz2) / 2;
			for (int jj = sst[ii].Sx1; jj >= sst[ii].Sx2; jj--)
				Ex[jj][sy][sz] = Ex[jj][sy][sz] - dt / eps0 * Pulse / abs(x_dir);
		}
		if (y_dir > 0)
		{
			sx = (sst[ii].Sx1 + sst[ii].Sx2) / 2;
			sz = (sst[ii].Sz1 + sst[ii].Sz2) / 2;
			for (int jj = sst[ii].Sy1; jj <= sst[ii].Sy2; jj++)
				Ey[sx][jj][sz] = Ey[sx][jj][sz] - dt / eps0 * Pulse / abs(y_dir);

		}
		else if (y_dir < 0)
		{
			sx = (sst[ii].Sx1 + sst[ii].Sx2) / 2;
			sz = (sst[ii].Sz1 + sst[ii].Sz2) / 2;
			for (int jj = sst[ii].Sy1; jj >= sst[ii].Sy2; jj--)
				Ey[sx][jj][sz] = Ey[sx][jj][sz] - dt / eps0 * Pulse / abs(y_dir);
		}

		if (z_dir > 0)
		{
			sx = (sst[ii].Sx1 + sst[ii].Sx2) / 2;
			sy = (sst[ii].Sy1 + sst[ii].Sy2) / 2;
			for (int jj = sst[ii].Sz1; jj <= sst[ii].Sz2; jj++)
				Ez[sx][sy][jj] = Ez[sx][sy][jj] - dt / eps0 * Pulse / abs(z_dir);

		}
		else if (z_dir < 0)
		{
			sx = (sst[ii].Sx1 + sst[ii].Sx2) / 2;
			sy = (sst[ii].Sy1 + sst[ii].Sy2) / 2;
			for (int jj = sst[ii].Sz1; jj >= sst[ii].Sz2; jj--)
				Ez[sx][sy][jj] = Ez[sx][sy][jj] - dt / eps0 * Pulse / abs(z_dir);
		}
	}

	// Edge
	//for (int ii = 0; ii < no_source; ii++)
	//{
	//	if (sst[ii].point_dir == 0)
	//	{
	//		Ex[sst[ii].Sx1][sst[ii].Sy1][sst[ii].Sz1] = Ex[sst[ii].Sx1][sst[ii].Sy1][sst[ii].Sz1] - dt / eps0 * Pulse;
	//	}
	//	if (sst[ii].point_dir == 1)
	//	{
	//		for (int xx = sst[ii].Sx1; xx < sst[ii].Sx2; xx++)
	//			Ex[xx][sst[ii].Sy1][sst[ii].Sz1] = Ex[xx][sst[ii].Sy1][sst[ii].Sz1] - dt / eps0 * Pulse / (sst[ii].Sx2 - sst[ii].Sx1 + 1);
	//	}
	//	if (sst[ii].point_dir == 2)
	//	{
	//		for (int yy = sst[ii].Sy1; yy < sst[ii].Sy2; yy++)
	//			Ey[sst[ii].Sx1][yy][sst[ii].Sz1] = Ey[sst[ii].Sx1][yy][sst[ii].Sz1] - dt / eps0 * Pulse / (sst[ii].Sy2 - sst[ii].Sy1 + 1);
	//	}
	//	if (sst[ii].point_dir == 3)
	//	{
	//		for (int zz = sst[ii].Sz1; zz < sst[ii].Sy2; zz++)
	//			Ez[sst[ii].Sx1][sst[ii].Sy1][zz] = Ez[sst[ii].Sx1][sst[ii].Sy1][zz] - dt / eps0 * Pulse / (sst[ii].Sy2 - sst[ii].Sy1 + 1);
	//	}
	//}
}

void TSource_Excitation::Inject_Pointsource()
{
	int sx, sy, sz;
	for (int ii = 0; ii < no_source; ii++)
	{
		sx = sst[ii].Sx1;
		sy = sst[ii].Sy1;
		sz = sst[ii].Sz1;
		if (sst[ii].Sx2 == 1)
		{
			Ex[sx][sy][sz] = Ex[sx][sy][sz] - dt / eps0 * Pulse / 2;
			Ex[sx - 1][sy][sz] = Ex[sx - 1][sy][sz] - dt / eps0 * Pulse / 2;
		}
		if (sst[ii].Sy2 == 1)
		{
			Ey[sx][sy][sz] = Ey[sx][sy][sz] - dt / eps0 * Pulse / 2;
			Ey[sx][sy - 1][sz] = Ey[sx][sy - 1][sz] - dt / eps0 * Pulse / 2;
		}
		if (sst[ii].Sz2 == 1)
		{
			Ez[sx][sy][sz] = Ez[sx][sy][sz] - dt / eps0 * Pulse / 2;
			Ez[sx][sy][sz - 1] = Ez[sx][sy][sz - 1] - dt / eps0 * Pulse / 2;
		}
	}
}

void TSource_Excitation::Inject_SourceE(int ss)
{
	switch (SourceType)
	{
	case 0: // voltage source
		break;
	case 1: // current source
		break;
	case 2: // edge source
		Inject_Edge_Source();
		break;
	case 3: // Plane wave (TF/SF)
		for (int ii = 0; ii < no_source; ii++)
		{
			Pulse = Get_Pulse(ss, Phase[ii]);
			planewave[ii]->TFSF_Electric_Field_OMP(Pulse);
			planewave[ii]->DApply_OMP();
		}
		break;
	case 4: // Hybrid 19년도 구현
		pHybrid->set_Timestep(ss);
		pHybrid->Calculate_Incident_Fiend();
		pHybrid->Apply_sourceE();
		break; 
	case 5: // Hybrid TFSF method
		for (int ii = 0; ii < no_source; ii++)
		{
			pHybrid_TFSF[ii]->set_current_steps(ss);
			pHybrid_TFSF[ii]->DApply_OMP();
			pHybrid_TFSF[ii]->HApply_OMP();
		}
		break;
	case 99: // point source
		Inject_Pointsource();
		break;
	}
}

void TSource_Excitation::Inject_SourceH(int ss)
{
	switch (SourceType)
	{
	case 0: // voltage source
		break;
	case 1: // current source
		break;
	case 2: // edge source
		break;
	case 3: // Plane wave (TF/SF)
		for (int ii = 0; ii < no_source; ii++)
		{
			planewave[ii]->TFSF_Magnetic_Field_OMP();
			planewave[ii]->HApply_OMP();
		}
		break;
	case 4: // Hybrid 19년도 구현
		break;
	case 5: // Hybrid TFSF method
		break;
	case 99: // point source
		break;
	}
}

void TSource_Excitation::Inject_Source_add(int n)
{
	if (SourceType == 4)
		switch (n)
		{
		case 1: // zero in E
			pHybrid->Set_In_zeroE();
			break;
		case 2: // source for H and zero in H
			pHybrid->Apply_sourceH();
			pHybrid->Set_In_zeroE();
			break;
		}
	else
		return;
}

double TSource_Excitation::Get_Pulse(int ss, double phi)
{
	switch (Opt_Pulse)
	{
	//case 0: // Sinewave
	//	Pulse = E0*sin(2.* pi *f0*dt*ss + phi);
	//	break;
	case 0: // Raised-Cosine-Ramped (DC-offset) Sinwave
		if (ss*dt < alpha*T_0)
			Pulse = E0*0.5*(1. - cos((pi *f0*dt*ss) / alpha))*sin(2.* pi *f0*dt*ss + phi);
		else
			Pulse = E0*sin(2.* pi *f0*dt*ss + phi);
		break;
	case 1:  // Gaussian-modulated sinewave
		Pulse = E0*sin(2.* pi *f0*dt*ss + phi)*exp(-(dt*ss - n_T*T)*(dt*ss - n_T*T) / (T*T));
		break;
		
		//case 3: // Raised-Cosine-Ramped (DC-offset) Constant
		//	if (ss*dt < alpha*T_0)
		//		Pulse = E0*0.5*(1.-cos(m_pi*f0*dt*ss/alpha))*E0;
		//	else
		//		Pulse = E0;
		//	break;
		//case 5: // Diff Gaussian-modulated sinewave
		//	Pulse = -E0*sin(2.*m_pi*f0*dt*ss)*(dt*ss-n_T*T)/T*exp(-(dt*ss-n_T*T)*(dt*ss-n_T*T)/(T*T));
		//	 break;
		//case 6: // Diff Gaussian
		//	Pulse = -E0*(dt*ss-n_T*T)/T*exp(-(dt*ss-n_T*T)*(dt*ss-n_T*T)/(T*T));
		//	break;
		//case 7:    // cosine-modulated sinewave
		//	Pulse = E0*cos(2.*m_pi*0.2*f0*dt*ss)*sin(2.*m_pi*f0*dt*ss);
		//	break;
		//case 8: // cosine sinewave
		//	Pulse = E0*cos(2.*m_pi*0.2*f0*dt*ss);
		//	break;
		//case 9:  // Slow-Smooth-Switched Sinewave
		//	int mm=30.;
		//	if (ss*dt < mm*T_0)
		//	{
		//		double xon=1.0-(mm*T_0-ss*dt)/(mm*T_0);
		//		Pulse = E0*(10.*xon*xon*xon - 15.*xon*xon*xon*xon + 6.*xon*xon*xon*xon*xon)*sin(2.*m_pi*f0*dt*ss);
		//          }
		//	else
		//		Pulse = E0*sin(2.*m_pi*f0*dt*ss);
		//	break;
	}
	return Pulse;
}

void TSource_Excitation::Set_Hybrid_19() 
{
	pHybrid = new TFDTD_Hybrid_mod_with_MOM(Ex, Ey, Ez, Hx, Hy, Hz);
	pHybrid->set_dt(dt);
	pHybrid->insert_parameters(pHsi);
	pHybrid->set_surface_current_data();
}

void TSource_Excitation::Set_Hybrid_TFSF()
{
	printf("Hybrid TFSF Source Setting...\n");
	printf("no_source = %d\n", no_source);
	pHybrid_TFSF = new TFDTD_Hybrid_TFSF * [no_source];
	for (int ii = 0 ; ii < no_source; ii++)
		pHybrid_TFSF[ii] = new TFDTD_Hybrid_TFSF(pFc , coordinate);

	for (int ii = 0; ii < no_source; ii++)
		pHybrid_TFSF[ii]->set_Hybrid_TFSF_box(pHybTFSF);

	// TF/SF에선 박스가 사전에 할당?되기 때문에 할당된 박스를 넘겨주는 작업이 필요하다면, 여기서는 source파일에서 직접 가져와야 한다.
	// TF/SF 코드와는 다른 형태를 취해야 할 듯 보인다...
	// box 설정은 생성자에서 전부 처리하도록 하려고 한다...
	
}

void TSource_Excitation::Set_TFSF()
{
	PbufferX = 0;
	PbufferY = 0;
	PbufferZ = 0;
	// Plane wave
	planewave = new TFDTD_Boundary_TFSF * [no_source];
	for (int ii = 0; ii < no_source; ii++)
		planewave[ii] = new TFDTD_Boundary_TFSF(pFc);

	for (int ii = 0; ii < no_source; ii++)
	{
		//planewave2 = new TFDTD_Boundary_TFSF(Ex, Ey, Ez, Hx, Hy, Hz);
		int pmlwbuffer_x1, pmlwbuffer_x2, pmlwbuffer_y1, pmlwbuffer_y2, pmlwbuffer_z1, pmlwbuffer_z2;
		pmlwbuffer_x1 = sst[0].Sx1;// PMLX + PbufferX;
		pmlwbuffer_x2 = sst[0].Sx2; // txsize; //txsize-(PMLX+PbufferX);

		pmlwbuffer_y1 = sst[0].Sy1; // PMLY + PbufferY;
		pmlwbuffer_y2 = sst[0].Sy2; // tysize; //tysize-(PMLY+PbufferY); 

		pmlwbuffer_z1 = sst[0].Sz1;// 0; // PMLZ + PbufferZ;
		pmlwbuffer_z2 = sst[0].Sz2;// planewave_pos; // -(PMLZ + PbufferZ);

		//planewave[ii]->insert_parameters(txsize, tysize, tzsize, min1D(dx), min1D(dy), min1D(dz), dt, f0, NSTEPS);
		//	planewave2->insert_parameters(txsize, tysize, tzsize, min1D(dx), min1D(dy), dz[pmlwbuffer_z2], dt, f0, NSTEPS);

		planewave[ii]->insert_inc_angle(PW_Theta[ii], PW_Fhi[ii], PW_Polar[ii]);
		//	planewave2->insert_inc_angle(PW_Theta, PW_Fhi2, PW_Polar);

		planewave[ii]->insert_TFSF_Box(pmlwbuffer_x1, pmlwbuffer_x2, pmlwbuffer_y1, pmlwbuffer_y2, pmlwbuffer_z1, pmlwbuffer_z2);
		//	planewave2->insert_TFSF_Box(pmlwbuffer_x1, pmlwbuffer_x2, pmlwbuffer_y1, pmlwbuffer_y2, pmlwbuffer_z1, pmlwbuffer_z2);
	}

}

double TSource_Excitation::Set_User_define_source(int ss)
{
	if (ss < 1000)
		return 0;
	else if (ss >= 1000 && ss < saved_Pulse_length)
		return pPulse[ss - 1000];
	else if (ss >= saved_Pulse_length)
		return 0;
}

void TSource_Excitation::Get_User_define_source()
{
	FILE *pf_user_define;
	pf_user_define = fopen("Pulse_user.dat", "rb");
	fread(&saved_Pulse_length, sizeof(double), 1, pf_user_define);
	pPulse = Tptr1<double>(saved_Pulse_length);
	fread(pPulse, sizeof(double), saved_Pulse_length, pf_user_define);
	fclose(pf_user_define);
}