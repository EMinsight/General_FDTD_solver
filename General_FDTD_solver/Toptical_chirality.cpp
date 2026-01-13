#include "Toptical_chirality.h"
#include <math.h>
#include <direct.h>

Toptical_chirality::Toptical_chirality()
{
	pi = acos(-1.0);
	eps0 = 1e-9 / 36 / pi;
	mu0 = 4 * pi*1e-7;
	c0 = 1 / sqrt(eps0*mu0);
	
	_mkdir("result\\optical_chirality");

	FILE *pf1;
	pf1 = fopen("CHIRALPBC.dat", "r");
	fscanf(pf1, "%d,%d\n", &aa1m1, &aa1m2);
	fclose(pf1);
}

void Toptical_chirality::Save(int Nststep_iter)
{
	int steady_state_offset = 0;
	complex <double> optical_Ex, optical_Ey, optical_Ez;
	complex <double> optical_Hx, optical_Hy, optical_Hz;
	complex <double> temp;
	double f1, f2;
	if (Nststep_iter == save_sta)
	{
		tt1 = Nststep_iter;
		for (int ii = 0; ii < txsize; ii++)
			for (int jj = 0; jj < tysize; jj++)
			{
				Single_Freq_Amp_Phase1(ii, jj);
			}
	}
	if (Nststep_iter == save_end)
	{
		// optical chirality 계산 가능.
		tt2 = Nststep_iter;
		for (int ii = 0; ii < txsize; ii++)
			for (int jj = 0; jj < tysize; jj++)
			{
				Single_Freq_Amp_Phase2(ii, jj);
			}

		for (int nn = 0; nn < Layer; nn++)
		{
			// optical chirliaty 를 구하는 평면
			for (int ii = 0; ii < txsize; ii++)
				for (int jj = 0; jj < tysize; jj++)
				{
					// 면적에 대해 크기와 페이즈 추출 후 optical chirality 계산.
					// optical_Ex, optical_Ey, optical_Ez
					// Ex
					f1 = pOptical_plane[nn].Ex_save_array_t1[ii][jj];
					f2 = pOptical_plane[nn].Ex_save_array_t2[ii][jj];
					//Phase = atan2(cos(2 * pi*f0*dt*tt2) - f2 / f1, sin(2 * pi*f0*tt2));
					Phase = atan2((f2 * sin(2 * pi * f0 * dt * tt1) - f1 * sin(2 * pi * f0 * dt * tt2)),
						(f1 * cos(2 * pi * f0 * dt * tt2) - f2 * cos(2 * pi * f0 * dt * tt1)));
					Amp = f1 / sin(2 * pi * f0 * dt * tt1 + Phase);
					Amp_Phase_Tunning();
					optical_Ex = complex<double>(Amp * cos(Phase), Amp * sin(Phase));

					// Ey
					f1 = pOptical_plane[nn].Ey_save_array_t1[ii][jj];
					f2 = pOptical_plane[nn].Ey_save_array_t2[ii][jj];
					Phase = atan2((f2 * sin(2 * pi * f0 * dt * tt1) - f1 * sin(2 * pi * f0 * dt * tt2)),
						(f1 * cos(2 * pi * f0 * dt * tt2) - f2 * cos(2 * pi * f0 * dt * tt1)));
					Amp = f1 / sin(2 * pi * f0 * dt * tt1 + Phase);
					Amp_Phase_Tunning();
					optical_Ey = complex<double>(Amp * cos(Phase), Amp * sin(Phase));

					// Ez
					f1 = pOptical_plane[nn].Ez_save_array_t1[ii][jj];
					f2 = pOptical_plane[nn].Ez_save_array_t2[ii][jj];
					Phase = atan2((f2 * sin(2 * pi * f0 * dt * tt1) - f1 * sin(2 * pi * f0 * dt * tt2)),
						(f1 * cos(2 * pi * f0 * dt * tt2) - f2 * cos(2 * pi * f0 * dt * tt1)));
					Amp = f1 / sin(2 * pi * f0 * dt * tt1 + Phase);
					Amp_Phase_Tunning();
					optical_Ez = complex<double>(Amp * cos(Phase), Amp * sin(Phase));

					// Hx
					f1 = pOptical_plane[nn].Hx_save_array_t1[ii][jj];
					f2 = pOptical_plane[nn].Hx_save_array_t2[ii][jj];
					Phase = atan2((f2 * sin(2 * pi * f0 * dt * tt1) - f1 * sin(2 * pi * f0 * dt * tt2)),
						(f1 * cos(2 * pi * f0 * dt * tt2) - f2 * cos(2 * pi * f0 * dt * tt1)));
					Amp = f1 / sin(2 * pi * f0 * dt * tt1 + Phase);
					Amp_Phase_Tunning();
					optical_Hx = complex<double>(Amp * cos(Phase), Amp * sin(Phase));

					// Hy
					f1 = pOptical_plane[nn].Hy_save_array_t1[ii][jj];
					f2 = pOptical_plane[nn].Hy_save_array_t2[ii][jj];
					Phase = atan2((f2 * sin(2 * pi * f0 * dt * tt1) - f1 * sin(2 * pi * f0 * dt * tt2)),
						(f1 * cos(2 * pi * f0 * dt * tt2) - f2 * cos(2 * pi * f0 * dt * tt1)));
					Amp = f1 / sin(2 * pi * f0 * dt * tt1 + Phase);
					Amp_Phase_Tunning();
					optical_Hy = complex<double>(Amp * cos(Phase), Amp * sin(Phase));

					// Hz
					f1 = pOptical_plane[nn].Hz_save_array_t1[ii][jj];
					f2 = pOptical_plane[nn].Hz_save_array_t2[ii][jj];
					Phase = atan2((f2 * sin(2 * pi * f0 * dt * tt1) - f1 * sin(2 * pi * f0 * dt * tt2)),
						(f1 * cos(2 * pi * f0 * dt * tt2) - f2 * cos(2 * pi * f0 * dt * tt1)));
					Amp = f1 / sin(2 * pi * f0 * dt * tt1 + Phase);
					Amp_Phase_Tunning();
					optical_Hz = complex<double>(Amp * cos(Phase), Amp * sin(Phase));

					temp = (optical_Ex * conj(optical_Hx) + optical_Ey * conj(optical_Hy) + optical_Ez * conj(optical_Hz));
					pOptical_plane[nn].pOptical_chirality[ii][jj] = -(eps0 * mu0) * (2 * pi * f0) / 2 * imag(temp);
				}
		}
		save_optical_chirality_Frequency();
	}

	//if (Nststep_iter >= steady_state_start_discrete && Nststep_iter < steady_state_start_discrete + 2 * T0_discrete)
		//if (Nststep_iter >= 1 && Nststep_iter < steady_state_start_discrete + 2 * T0_discrete)
	// if (0)
	// {
	// 	store_2D_field_frequency(Nststep_iter);
	// 	steady_state_offset = Nststep_iter - save_sta;
	// 	if (steady_state_offset % T0_stepSize == 0)
	// 	{
	// 		optical_chirality_caclulate();
	// 		for (int LL = 0; LL < Layer; LL++)
	// 		{
	// 			sprintf(Fname_CC, "./result/optical_chirality/Optical_chirality_[%d][%d].dat", LL, Nststep_iter);
	// 			pF_CC = fopen(Fname_CC, "wb");
	// 			for (int ii = 0; ii < txsize; ii++)
	// 			{
	// 				fwrite(pOptical_chirality[LL][ii], tysize, sizeof(double), pF_CC);
	// 			}
	// 			fclose(pF_CC);
	// 		}
	// 	}
	// }
}

void Toptical_chirality::Prepare2save()
{
	double steady_points = (0.01 / pFc->f0 / pFc->dt) / (double)(pFc->NSTEPS);
	set_steady_start(0.85);
	set_steady_end(0.85 + steady_points);
	Initiate_optical_option_Frequency();
	Exctract_chiral_surface();
	MemoryAssign_for_calculate_optical_chirality();
}

void Toptical_chirality::Initiate_optical_option_Frequency()
{
	ifstream ifs(JsonFileName);
	IStreamWrapper isw(ifs);
	Jdocument.ParseStream(isw);

	this->f0 = Jdocument["SIMULATION_ENVIORONMENT"]["SOURCE"]["f0"].GetDouble();
	this->f0 = f0;
	this->sfre = f0;
	this->tfre = f0;
	this->sampling = 1;
	this->df = round((tfre - sfre) / (double)sampling);
	// this->source_pos = Jdocument
	ifs.close();
}

void Toptical_chirality::Exctract_chiral_surface()
{
	// int Air = 0;
	// int PDMS = 2;
	// int Au = 3;
	// int LL = 0;
	// int Water = 1;
	// Get_Pulse_infomation();
	// Layer = no_composit - 2 + 3;
	// 굴곡 경계 추출 알고리즘.
	
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
		return;
	}
	
	Layer = 2 * ActualNum + 3; // pivot_min-1, max+1, middle
	// up, down, middle
	pivot_min = pow(2, 20);
	pivot_max = 0;
	// Num. Layer check;
	pOptical_chirality			 = Tptr3<double>(ActualNum, txsize, tysize);
	pOptical_chirality_frequency = Tptr3<double>(ActualNum, txsize, tysize);
	// Plane_save_name					 = Tptr2<char>(Layer, 256);
	int nn = 0;
	for (int nn1 = 0; nn1 < ActualNum; nn1++)
	{
		for (int nn2 = 0; nn2 < 2; nn2++)
		{
//			nn = nn1 * 2 + nn2;
			for (int ii = 0; ii < txsize; ii++)
			for (int jj = 0; jj < tysize; jj++)
			{
				if (psEnvelop_Material[nn1].pEnvelop[nn2][ii][jj] == 0)
					continue;
				else if (pivot_min > psEnvelop_Material[nn1].pEnvelop[nn2][ii][jj])
					pivot_min = psEnvelop_Material[nn1].pEnvelop[nn2][ii][jj];

				if (psEnvelop_Material[nn1].pEnvelop[nn2][ii][jj] == 0)
					continue;
				else if (pivot_max < psEnvelop_Material[nn1].pEnvelop[nn2][ii][jj])
					pivot_max = psEnvelop_Material[nn1].pEnvelop[nn2][ii][jj];
			}
			// UD기준 middle 위치 필요.
		}
	}

	// optical chirality를 계산할 위치 : pivot - Au length
	// 1. pivot - 3;
	// 2. pSOC[3]
	char LayerName_up[256];
	char LayerName_down[256];
	FILE* p2D_optical_chirality;
	for (int nn1 = 0; nn1 < ActualNum; nn1++)
	for (int nn2 = 0; nn2 < 2; nn2++)
	{
		sprintf(LayerName_up, "./result/2D_chiral[%s]_up.dat", psEnvelop_Material[nn1].Name);
		sprintf(LayerName_down, "./result/2D_chiral[%s]_down.dat", psEnvelop_Material[nn1].Name);
		
		p2D_optical_chirality= fopen(LayerName_up, "wb");
		for (int ii = 0; ii < txsize; ii++)
			fwrite(psEnvelop_Material[nn1].pEnvelop[nn2][ii], sizeof(int), tysize, p2D_optical_chirality);
		fclose(p2D_optical_chirality);

		p2D_optical_chirality = fopen(LayerName_down, "wb");
		for (int ii = 0; ii < txsize; ii++)
			fwrite(psEnvelop_Material[nn1].pEnvelop[nn2][ii], sizeof(int), tysize, p2D_optical_chirality);
		fclose(p2D_optical_chirality);
	}
	
}

void Toptical_chirality::Single_Freq_Amp_Phase1(int xx, int yy)
{
	int zz = 0, ff = 0, nn = 0;
	//#pragma omp parallel private(zz, nn, ff)
	for (int nn1 = 0; nn1 < ActualNum; nn1++)
	{
		for (int nn2 = 0; nn2 < 2; nn2++)
		{
			nn = nn1*2 + nn2;
			zz = psEnvelop_Material[nn1].pEnvelop[nn2][xx][yy];
			
			Voxel_center_average(xx, yy, zz, "optical_C_p1");
			pOptical_plane[nn].Ex_save_array_t1[xx][yy] = ExV;
			pOptical_plane[nn].Ey_save_array_t1[xx][yy] = EyV;
			pOptical_plane[nn].Ez_save_array_t1[xx][yy] = EzV;
											 
			pOptical_plane[nn].Hx_save_array_t1[xx][yy] = HxV;
			pOptical_plane[nn].Hy_save_array_t1[xx][yy] = HyV;
			pOptical_plane[nn].Hz_save_array_t1[xx][yy] = HzV;
		}
	}
	// Layer 나머지, min, max, middle
	int zz1[3];
	zz1[0] = pivot_min;
	zz1[1] = pivot_max;
	zz1[2] = (pivot_min + pivot_max) / 2;
	for (int nn1 = nn+1; nn1 < Layer; nn1++)
	{
		int	nn2 = nn1 - (nn+1);
		zz = zz1[nn2];
		Voxel_center_average(xx, yy, zz, "Layer Position P1");
		pOptical_plane[nn1].Ex_save_array_t1[xx][yy] = ExV;
		pOptical_plane[nn1].Ey_save_array_t1[xx][yy] = EyV;
		pOptical_plane[nn1].Ez_save_array_t1[xx][yy] = EzV;

		pOptical_plane[nn1].Hx_save_array_t1[xx][yy] = HxV;
		pOptical_plane[nn1].Hy_save_array_t1[xx][yy] = HyV;
		pOptical_plane[nn1].Hz_save_array_t1[xx][yy] = HzV;
		
	}
}

void Toptical_chirality::Single_Freq_Amp_Phase2(int xx, int yy)
{
	int zz = 0, ff = 0, nn = 0;
	//#pragma omp parallel private(zz, nn, ff)
	for (int nn1 = 0; nn1 < ActualNum; nn1++)
	{
		for (int nn2 = 0; nn2 < 2; nn2++)
		{
			nn = nn1 * 2 + nn2;
			zz = psEnvelop_Material[nn1].pEnvelop[nn2][xx][yy];
			Voxel_center_average(xx, yy, zz, "optical_C_p2");
			pOptical_plane[nn].Ex_save_array_t2[xx][yy] = ExV;
			pOptical_plane[nn].Ey_save_array_t2[xx][yy] = EyV;
			pOptical_plane[nn].Ez_save_array_t2[xx][yy] = EzV;
											  
			pOptical_plane[nn].Hx_save_array_t2[xx][yy] = HxV;
			pOptical_plane[nn].Hy_save_array_t2[xx][yy] = HyV;
			pOptical_plane[nn].Hz_save_array_t2[xx][yy] = HzV;
		}
	}
	// Layer 나머지, min, max, middle
	// Layer 나머지, min, max, middle
	int zz1[3];
	zz1[0] = pivot_min;
	zz1[1] = pivot_max;
	zz1[2] = (pivot_min + pivot_max) / 2;
	for (int nn1 = nn + 1; nn1 < Layer; nn1++)
	{
		int	nn2 = nn1 - (nn + 1);
		zz = zz1[nn2];
		Voxel_center_average(xx, yy, zz, "Layer Position P2");
		pOptical_plane[nn1].Ex_save_array_t2[xx][yy] = ExV;
		pOptical_plane[nn1].Ey_save_array_t2[xx][yy] = EyV;
		pOptical_plane[nn1].Ez_save_array_t2[xx][yy] = EzV;
											   
		pOptical_plane[nn1].Hx_save_array_t2[xx][yy] = HxV;
		pOptical_plane[nn1].Hy_save_array_t2[xx][yy] = HyV;
		pOptical_plane[nn1].Hz_save_array_t2[xx][yy] = HzV;
		
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Frequency domain
// void Toptical_chirality::store_2D_field_frequency(int Nststep_iter)
// {
// 	//#pragma omp parallel
// 	//{
// 	//#pragma omp for schedule(dynamic)
// 	for (int ii = 0; ii < txsize; ii++)
// 	for (int jj = 0; jj < tysize; jj++)
// 	{ 
// 		DFT_chiral_surface_wide(ii, jj, Nststep_iter);
// 	}
// }

void Toptical_chirality::MemoryAssign_for_calculate_optical_chirality()
{
	// 5개의 면,
	// 1. 카이랄 평면 2등분해서 처음 가운데 끝
	// 2. 투과면 특정지점
	// 3. 카이랄 경계면
	pOptical_plane = Tptr1<struct_optical_chirality>(Layer);
	for(int nn = 0 ; nn < Layer; nn++)
	{
		pOptical_plane[nn].pOptical_chirality = Tptr2<double>(txsize, tysize);
		pOptical_plane[nn].Ex_save_array_t1 = Tptr2<double>(txsize, tysize); pOptical_plane[nn].Ex_save_array_t2 = Tptr2<double>(txsize, tysize);
		pOptical_plane[nn].Ey_save_array_t1 = Tptr2<double>(txsize, tysize); pOptical_plane[nn].Ey_save_array_t2 = Tptr2<double>(txsize, tysize);
		pOptical_plane[nn].Ez_save_array_t1 = Tptr2<double>(txsize, tysize); pOptical_plane[nn].Ez_save_array_t2 = Tptr2<double>(txsize, tysize);
		pOptical_plane[nn].Hx_save_array_t1 = Tptr2<double>(txsize, tysize); pOptical_plane[nn].Hx_save_array_t2 = Tptr2<double>(txsize, tysize);
		pOptical_plane[nn].Hy_save_array_t1 = Tptr2<double>(txsize, tysize); pOptical_plane[nn].Hy_save_array_t2 = Tptr2<double>(txsize, tysize);
		pOptical_plane[nn].Hz_save_array_t1 = Tptr2<double>(txsize, tysize); pOptical_plane[nn].Hz_save_array_t2 = Tptr2<double>(txsize, tysize);
	}
}

//void Toptical_chirality::DFT_chiral_surface_wide(int xx, int yy, int Nststep_iter)
//{
//	for (int nn1 = 0; nn1 < ActualNum; nn1++)
//	for (int nn2 = 0; nn2 < 2; nn2++)
//	{
//
//		int zz = pExtract_Optical_Chirality_plane[nn][xx][yy];
//		Voxel_center_average(xx, yy, zz);
//
//		pFourier->onthefly_dft(ExV, Nststep_iter);
//		pSave_Temp_re = pFourier->Get_Frequency_Domain_Response_Real();
//		pSave_Temp_im = pFourier->Get_Frequency_Domain_Response_Imag();
//
//		for (int ff = 0; ff < sampling; ff++)
//		{
//			pOptical_plane[nn][xx][yy].Ex_save_array_re[ff] = pOptical_plane[nn][xx][yy].Ex_save_array_re[ff] + pSave_Temp_re[ff];
//			pOptical_plane[nn][xx][yy].Ex_save_array_im[ff] = pOptical_plane[nn][xx][yy].Ex_save_array_im[ff] + pSave_Temp_im[ff];
//		}
//		pFourier->onthefly_dft(EyV, Nststep_iter);
//		pSave_Temp_re = pFourier->Get_Frequency_Domain_Response_Real();
//		pSave_Temp_im = pFourier->Get_Frequency_Domain_Response_Imag();
//
//		for (int ff = 0; ff < sampling; ff++)
//		{
//			pOptical_plane[nn][xx][yy].Ey_save_array_re[ff] = pOptical_plane[nn][xx][yy].Ey_save_array_re[ff] + pSave_Temp_re[ff];
//			pOptical_plane[nn][xx][yy].Ey_save_array_im[ff] = pOptical_plane[nn][xx][yy].Ey_save_array_im[ff] + pSave_Temp_im[ff];
//		}
//
//		pFourier->onthefly_dft(EzV, Nststep_iter);
//		pSave_Temp_re = pFourier->Get_Frequency_Domain_Response_Real();
//		pSave_Temp_im = pFourier->Get_Frequency_Domain_Response_Imag();
//
//		for (int ff = 0; ff < sampling; ff++)
//		{
//			pOptical_plane[nn][xx][yy].Ez_save_array_re[ff] = pOptical_plane[nn][xx][yy].Ez_save_array_re[ff] + pSave_Temp_re[ff];
//			pOptical_plane[nn][xx][yy].Ez_save_array_im[ff] = pOptical_plane[nn][xx][yy].Ez_save_array_im[ff] + pSave_Temp_im[ff];
//		}
//
//		//
//		pFourier->onthefly_dft(HxV, Nststep_iter);
//		pSave_Temp_re = pFourier->Get_Frequency_Domain_Response_Real();
//		pSave_Temp_im = pFourier->Get_Frequency_Domain_Response_Imag();
//
//		for (int ff = 0; ff < sampling; ff++)
//		{
//			pOptical_plane[nn][xx][yy].Hx_save_array_re[ff] = pOptical_plane[nn][xx][yy].Hx_save_array_re[ff] + pSave_Temp_re[ff];
//			pOptical_plane[nn][xx][yy].Hx_save_array_im[ff] = pOptical_plane[nn][xx][yy].Hx_save_array_im[ff] + pSave_Temp_im[ff];
//		}
//
//		pFourier->onthefly_dft(HyV, Nststep_iter);
//		pSave_Temp_re = pFourier->Get_Frequency_Domain_Response_Real();
//		pSave_Temp_im = pFourier->Get_Frequency_Domain_Response_Imag();
//
//		for (int ff = 0; ff < sampling; ff++)
//		{
//			pOptical_plane[nn][xx][yy].Hy_save_array_re[ff] = pOptical_plane[nn][xx][yy].Hy_save_array_re[ff] + pSave_Temp_re[ff];
//			pOptical_plane[nn][xx][yy].Hy_save_array_im[ff] = pOptical_plane[nn][xx][yy].Hy_save_array_im[ff] + pSave_Temp_im[ff];
//		}
//
//		pFourier->onthefly_dft(HzV, Nststep_iter);
//		pSave_Temp_re = pFourier->Get_Frequency_Domain_Response_Real();
//		pSave_Temp_im = pFourier->Get_Frequency_Domain_Response_Imag();
//
//		for (int ff = 0; ff < sampling; ff++)
//		{
//			pOptical_plane[nn][xx][yy].Hz_save_array_re[ff] = pOptical_plane[nn][xx][yy].Hz_save_array_re[ff] + pSave_Temp_re[ff];
//			pOptical_plane[nn][xx][yy].Hz_save_array_im[ff] = pOptical_plane[nn][xx][yy].Hz_save_array_im[ff] + pSave_Temp_im[ff];
//		}
//	}
//}

//void Toptical_chirality::Calculate_optical_chirality_Frequency()
//{
//	int nn = 0;
//	for (int nn1 = 0; nn1 < ActualNum; nn1++)
//	for (int nn2 = 0; nn2 < 2; nn2++)
//	{
//		nn = nn1 * nn2 + nn2;
//		// optical chirliaty 를 구하는 평면
//		for (int ii = 0; ii < txsize; ii++)
//		for (int jj = 0; jj < tysize; jj++)
//		{
//				pOptical_plane[nn].pOptical_chirality[ii][jj] = -(eps0*mu0)*(2 * pi*f0) / 2 *
//					   (pOptical_plane[nn].Ex_save_array_t1[ii][jj] * (-1) * pOptical_plane[nn].Hx_save_array_t1[ii][jj] +
//						pOptical_plane[nn].Ex_save_array_t2[ii][jj] * (+1) * pOptical_plane[nn].Hx_save_array_t1[ii][jj] +
//						pOptical_plane[nn].Ey_save_array_t1[ii][jj] * (-1) * pOptical_plane[nn].Hy_save_array_t2[ii][jj] +
//						pOptical_plane[nn].Ey_save_array_t2[ii][jj] * (+1) * pOptical_plane[nn].Hy_save_array_t1[ii][jj] +
//						pOptical_plane[nn].Ez_save_array_t1[ii][jj] * (-1) * pOptical_plane[nn].Hz_save_array_t2[ii][jj] +
//						pOptical_plane[nn].Ez_save_array_t2[ii][jj] * (+1) * pOptical_plane[nn].Hz_save_array_t1[ii][jj]);
//
//		}
//	}
//}

// time marching 이후 저장. (sine wave)
void Toptical_chirality::save_optical_chirality_Frequency()
{
	int nn = 0;
	for (int nn1 = 0; nn1 < ActualNum; nn1++)
	for (int nn2 = 0; nn2 < 2; nn2++)
	{
		nn = nn1 * 2 + nn2;
		if (nn2 == 0)
			sprintf(Fname_optical_chirality, "./result/optical_chirality/frequency_[%s]_up", psEnvelop_Material[nn1].Name);
		else if (nn2 == 1)
			sprintf(Fname_optical_chirality, "./result/optical_chirality/frequency_[%s]_down", psEnvelop_Material[nn1].Name);
		_mkdir(Fname_optical_chirality);

		if (nn2 == 0)
			sprintf(Fname_optical_chirality, "./result/optical_chirality/frequency_[%s]_up/Optical_chirality_[%1.2le].dat", psEnvelop_Material[nn1].Name, f0);
		else if (nn2 == 1)
			sprintf(Fname_optical_chirality, "./result/optical_chirality/frequency_[%s]_down/Optical_chirality_[%1.2le].dat", psEnvelop_Material[nn1].Name, f0);
		pFname_optical_chirality = fopen(Fname_optical_chirality, "wb");

		for (int ii = 0; ii < txsize; ii++)
		{
			fwrite(pOptical_plane[nn].pOptical_chirality[ii], sizeof(double), tysize, pFname_optical_chirality);
		}
		fclose(pFname_optical_chirality);
	}

	sprintf(Fname_optical_chirality, "./result/optical_chirality/frequency_[min]");
	_mkdir(Fname_optical_chirality);
	sprintf(Fname_optical_chirality, "./result/optical_chirality/frequency_[min]/Optical_chirality_[%1.2le].dat", f0);
	pFname_optical_chirality = fopen(Fname_optical_chirality, "wb");
	for (int ii = 0; ii < txsize; ii++)
	{
		fwrite(pOptical_plane[nn].pOptical_chirality[ii], sizeof(double), tysize, pFname_optical_chirality);
	}
	fclose(pFname_optical_chirality);

	sprintf(Fname_optical_chirality, "./result/optical_chirality/frequency_[max]");
	_mkdir(Fname_optical_chirality);
	sprintf(Fname_optical_chirality, "./result/optical_chirality/frequency_[max]/Optical_chirality_[%1.2le].dat", f0);
	pFname_optical_chirality = fopen(Fname_optical_chirality, "wb");
	for (int ii = 0; ii < txsize; ii++)
	{
		fwrite(pOptical_plane[nn+1].pOptical_chirality[ii], sizeof(double), tysize, pFname_optical_chirality);
	}
	fclose(pFname_optical_chirality);

	sprintf(Fname_optical_chirality, "./result/optical_chirality/frequency_[mid]/");
	_mkdir(Fname_optical_chirality);
	sprintf(Fname_optical_chirality, "./result/optical_chirality/frequency_[mid]/Optical_chirality_[%1.2le].dat", f0);
	pFname_optical_chirality = fopen(Fname_optical_chirality, "wb");
	for (int ii = 0; ii < txsize; ii++)
	{
		fwrite(pOptical_plane[nn+2].pOptical_chirality[ii], sizeof(double), tysize, pFname_optical_chirality);
	}
	fclose(pFname_optical_chirality);

}

void Toptical_chirality::optical_chirality_caclulate()
{
	double optical_chirality_E, optical_chirality_H;
	int ii, jj, LL, kk;
	ii = 0; jj = 0; kk = 0; LL = 0;
#pragma omp parallel private(ii, jj, kk, LL)
	{
		for (LL = 0; LL < Layer; LL++)
		{
			for (ii = 1; ii < txsize - 1; ii++)
			{
#pragma omp for schedule(dynamic)
				for (jj = 1; jj < tysize - 1; jj++)
				{
					kk = pExtract_Optical_Chirality_plane[LL][ii][jj];
					optical_chirality_E =
						E_dot_p(ii, jj, kk, 'x') * curl_E(ii, jj, kk, 'x')
						+ E_dot_p(ii, jj, kk, 'y') * curl_E(ii, jj, kk, 'y')
						+ E_dot_p(ii, jj, kk, 'z') * curl_E(ii, jj, kk, 'z');

					optical_chirality_H =
						H_dot_p(ii, jj, kk, 'x') * curl_H(ii, jj, kk, 'x')
						+ H_dot_p(ii, jj, kk, 'y') * curl_H(ii, jj, kk, 'y')
						+ H_dot_p(ii, jj, kk, 'z') * curl_H(ii, jj, kk, 'z');
					pOptical_chirality[LL][ii][jj] = eps0 / 2.0 * optical_chirality_E + mu0 / 2.0 * optical_chirality_H;
				}
			}
		}
	}
}

double Toptical_chirality::curl_E(int ii, int jj, int kk, char comp)
{
	double E_output1, E_output2, E_output3, E_output4;
	double E_output;
	if ((toupper(comp)) == 'X')
	{
		E_output1 = EOz[ii + 1][jj + 1][kk] - EOz[ii + 1][jj][kk];
		E_output2 = EOz[ii][jj + 1][kk] - EOz[ii][jj][kk];
		E_output3 = EOy[ii + 1][jj][kk + 1] - EOy[ii + 1][jj][kk];
		E_output4 = EOy[ii][jj][kk + 1] - EOy[ii][jj][kk];
		E_output = (E_output1 + E_output2) / 2.0 / dy[jj] - (E_output3 + E_output4) / 2.0 / dz[kk];
	}
	if ((toupper(comp)) == 'Y')
	{
		E_output1 = EOx[ii][jj + 1][kk + 1] - EOx[ii][jj + 1][kk];
		E_output2 = EOx[ii][jj][kk + 1] - EOx[ii][jj][kk];
		E_output3 = EOz[ii + 1][jj + 1][kk] - EOz[ii][jj + 1][kk];
		E_output4 = EOz[ii + 1][jj][kk] - EOz[ii][jj][kk];
		E_output = (E_output1 + E_output2) / 2.0 / dz[kk] - (E_output3 + E_output4) / 2.0 / dx[ii];
	}
	if ((toupper(comp)) == 'Z')
	{
		E_output1 = EOy[ii + 1][jj][kk + 1] - EOy[ii][jj][kk + 1];
		E_output2 = EOy[ii + 1][jj][kk] - EOy[ii][jj][kk];
		E_output3 = EOx[ii][jj + 1][kk + 1] - EOx[ii][jj][kk + 1];
		E_output4 = EOx[ii][jj + 1][kk] - EOx[ii][jj][kk];
		E_output = (E_output1 + E_output2) / 2.0 / dx[ii] - (E_output3 + E_output4) / 2.0 / dy[jj];
	}
	return E_output;
}

double Toptical_chirality::E_dot_p(int ii, int jj, int kk, char comp)
{
	double E_output = 0;;
	if ((toupper(comp)) == 'X')
	{
		E_output = EOx[ii][jj][kk] + EOx[ii][jj + 1][kk] + EOx[ii][jj][kk + 1] + EOx[ii][jj + 1][kk + 1];
	}
	else if ((toupper(comp)) == 'Y')
	{
		E_output = EOy[ii][jj][kk] + EOy[ii + 1][jj][kk] + EOy[ii][jj][kk + 1] + EOy[ii + 1][jj][kk + 1];
	}
	else if ((toupper(comp)) == 'Z')
	{
		E_output = EOz[ii][jj][kk] + EOz[ii + 1][jj][kk] + EOz[ii][jj + 1][kk] + EOz[ii + 1][jj + 1][kk];
	}
	return E_output / 4.0;
}

double Toptical_chirality::curl_H(int ii, int jj, int kk, char comp)
{
	double H_output1, H_output2, H_output3, H_output4;
	double H_output;
	if ((toupper(comp)) == 'X')
	{
		H_output1 = (HOz[ii][jj + 1][kk] - HOz[ii][jj - 1][kk]) / 4.0 / dy[jj];
		H_output2 = (HOz[ii][jj + 1][kk + 1] - HOz[ii][jj - 1][kk + 1]) / 4.0 / dy[jj];
		H_output3 = (HOy[ii][jj][kk + 1] - HOy[ii][jj][kk - 1]) / 4.0 / dz[kk];
		H_output4 = (HOy[ii][jj + 1][kk + 1] - HOy[ii][jj + 1][kk - 1]) / 4.0 / dz[kk];
		H_output = (H_output1 + H_output2) - (H_output3 + H_output4);
	}
	if ((toupper(comp)) == 'Y')
	{
		H_output1 = (HOx[ii][jj][kk + 1] - HOx[ii][jj][kk - 1]) / 4.0 / dz[kk];
		H_output2 = (HOx[ii + 1][jj][kk + 1] - HOx[ii + 1][jj][kk - 1]) / 4.0 / dz[kk];
		H_output3 = (HOz[ii + 1][jj][kk] - HOz[ii - 1][jj + 1][kk]) / 4.0 / dx[ii];
		H_output4 = (HOz[ii + 1][jj][kk + 1] - HOz[ii - 1][jj - 1][kk + 1]) / 4.0 / dx[ii];
		H_output = (H_output1 + H_output2) - (H_output3 + H_output4);
	}
	if ((toupper(comp)) == 'Z')
	{
		H_output1 = (HOy[ii + 1][jj][kk] - HOy[ii - 1][jj][kk]) / 4.0 / dx[ii];
		H_output2 = (HOy[ii + 1][jj + 1][kk] - HOy[ii - 1][jj + 1][kk]) / 4.0 / dx[ii];
		H_output3 = (HOx[ii][jj + 1][kk] - HOx[ii][jj - 1][kk]) / 4.0 / dy[jj];
		H_output4 = (HOx[ii + 1][jj + 1][kk] - HOx[ii + 1][jj - 1][kk]) / 4.0 / dy[jj];
		H_output = (H_output1 + H_output2) - (H_output3 + H_output4);
	}
	return H_output;
}

double Toptical_chirality::H_dot_p(int ii, int jj, int kk, char comp)
{
	double H_output = 0;
	if ((toupper(comp)) == 'X')
	{
		H_output = HOx[ii][jj][kk] + HOx[ii + 1][jj][kk];
	}
	else if ((toupper(comp)) == 'Y')
	{
		H_output = HOy[ii][jj][kk] + HOy[ii][jj + 1][kk];
	}
	else if ((toupper(comp)) == 'Z')
	{
		H_output = HOz[ii][jj][kk] + HOz[ii][jj][kk + 1];
	}
	return H_output / 2.0;
}

void Toptical_chirality::Amp_Phase_Tunning()
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