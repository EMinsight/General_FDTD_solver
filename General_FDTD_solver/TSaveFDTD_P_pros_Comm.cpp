#include "TSaveFDTD_P_pros_Comm.h"

void TSaveFDTD_P_pros_Comm::set_TimeDiff_NstepMax(double dt, int Nststep_max)
{
	pi = acos(-1);
	eps_0 = 1e-9 / 36.0 / pi;
	this->dt = dt;
	if (this->Nststep_Max != Nststep_max)
		this->Nststep_Max = Nststep_max;
	if (Nststep_max > 100)
		Nststep_save = Nststep_max * 0.01;
	else
		Nststep_save = Nststep_max;
}

void TSaveFDTD_P_pros_Comm::set_postprocessing_info(struct_Field_carrier* pFc, Material_properties* pMp)
{
	this->pFc = pFc;
	RotationZXZp = Identity_matrix_2D<double>(3, 3);
	this->EOx = pFc->Ex_address;
	this->EOy = pFc->Ey_address;
	this->EOz = pFc->Ez_address;
	this->HOx = pFc->Hx_address;
	this->HOy = pFc->Hy_address;
	this->HOz = pFc->Hz_address;
	this->dx = pFc->dx;
	this->dy = pFc->dy;
	this->dz = pFc->dz;
	this->dt = pFc->dt;
	this->Media = pFc->media;
	this->txsize = pFc->txsize;
	this->tysize = pFc->tysize;
	this->tzsize = pFc->tzsize;
	this->pMp = pMp;
	this->BW = pFc->BW;
	this->f0 = pFc->f0;
	

	strcpy(JsonFileName, "./IEM_solver_input_format.json");
	set_TimeDiff_NstepMax(dt, pFc->NSTEPS);
	
}

void TSaveFDTD_P_pros_Comm::set_steady_start(double percent)
{
	this->save_sta = percent * Nststep_Max;
}
void TSaveFDTD_P_pros_Comm::set_steady_end(double percent)
{
	this->save_end = percent * Nststep_Max;
}

void TSaveFDTD_P_pros_Comm::set_T0(double T0)
{
	this->T0 = T0;
	T0_discrete = round(T0 / dt);
	if (T0_discrete > 100)
		T0_stepSize = 1; // round(T0_discrete / 100.0);
	else
		T0_stepSize = 1;
}

void TSaveFDTD_P_pros_Comm::Voxel_face_average(int xx, int yy, int zz, string who)
{
	//. PBC 경계는 나중에 고려....
	// E -fields
	if (xx + 1 < txsize && yy + 1 < tysize && zz + 1 < tzsize)
	{
		ExV = (EOx[xx][yy][zz] + EOx[xx][yy + 1][zz] + EOx[xx][yy][zz + 1] + EOx[xx][yy + 1][zz + 1]+
			  EOx[xx + 1][yy][zz] + EOx[xx + 1][yy + 1][zz] + EOx[xx + 1][yy][zz + 1] + EOx[xx + 1][yy + 1][zz + 1]) / 8.0;
		EyV = (EOy[xx][yy][zz] + EOy[xx + 1][yy][zz] + EOy[xx][yy][zz + 1] + EOy[xx + 1][yy][zz + 1] +
			   EOy[xx][yy + 1][zz] + EOy[xx + 1][yy + 1][zz] + EOy[xx][yy + 1][zz + 1] + EOy[xx + 1][yy + 1][zz + 1]) / 8.0;
		EzV = (EOz[xx][yy][zz] + EOz[xx + 1][yy][zz] + EOz[xx][yy + 1][zz] + EOz[xx + 1][yy + 1][zz] + 
			   EOz[xx][yy][zz + 1] + EOz[xx + 1][yy][zz + 1] + EOz[xx][yy + 1][zz + 1] + EOz[xx + 1][yy + 1][zz + 1]) / 8.0;
	}

	// J -fields
	if ( (xx + 1) < txsize && (yy + 1) < tysize && (zz + 1) < tzsize)
	{
		JxV = 0;		JyV = 0;		JzV = 0;
		for (int pp = 0 ; pp < pFc->total_pol ; pp++)
		{
			JxV += (pFc->Jx_address[xx    ][yy    ][zz    ][pp] + 
				    pFc->Jx_address[xx    ][yy + 1][zz    ][pp] + 
				    pFc->Jx_address[xx    ][yy    ][zz + 1][pp] + 
				    pFc->Jx_address[xx    ][yy + 1][zz + 1][pp] +
				    pFc->Jx_address[xx + 1][yy    ][zz    ][pp] + 
				    pFc->Jx_address[xx + 1][yy + 1][zz    ][pp] + 
				    pFc->Jx_address[xx + 1][yy    ][zz + 1][pp] + 
				    pFc->Jx_address[xx + 1][yy + 1][zz + 1][pp]) / 8.0;
			
			JyV += (pFc->Jy_address[xx    ][yy    ][zz    ][pp]+ 
				    pFc->Jy_address[xx + 1][yy    ][zz    ][pp]+ 
				    pFc->Jy_address[xx    ][yy    ][zz + 1][pp]+ 
				    pFc->Jy_address[xx + 1][yy    ][zz + 1][pp]+
				    pFc->Jy_address[xx    ][yy + 1][zz    ][pp]+ 
				    pFc->Jy_address[xx + 1][yy + 1][zz    ][pp]+ 
				    pFc->Jy_address[xx    ][yy + 1][zz + 1][pp]+ 
				    pFc->Jy_address[xx + 1][yy + 1][zz + 1][pp]) / 8.0;

			JzV += (pFc->Jz_address[xx    ][yy    ][zz    ][pp]+ 
				    pFc->Jz_address[xx + 1][yy    ][zz    ][pp]+ 
				    pFc->Jz_address[xx    ][yy + 1][zz    ][pp]+ 
				    pFc->Jz_address[xx + 1][yy + 1][zz    ][pp]+
				    pFc->Jz_address[xx    ][yy    ][zz + 1][pp]+ 
				    pFc->Jz_address[xx + 1][yy    ][zz + 1][pp]+ 
				    pFc->Jz_address[xx    ][yy + 1][zz + 1][pp]+ 
				    pFc->Jz_address[xx + 1][yy + 1][zz + 1][pp]) / 8.0;
		}
	}

}

void TSaveFDTD_P_pros_Comm::Voxel_center_average(int xx, int yy, int zz, string who)
{
	ExV = (EOx[xx][yy][zz] +
		EOx[xx][yy + 1][zz] +
		EOx[xx][yy][zz + 1] +
		EOx[xx][yy + 1][zz + 1]) / 4.0;
	EyV = (EOy[xx][yy][zz] +
		EOy[xx + 1][yy][zz] +
		EOy[xx][yy][zz + 1] +
		EOy[xx + 1][yy][zz + 1]) / 4.0;
	EzV = (EOz[xx][yy][zz] +
		EOz[xx + 1][yy][zz] +
		EOz[xx][yy + 1][zz] +
		EOz[xx + 1][yy + 1][zz]) / 4.0;

	HxV = (HOx[xx][yy][zz] + HOx[xx + 1][yy][zz]) / 2.0;
	HyV = (HOy[xx][yy][zz] + HOy[xx][yy + 1][zz]) / 2.0;
	HzV = (HOz[xx][yy][zz] + HOz[xx][yy][zz + 1]) / 2.0;
}
void TSaveFDTD_P_pros_Comm::Amp_Phase_Tunning()
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
void TSaveFDTD_P_pros_Comm::makeRotation(double fhi, double theta, double psi)
{
	memfree(&RotationZXZp, 3);
	double** RotationZ  = Identity_matrix_2D<double>(3, 3);
	double** RotationX  = Identity_matrix_2D<double>(3, 3);
	double** RotationZp = Identity_matrix_2D<double>(3, 3);
	fhi = fhi * acos(-1) / 180;
	theta = theta * acos(-1) / 180;
	psi = psi * acos(-1) / 180;

	RotationZ[0][0] = cos(fhi);	    RotationZ[0][1] = sin(fhi); RotationZ[0][2] = 0;
	RotationZ[1][0] = -sin(fhi); 	RotationZ[1][1] = cos(fhi); RotationZ[1][2] = 0;
	RotationZ[2][0] = 0;		    RotationZ[2][1] = 0;		RotationZ[2][2] = 1;

	cout << "Rotation Z - axis "<< fhi << " deg" << endl;
	for (int ii = 0; ii < 3; ii++)
	{
		for (int jj = 0; jj < 3; jj++)
			cout << RotationZ[ii][jj] << "\t";
		cout << endl;
	}

	RotationX[0][0] = 1;			RotationX[0][1] = 0; RotationX[0][2] = 0;
	RotationX[1][1] = cos(theta);	RotationX[1][1] = 0; RotationX[1][2] = sin(theta);
	RotationX[2][1] = -sin(theta);	RotationX[2][1] = 0; RotationX[2][2] = cos(theta);

	cout << "Rotation X - axis " << theta << " deg" << endl;
	for (int ii = 0; ii < 3; ii++)
	{
		for (int jj = 0; jj < 3; jj++)
			cout << RotationX[ii][jj] << "\t";
		cout << endl;
	}

	RotationZp[0][0] = cos(psi);	RotationZp[0][1] = sin(psi); RotationZp[0][2] = 0;
	RotationZp[1][0] = -sin(psi);	RotationZp[1][1] = cos(psi); RotationZp[1][2] = 0;
	RotationZp[2][2] = 1;			RotationZp[2][1] = 0;		 RotationZp[2][2] = 0;
	cout << "Rotation X - axis " << psi << " deg" << endl;
	cout << "Rotation Zp" << endl;
	for (int ii = 0; ii < 3; ii++)
	{
		for (int jj = 0; jj < 3; jj++)
			cout << RotationZp[ii][jj] << "\t";
		cout << endl;
	}

	cout << "RotationZX" << endl;
	double** RotationZX = Identity_matrix_2D<double>(3, 3); // R4
	for (int ii = 0; ii < 3; ii++)
	{
		for (int jj = 0; jj < 3; jj++)
		{
			double temp = 0;
			for (int ll = 0; ll < 3; ll++)
				temp += RotationZ[ii][ll] * RotationX[ll][jj];
			RotationZX[ii][jj] = temp;
			cout << RotationZX[ii][jj] << "\t";
		}
		cout << endl;
	}

	cout << "RotationZXZp" << endl;
	RotationZXZp = Tptr2<double>(3, 3);
	for (int ii = 0; ii < 3; ii++)
	{
		for (int jj = 0; jj < 3; jj++)
		{
			double temp = 0;
			for (int ll = 0; ll < 3; ll++)
			{
				temp += RotationZX[ii][ll] * RotationZp[ll][jj];
			}
			//cout << temp << "\t";
			RotationZXZp[ii][jj] = temp;
			cout << RotationZXZp[ii][jj] << "\t";
		}
		cout << endl;
		//cout << endl;
	}

}