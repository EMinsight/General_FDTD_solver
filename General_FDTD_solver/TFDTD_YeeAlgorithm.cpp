//---------------------------------------------------------------------------
#pragma hdrstop
//#include <boost/thread/thread.hpp>
//#include <boost/fusion/algorithm/transformation/push_back.hpp>
//#include <boost/fusion/include/push_back.hpp>
#include "TFDTD_YeeAlgorithm.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)

__fastcall TFDTD_YeeAlgorithm::TFDTD_YeeAlgorithm()
{
    isFDTD=1;
	this->Ex = Ex;
	this->Ey = Ey;
	this->Ez = Ez; 

	this->Hx = Hx;
	this->Hy = Hy;
	this->Hz = Hz;
}

TFDTD_YeeAlgorithm::~TFDTD_YeeAlgorithm()
{

}

void TFDTD_YeeAlgorithm::set_FDTDmod()
{
    FDTDmod = 10;
}

void TFDTD_YeeAlgorithm::set_f0(double f0)
{
	this->f0 = f0;
}
void TFDTD_YeeAlgorithm::set_cont(int cont){this->cont = cont;}
void TFDTD_YeeAlgorithm::YEE_Initiate(int mod)
{
}

void TFDTD_YeeAlgorithm::Field_Initiate(double ***Ex, double ***Ey, double ***Ez, double ***Hx, double ***Hy, double ***Hz)
{
	this->Ex = Ex;
	this->Ey = Ey;
	this->Ez = Ez;

	this->Hx = Hx;
	this->Hy = Hy;
	this->Hz = Hz;
}

void TFDTD_YeeAlgorithm::E_field_update()
{
	for (i = 0; i < txsize; i++)
	for (j = 1; j < tysize; j++)
	for (k = 1; k < tzsize; k++)
	{
		mx = Media_Ex[i][j][k];
		Ex[i][j][k] = C1[mx] * Ex[i][j][k] + Cb[mx] / dy[k] * (Hz[i][j][k] - Hz[i][j - 1][k]) 
										   - Cb[mx] / dz[k] * (Hy[i][j][k] - Hy[i][j][k - 1]);
	}

	for (i = 1; i < txsize; i++)
	for (j = 0; j < tysize; j++)
	for (k = 1; k < tzsize; k++)
	{
		my = Media_Ey[i][j][k];
		Ey[i][j][k] = C1[my] * Ey[i][j][k] + Cb[my] / dz[j] * (Hx[i][j][k] - Hx[i][j][k - 1]) 
										   - Cb[my] / dx[j] * (Hz[i][j][k] - Hz[i - 1][j][k]);
	}

	for (i = 1; i<txsize; i++)
	for (j = 1; j < tysize; j++)
	for (k = 0; k < tzsize; k++)
	{
		mz = Media_Ez[i][j][k];
		Ez[i][j][k] = C1[mz] * Ez[i][j][k] + Cb[mz] / dx[i] * (Hy[i][j][k] - Hy[i - 1][j][k]) 
										   - Cb[mz] / dy[i] * (Hx[i][j][k] - Hx[i][j - 1][k]);
	}
}

void TFDTD_YeeAlgorithm::PBC_E_field_update()
{
    for(i=0;i<txsize;i++)
    for(j=0;j<tysize+1;j++)
    for(k=1;k<tzsize;k++)
    {
       mx = Media_Ex[i][j][k];
	   Ex[i][j][k] = C1[mx] * Ex[i][j][k] + Cb[mx] / dy[k] * (Hz[i][j][k] - Hz[i][j - 1][k])
		   - Cb[mx] / dz[k] * (Hy[i][j][k] - Hy[i][j][k - 1]);
    }
    for(i=0;i<txsize+1;i++)
    for(j=0;j<tysize;j++)
    for(k=1;k<tzsize;k++)
    {
        my = Media_Ey[i][j][k];
		Ey[i][j][k] = C1[my] * Ey[i][j][k] + Cb[my] / dz[j] * (Hx[i][j][k] - Hx[i][j][k - 1])
			- Cb[my] / dx[j] * (Hz[i][j][k] - Hz[i - 1][j][k]);
    }
    for(i=0;i<txsize+1;i++)
    for(j=0;j<tysize+1;j++)
    for(k=1;k<tzsize;k++)
    {
        mz = Media_Ez[i][j][k];
		Ez[i][j][k] = C1[mz] * Ez[i][j][k] + Cb[mz] / dx[i] * (Hy[i][j][k] - Hy[i - 1][j][k])
			- Cb[mz] / dy[i] * (Hx[i][j][k] - Hx[i][j - 1][k]);
    }
}