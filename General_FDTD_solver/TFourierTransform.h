#pragma once
#include <memory>
#include <complex>
#include "TFunctions.h"
// using namespace std;

class TFourierTransform
{
private:
	double pi;
	double tFre;
	double sFre;
	double f0;
public:
	TFourierTransform();
	~TFourierTransform();

private:
	int Frequency_len;
	int Time_len;
	int dir;

	double *pTime;
	double *pTData;
	double *pFaxis;
	double *pRe;
	double *pIm;
	double df;
	double dt;

	double *pAbs;
	double *pAngle;
	complex <double> *Output;
public:
	double *GetFrequencyAxisData();
	double *GetFrequencyDataArray();
	double *GetPhaseDataArray();
	complex<double>*GetFourierResults();
	void setParameters(double sFre, double tFre, int sampling, double dt);

public:
	void FourierTransform(double *tData, int Time_len);
	double * Get_Frequency_Domain_Response_Imag() { return pRe; };
	double * Get_Frequency_Domain_Response_Real() { return pIm; };
	void onthefly_dft(double tData, int ss);
	void onthefly_dft_single(double tData, int ss);
	// 		
public:
	long  PowerValueOf2(int n);
	bool  IsPowerOf2(int n);

	//void DataReArrangement(double *re, double *im, int nLength);
	//void zeropadding();
	//void FastFourierTransform();
};

//
//class TEffectiveConstant
//{
//public:
//	TEffectiveConstant();
//	~TEffectiveConstant();
//public:
//	FILE *pRead;
//	complex <double> *pFreq;
//	complex <double> k0;
//	complex <double> deff;
//	complex <double> *pS11;
//	complex <double> *pS21;
//	complex <double> *pZ_eff;
//	complex <double> *pNatural_e;
//	complex <double> *pNeffKD;
//	complex <double> *pKeffKD;
//	complex <double> *pNKD;
//
//	double *Eeff;
//	double *Meff;
//	int array_length;
//
//public:
//	//void Set_data_array(double *freq, complex <double> *S11, complex <double> *S21);
//	//void Set_frelength(int len);
//	//void KramersKronig();
//	//double *Get_Effective_permittivity();
//	//double *Get_Effective_permeability();
//
//private:
//	TFourierTransform *pFourier;
//};
