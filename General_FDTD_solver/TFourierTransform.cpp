#include "TFourierTransform.h"

TFourierTransform::TFourierTransform() 
{
	pi = acos(-1);
}
TFourierTransform::~TFourierTransform()
{
	memfree(&pFaxis);
	memfree(&pAbs);
	memfree(&pAngle);
	memfree(&Output);
}

void TFourierTransform::setParameters(double sFre, double tFre, int sampling, double dt)
{
	this->f0 = f0;
	this->sFre = sFre;
	this->tFre = tFre;
	this->dt = dt;
	Frequency_len = sampling;
	this->df = (tFre - sFre) / sampling;

	pFaxis = Tptr1<double>(Frequency_len);
	pAbs = Tptr1<double>(Frequency_len);
	pAngle = Tptr1<double>(Frequency_len);
	pRe = Tptr1<double>(Frequency_len);
	pIm = Tptr1<double>(Frequency_len);
	Output = Tptr1<complex <double> >(Frequency_len);
}

void TFourierTransform::onthefly_dft_single(double tData, int ss)
{
	int i;
	double tr = tData;
	double ti = 0;
	memset(pRe, 0, sizeof(double));
	memset(pIm, 0, sizeof(double));
	double sum_re, sum_im, temp, ctemp, stemp;
	sum_re = sum_im = 0;
	temp = 0.0; ctemp = 0.0; stemp = 0.0;

	temp = -2 * pi*((double)f0*(ss*dt));
	ctemp = cos(temp); stemp = sin(temp);
	sum_re = sum_re + (tr * ctemp + ti * stemp);
	sum_im = sum_im + (-ti * ctemp + tr * stemp);

	pRe[0] = sum_re;
	pIm[0] = sum_im;
}

void TFourierTransform::onthefly_dft(double tData, int ss)
{
	int i;
	double tr = tData;
	double ti = 0;
	memset(pRe,     0, sizeof(double)*Frequency_len);
	memset(pIm,     0, sizeof(double)*Frequency_len);
	double sum_re, sum_im, temp, ctemp, stemp;
	for (i = 0; i < Frequency_len; i++)
	{
		sum_re = sum_im = 0;
		temp = 0.0; ctemp = 0.0; stemp = 0.0;

		temp = -2 * pi*((double)(sFre + df*i)*(ss*dt));
		ctemp = cos(temp); stemp = sin(temp);
		sum_re = sum_re + (tr * ctemp + ti * stemp);
		sum_im = sum_im + (-ti * ctemp + tr * stemp);

		pRe[i] = sum_re;
		pIm[i] = sum_im;
	}
}

void TFourierTransform::FourierTransform(double *tData, int Time_len)
{
	int i, x;
	double *tr; tr = Tptr1<double>(Time_len);
	double *ti; ti = Tptr1<double>(Time_len);
	memcpy(tr, tData, sizeof(double)*Time_len);
	memset(ti, 0, sizeof(double)*Time_len);
	memset(pRe, 0, sizeof(double)*Frequency_len);
	memset(pIm, 0, sizeof(double)*Frequency_len);
	double sum_re, sum_im, temp, ctemp, stemp;
	for (i = 0; i < Frequency_len; i++)
	{
		sum_re = sum_im = 0;
		temp = 0.0; ctemp = 0.0; stemp = 0.0;
		for (x = 0; x < Time_len; x++)
		{
			temp = 2 * dir*pi*((double)(sFre + df*i)*(x*dt));
			ctemp = cos(temp); stemp = sin(temp);
			sum_re = sum_re + (tr[x] * ctemp + ti[x] * stemp);
			sum_im = sum_im + (-ti[x] * ctemp + tr[x] * stemp);
		}
		pRe[i] = sum_re;
		pIm[i] = sum_im;
	}
	delete[] tr;
	delete[] ti;
}

double *TFourierTransform::GetFrequencyAxisData()
{
	for (int ii = 0; ii < Frequency_len; ii++)
	{
		pFaxis[ii] = sFre + df*ii;
	}
	return pFaxis;
}

double *TFourierTransform::GetFrequencyDataArray()
{

	for (int ii = 0; ii < Frequency_len; ii++)
	{
		pAbs[ii] = sqrt(pRe[ii] * pRe[ii] + pIm[ii] * pIm[ii]);
	}
	return pAbs;
}

double *TFourierTransform::GetPhaseDataArray()
{
	for (int ii = 0; ii < Frequency_len; ii++)
	{
		pAngle[ii] = atan(pIm[ii] / pRe[ii]); // radian
	}
	return pAngle;
}

complex<double>* TFourierTransform::GetFourierResults()
{
	memset(Output, 0, sizeof(complex<double>)*Frequency_len);
	for (int ii = 0; ii < Frequency_len; ii++)
	{
		Output[ii] = complex<double>(pRe[ii], pIm[ii]);
	}
	return Output;
}
/*
TEffectiveConstant::TEffectiveConstant() {}

TEffectiveConstant::~TEffectiveConstant()
{
	memfree(&Eeff);
	memfree(&Meff);
	memfree(&pFreq);
	memfree(&pS11);
	memfree(&pS21);
	memfree(&pNeffKD);
	memfree(&pKeffKD);
	memfree(&pNatural_e);
	memfree(&pNKD);
}

void TEffectiveConstant::Set_data_array(double *freq, complex <double> *S11, complex <double> *S21)
{
	memcpy(pS11, S11, sizeof(complex <double>)*array_length);
	memcpy(pS21, S21, sizeof(complex <double>)*array_length);
	for (int ii = 0; ii < array_length; ii++)
		pFreq[ii] = complex<double>(*(freq + ii), 0);

	deff = complex <double>(0.017e-3, 0);
}

void TEffectiveConstant::Set_frelength(int len)
{
	array_length = len;
	pFreq = Tptr1< complex <double> >(array_length);
	pS11 = Tptr1< complex <double> >(array_length);
	pS21 = Tptr1< complex <double> >(array_length);
	pZ_eff = Tptr1< complex <double> >(array_length);
	pNeffKD = Tptr1< complex <double> >(array_length);
	pKeffKD = Tptr1< complex <double> >(array_length);
	pNatural_e = Tptr1< complex <double> >(array_length);
	pNKD = Tptr1< complex <double> >(array_length);
}

void TEffectiveConstant::KramersKronig()
{
	complex <double> t_numerator, t_denominator;
	complex <double> R01(0, 0);
	complex <double> complex_re(1, 0);
	complex <double> complex_im(0, 1);
	complex <double> KD(0, 0);
	complex <double> temp1(0, 0);
	complex <double> temp2(0, 0);
	complex <double> Psi_ij1(0, 0);
	complex <double> Psi_ij2(0, 0);
	complex <double> Psi_ijN(0, 0);
	complex <double> Psi_ijD(0, 0);

	for (int ii = 0; ii < array_length; ii++)
	{
		t_numerator = (complex_re + pS11[ii])*conj(complex_re + pS11[ii]) - pS21[ii] * conj(pS21[ii]);
		t_denominator = (complex_re - pS11[ii])*conj(complex_re - pS11[ii]) - pS21[ii] * conj(pS21[ii]);
		pZ_eff[ii] = sqrt(t_numerator / t_denominator);
	}
	for (int ii = 0; ii < array_length; ii++)
	{
		R01 = (pZ_eff[ii] - complex_re) / (pZ_eff[ii] + complex_re);
		t_numerator = pS21[ii];
		t_denominator = (complex_re - pS11[ii] * R01);
		pNatural_e[ii] = sqrt(t_numerator / t_denominator);
		if (pNatural_e[ii].real() == 0)
			continue;
		pNKD[ii] = log(pNatural_e[ii]);
		//k0 = complex<double>(2.0,0) * pi*pFreq[ii].real() / complex<double>(3e8,0);
		k0 = 2.0 * pi*pFreq[ii].real() / complex<double>(3e8, 0);
		KD = k0*deff;
		pKeffKD[ii] = -pNKD[ii] / KD;

		temp1 = complex <double>(0, 0);
		for (int jj = 0; jj < ii - 2; jj++)
		{
			Psi_ijN = (2 * pi*pFreq[jj].real() * pKeffKD[jj].real());
			Psi_ijD = (2 * pi*pFreq[jj].real())*(2 * pi*pFreq[jj].real()) -
				(2 * pi*pFreq[ii].real())*(2 * pi*pFreq[ii].real());
			if (Psi_ijD.real() == 0)// && Psi_ijD.real() == 0)
				continue;
			Psi_ij1 = Psi_ijN / Psi_ijD;

			Psi_ijN = 2 * pi*pFreq[jj + 1].real() * pKeffKD[jj + 1].real();
			Psi_ijD = (2 * pi*pFreq[jj + 1].real())*(2 * pi*pFreq[jj + 1].real()) -
				(2 * pi*pFreq[ii].real())*(2 * pi*pFreq[ii].real());
			if (Psi_ijD.real() == 0)// && Psi_ijD.real() == 0)
				continue;
			Psi_ij2 = Psi_ijN / Psi_ijD;

			temp1 = temp1 + Psi_ij1 + Psi_ij2;

		}
		temp2 = complex <double>(0, 0);
		for (int jj = ii; jj < array_length - 1; jj++)
		{
			Psi_ijN = 2 * pi*pFreq[jj].real() * pKeffKD[jj].real();
			Psi_ijD = (2 * pi*pFreq[jj].real()) * (2 * pi*pFreq[jj].real()) - (2 * pi*pFreq[ii].real()) *(2 * pi*pFreq[ii].real());
			if (Psi_ijD.real() == 0 && Psi_ijD.real() == 0)
				continue;
			Psi_ij1 = Psi_ijN / Psi_ijD;

			Psi_ijN = 2 * pi*pFreq[jj + 1].real() * pKeffKD[jj + 1].real();
			Psi_ijD = (2 * pi*pFreq[jj + 1].real()) * (2 * pi*pFreq[jj + 1].real()) - (2 * pi*pFreq[ii].real()) * (2 * pi*pFreq[ii].real());
			if (Psi_ijD.real() == 0 && Psi_ijD.real() == 0)
				continue;
			Psi_ij2 = Psi_ijN / Psi_ijD;

			temp2 = temp2 + Psi_ij1 + Psi_ij2;
		}
		pNeffKD[ii] = complex <double>(1, 0) + complex <double>(2, 0)*(pFreq[2] - pFreq[1]) * (temp1.real() + temp2.real());
		printf("%lf\n", pNeffKD[ii].real());
	}
}

double *TEffectiveConstant::Get_Effective_permittivity()
{
	Eeff = Tptr1<double>(array_length);

	for (int ii = 0; ii < array_length; ii++)
	{
		Eeff[ii] = pNeffKD[ii].real() / pZ_eff[ii].real();
	}
	return Eeff;
}

double *TEffectiveConstant::Get_Effective_permeability()
{
	Meff = Tptr1<double>(array_length);

	for (int ii = 0; ii < array_length; ii++)
	{
		Meff[ii] = pNeffKD[ii].real() * pZ_eff[ii].real();
	}
	return Meff;
}

*/