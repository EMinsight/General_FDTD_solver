//---------------------------------------------------------------------------
#ifndef TFDTDparametersH
#define TFDTDparametersH
//---------------------------------------------------------------------------
#include "TFunctions.h"
#include <memory.h>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <string>
#include <complex>
#include <omp.h>

// using namespace std;

void swap(Material_properties* A, Material_properties* B);
inline void quicksort(Material_properties* Array, int lo, int hi);
inline int partition(Material_properties* Array, int lo, int hi);

// --------------------------------------------------------------------------- //
class TFDTDparameters 
{
public:
	__fastcall TFDTDparameters();
//    __fastcall TFDTDparameters(bool CreateSuspended);
	virtual ~TFDTDparameters();
	void Initial_Setting();
	void Initial_Setting(struct_Field_carrier* pFc);
	void Set_Kappa_coefficient(double alphaT, double kappaT);
protected:
	double m_pi;
	char JsonFileName[1024];
	char MediaFileLocation[1024];
public:
	// ����� ���� ��ȯ
	int set_planewave_buffer() { return planewave_buffer; };
	void set_PMLXYZ(int PMLX, int PMLY, int PMLZ);

	void set_PMLX(int PMLX) { this->PMLX =  PMLX;};
	void set_PMLY(int PMLY) { this->PMLY =  PMLY;};
	void set_PMLZ(int PMLZ) { this->PMLZ =  PMLZ;};

	int get_total_pol() { return total_pol; };
    int get_FDTDmod(){ return FDTDmod;};
	Material_properties* get_material_properties() { return pMp; };
protected:
	Document Jdocument;
	void Read_Media_Data();
	void Read_Material_Data();
	//	void virtual E_field_update();
	void Initiate_Fields();
	void Initiate_MediaSpace();
	//	void Modeling_YeeCell();
	void Spatial_Size_With_Kappa();

    int FDTDmod;
    int isFDTD;
    virtual void set_FDTDmod(){};
	
public:
	void H_field_update();
	void Hx_field_update();
	void Hy_field_update();
	void Hz_field_update();
	virtual void Modeling_YeeCell() {};

protected:
	double ***Ex, ***Ey, ***Ez;
	double ***Hx, ***Hy, ***Hz;

protected:
	struct_Field_carrier* Sfc;
	// ����
	int PX, PY, PZ;
	int PMLX, PMLY, PMLZ;
	int txsize, tysize, tzsize;
	double* dx, * dy, * dz;
	double* E_dx, * E_dy, * E_dz;
	double* H_dx, * H_dy, * H_dz;
	double dt;
	double eps_0, mu_0;
	double c_0, eta_0;
	// ����� ����
	int planewave_buffer;

	// ����
	Material_properties* pMp;
	int total_pol = 0;
	double *eps_r_inf;
	double *C1, *C2; // num. of materials
	double** Ca_temp, ** Cb_temp; // Ex - update, J coefficient.
	double *Cb;
	double Db;

	int disp_s;
	int num_material;
	int ***Media;
	int ***Media_Ex, ***Media_Ey, ***Media_Ez, ***Media_Hx, ***Media_Hy, ***Media_Hz;
	int i, j, k, pol;
	int mx, my, mz;
	// kappa
    void Kappa_parameter_initiateX();
    void Kappa_parameter_initiateY();
    void Kappa_parameter_initiateZ();

    double kappa_max;
    int m_kappa;
    double *kappa_e_xprof;
    double *kappa_h_xprof;

    double *kappa_e_yprof;
    double *kappa_h_yprof;

    double *kappa_e_zprof;
    double *kappa_h_zprof;

protected:
    int n_Thread;
    int n_bind;
};

inline void quicksort(Material_properties* Array, int lo, int hi)
{
	double Pivot;
	if (lo < hi)
	{
		Pivot = partition(Array, lo, hi);
		quicksort(Array, lo, Pivot);
		quicksort(Array, Pivot + 1, hi);
	}

}
inline int partition(Material_properties* Array, int lo, int hi)
{
	double	pivot = (Array + lo)->loss_tan;
	int i = lo - 1;
	int j = hi + 1;
	while (1)
	{
		do
		{
			i++;
		} while ((Array + i)->loss_tan < pivot);

		do
		{
			j--;
		} while ((Array + j)->loss_tan > pivot);

		if (i >= j)
			return j;

		swap((Array + j), (Array + i));
	}
}
inline void swap(Material_properties* A, Material_properties* B)
{
	Material_properties* temp = new Material_properties;
	memcpy(temp, A, sizeof(Material_properties));
	memcpy(A, B, sizeof(Material_properties));
	memcpy(B, temp, sizeof(Material_properties));
	delete temp;
}

#endif
