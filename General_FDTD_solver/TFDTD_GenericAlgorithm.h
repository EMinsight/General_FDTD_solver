//---------------------------------------------------------------------------
#ifndef TFDTD_GENERICALGORITHMH
#define TFDTD_GENERICALGORITHMH
//---------------------------------------------------------------------------
#include "TFDTDparameters.h"
#include "TFunctions.h"
#include <math.h>
#include <omp.h>
//#include <openblas/cblas.h>
class TFDTD_GenericAlgorithm : public TFDTDparameters
{
public:
	TFDTD_GenericAlgorithm();
	~TFDTD_GenericAlgorithm();
	void Field_Initiate(double*** Ex, double*** Ey, double*** Ez,
		double**** Jx, double**** Jy, double**** Jz, 
		double*** Hx, double*** Hy, double*** Hz);
	void Field_Initiate(struct_Field_carrier *sFc);
protected:
	int num_thread;
	int num_thread_pef_Fnc;
	int mx, my, mz; // �̵�� ��������
	double*** Ex_p, *** Ey_p, *** Ez_p;
	double*** Ex_pp, *** Ey_pp, *** Ez_pp;
	double**** Jx, **** Jy, **** Jz;
	double**** Jx_p, **** Jy_p, **** Jz_p;
	double tempx, tempy, tempz;
	double*** temp1, *** temp2, *** temp3, *** temp4, *** temp5, *** temp6;
	double f0;
    int index;
	int	select_mod;
public:
	void E_field_update();
	void Ex_previous();
	void Ey_previous();
	void Ez_previous();
	void Ex_field_update();
	void Ey_field_update();
	void Ez_field_update();
	void J_field_update();
	void Jx_field_update();
	void Jy_field_update();
	void Jz_field_update();
	struct_Field_carrier get_Field_address(struct_Field_carrier* Fc);
	//void PBC_E_field_update();
	
private: // �⺻ ���� update �� �̵�� ��� ����.
	void Modeling_YeeCell();
    void Modeling_Meta_Cell();
private:
	// Duke QCRF ��
	// void Derive_Complex_eps();

	double temp_save;
	//double temp1, temp2, temp3, temp4, temp5, temp6;
};
#endif

