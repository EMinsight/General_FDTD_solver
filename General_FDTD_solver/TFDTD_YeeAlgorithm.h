//---------------------------------------------------------------------------
#ifndef TFDTD_YeeAlgorithmH
#define TFDTD_YeeAlgorithmH
//---------------------------------------------------------------------------
#include <memory.h>
#include "TFDTDparameters.h"

class TFDTD_YeeAlgorithm : public TFDTDparameters
{
public:
	__fastcall TFDTD_YeeAlgorithm();
	~TFDTD_YeeAlgorithm();
	void Field_Initiate(double*** Ex, double*** Ey, double*** Ez, double*** Hx, double*** Hy, double*** Hz);
protected:
	int i, j, k;
	int mx, my, mz; // 미디아 참조변수
    int index;
	double f0;
    int cont;
public:
	void YEE_Initiate(int mod);
	virtual void E_field_update();
	virtual void J_field_update() {};
    void PBC_E_field_update();
	void set_f0(double f0);
    void set_cont(int cont);

private: // 물질의 loss tangent 를 고려한 sorting
		 // The entire array is sorted by quicksort(A, 0, length(A)-1).
	void Sorting_Material() { quicksort(pMp, 0, num_material - 1); };
    void set_FDTDmod();
protected:
    virtual void E_field_update_for_cavity_contents(){};

};

#endif
