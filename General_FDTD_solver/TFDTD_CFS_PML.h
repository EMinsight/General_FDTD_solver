//---------------------------------------------------------------------------
#ifndef TFDTD_CFS_PMLH
#define TFDTD_CFS_PMLH
//---------------------------------------------------------------------------
#define _USE_MATH_DEFINES
#include <iostream>
#include <math.h>
#include "TFunctions.h"
#include "TFDTD_PML.h"

class TFDTD_CFS_PML : public TFDTD_PML
{
	public:
        TFDTD_CFS_PML(struct_Field_carrier* pFc);
        ~TFDTD_CFS_PML();
    private:
        double alpha_max;
        double kappa_max;
        double sigma_max;
        double temp;
        double m_sigma;
        double m_kappa;
        double eta_0;
        double m;
        double eps_r4pml;
        double mu_r;
        int i, j, k;
        int ai, aj, ak;
        int bi, bj, bk;

    private:
    	void X_NFpml();
    	void X_NGpml();
    	void X_PFpml();
    	void X_PGpml();

    	void Y_NFpml();
    	void Y_NGpml();
    	void Y_PFpml();
    	void Y_PGpml();

		void Z_NFpml();
    	void Z_NGpml();
    	void Z_PFpml();
    	void Z_PGpml();

        void Set_Structure_and_Electric_property();
        void PML_parameter_initiateX();
        void PML_parameter_initiateY();
        void PML_parameter_initiateZ();
};
#endif
