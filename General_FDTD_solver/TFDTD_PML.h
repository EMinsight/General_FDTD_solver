//---------------------------------------------------------------------------
#ifndef TFDTD_PMLH
#define TFDTD_PMLH
//---------------------------------------------------------------------------
#define _USE_MATH_DEFINES
#include <iostream>
#include <math.h>
#include "TFunctions.h"
#include "TFDTDparameters.h"

class TFDTD_PML : public TFDTDparameters
{
    public:
        __fastcall TFDTD_PML();
        virtual ~TFDTD_PML();
    public:
        virtual void Set_Structure_and_Electric_property() {};
        
        void Set_Kappa_coefficient(double alphaT, double kappaT);
        void PML_Initial_Setting();

    protected:
        struct_Field_carrier* pFc;
        int SC, CFS;
        double alpha_max;
        double kappa_max;
    	
        double *C_d, C_h;
        double eps_0, mu_0;
        // PML필드값을 일시적으로 저장.
        double temp_field;
    protected:
        // Field pointer
        double ***Ex, ***Ey, ***Ez;
        double ***Hx, ***Hy, ***Hz;

        // PML field
        double ***fxy1, ***fxz1;
        double ***fxy2, ***fxz2;
        double ***fyz1, ***fyx1;
        double ***fyz2, ***fyx2;
        double ***fzx1, ***fzy1;
        double ***fzx2, ***fzy2;
        double ***gxy1, ***gxz1;
        double ***gxy2, ***gxz2;
        double ***gyz1, ***gyx1;
        double ***gyz2, ***gyx2;
        double ***gzx1, ***gzy1;
        double ***gzx2, ***gzy2;

        // PML parameter
        double *Cxa1, *Cya1, *Cza1;
        double *Cxa2, *Cya2, *Cza2;
        double *Cxb1, *Cyb1, *Czb1;
        double *Cxb2, *Cyb2, *Czb2;

        double *Dxa1, *Dya1, *Dza1;
        double *Dxa2, *Dya2, *Dza2;
        double *Dxb1, *Dyb1, *Dzb1;
        double *Dxb2, *Dyb2, *Dzb2;

    protected: // Field initiate
        void PML_field_initiate();
        void PML_parameter_initiate();
        virtual void PML_parameter_initiateX(){};
        virtual void PML_parameter_initiateY(){};
        virtual void PML_parameter_initiateZ(){};

    protected: // Field update
        virtual void X_NFpml(){};
        virtual void X_NGpml(){};
        virtual void X_PFpml(){};
        virtual void X_PGpml(){};

        virtual void Y_NFpml(){};
        virtual void Y_NGpml(){};
        virtual void Y_PFpml(){};
        virtual void Y_PGpml(){};

        virtual void Z_NFpml(){};
        virtual void Z_NGpml(){};
        virtual void Z_PFpml(){};
        virtual void Z_PGpml(){};

    public:
        void Fx_update();
        void Gx_update();
        void Fy_update();
        void Gy_update();
        void Fz_update();
        void Gz_update();

};
#endif
