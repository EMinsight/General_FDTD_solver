//---------------------------------------------------------------------------

#ifndef TFDTD_2PML_PBCH
#define TFDTD_2PML_PBCH
//---------------------------------------------------------------------------
#include "TFDTD_2PML.h"
class TFDTD_2PML_PBC: public TFDTD_2PML
{
    public:
        TFDTD_2PML_PBC() {};
        ~TFDTD_2PML_PBC();
    protected:
        // struct_Field_carrier pFc;
        double*** Ex_p, *** Ey_p, *** Ez_p;
        double*** Ex_pp, *** Ey_pp, *** Ez_pp;
        double**** Jx, **** Jy, **** Jz;
        double**** Jx_p, **** Jy_p, **** Jz_p;

    public:
		void PBC_Z_F_PML_minus();
		void PBC_Z_G_PML_minus();
		void PBC_Z_F_PML_plus();
		void PBC_Z_G_PML_plus();
       // void PBC_Chiral_Z_F_PML();
       // void PBC_Chiral_Z_G_PML();
};
#endif
