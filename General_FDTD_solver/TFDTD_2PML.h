//---------------------------------------------------------------------------
#ifndef TFDTD_2PMLH
#define TFDTD_2PMLH
//---------------------------------------------------------------------------
#include "TFunctions.h"
#include "TFDTD_PML.h"

class TFDTD_2PML : public TFDTD_PML
{
	public:
        TFDTD_2PML(struct_Field_carrier* pFc);
        TFDTD_2PML();
		~TFDTD_2PML();

    protected:
        int i, j, k;    // geometric index
        int ai, aj, ak; // offset index
        int bi, bj, bk;
        double temp;  // temporary memory for prev. F value

    protected: // Field initiate
        //void PML_field_initiate();
        void Set_Structure_and_Electric_property();
        void PML_parameter_initiateX();
        void PML_parameter_initiateY();
        void PML_parameter_initiateZ();
    private: // Field update
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
//    public:
//        void PBC_Z_F_PML();
//        void PBC_Z_G_PML();
//        void PBC_Chiral_Z_F_PML();
//        void PBC_Chiral_Z_G_PML();
    public:
        // PML 영역 테스트
        double ****get_PML_Fxy1() { return &fxy1; };
        double ****get_PML_Fxy2() { return &fxy2; };

        double ****get_PML_Fxz1() { return &fxz1; };
        double ****get_PML_Fxz2() { return &fxz2; };

        double ****get_PML_Fyz1() { return &fyz1; };
        double ****get_PML_Fyz2() { return &fyz2; };

        double ****get_PML_Fyx1() { return &fyx1; };
        double ****get_PML_Fyx2() { return &fyx2; };

        double ****get_PML_Fzx1() { return &fzx1; };
        double ****get_PML_Fzx2() { return &fzx2; };

        double ****get_PML_Fzy1() { return &fzy1; };
        double ****get_PML_Fzy2() { return &fzy2; };

        double ****get_PML_Gxy1() { return &gxy1; };
        double ****get_PML_Gxy2() { return &gxy2; };

        double ****get_PML_Gxz1() { return &gxz1; };
        double ****get_PML_Gxz2() { return &gxz2; };

        double ****get_PML_Gyz1() { return &gyz1; };
        double ****get_PML_Gyz2() { return &gyz2; };

        double ****get_PML_Gyx1() { return &gyx1; };
        double ****get_PML_Gyx2() { return &gyx2; };

        double ****get_PML_Gzx1() { return &gzx1; };
        double ****get_PML_Gzx2() { return &gzx2; };

        double ****get_PML_Gzy1() { return &gzy1; };
        double ****get_PML_Gzy2() { return &gzy2; };
};
#endif

