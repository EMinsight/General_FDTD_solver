//---------------------------------------------------------------------------

#ifndef TFDTD_Chiral_ContentsH
#define TFDTD_Chiral_ContentsH
//---------------------------------------------------------------------------
#include "TFDTD_GenericAlgorithm.h"
#include "TFunctions.h"
class TFDTD_Generic_PBC : public TFDTD_GenericAlgorithm
{
public:
    TFDTD_Generic_PBC();
    ~TFDTD_Generic_PBC();

public:
    void PBC_Save_cal();
    void E_field_update();
    void J_field_update();
    void H_field_update();
    

private:
    void Modeling_YeeCell();
    void set_PBC_Screwed_Parameters();
    double Max_dz, Air_Region_Factor;
    int use_Graded;

    double graded_dz;
    double* g_dz;
    double* g_cb_dz;
    double* g_db_dz;

    double temp_save;
    int mx, my, mz;
    int aa1m1, aa1m2;

private:
    double ch_pml_dz;
    int isFDTD;
};
#endif
