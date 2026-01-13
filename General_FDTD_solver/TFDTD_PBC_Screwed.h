#pragma once
#include "TFDTD_PBC_Boundary_normal.h"
#include "TFunctions.h"

class TFDTD_PBC_Screwed : public TFDTD_PBC_Boundary_normal
{
public:
    TFDTD_PBC_Screwed(struct_Field_carrier* pFc);
    TFDTD_PBC_Screwed();
    ~TFDTD_PBC_Screwed() {};

public:
    void FDTD_PBC_Chiral_BC1_omp();
    void FDTD_PBC_Chiral_BC2_omp();

private:

    void set_PBC_Screwed_Parameters();
    double Max_dz, Air_Region_Factor;
    int use_Graded;
   
    double graded_dz;
    double* g_dz;
    double* g_cb_dz;
    double* g_db_dz;

    double temp_save;
    int mx, my, mz;

private:
    double ch_pml_dz;
    int isFDTD;

};
