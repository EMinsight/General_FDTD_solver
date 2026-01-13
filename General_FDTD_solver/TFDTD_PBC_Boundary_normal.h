//---------------------------------------------------------------------------
#ifndef TFDTD_PBC_Boundary_normalH
#define TFDTD_PBC_Boundary_normalH
//---------------------------------------------------------------------------
#include "TFDTD_2PML_PBC.h"
#include "TFunctions.h"

class TFDTD_PBC_Boundary_normal : public TFDTD_2PML_PBC
{
public:
	TFDTD_PBC_Boundary_normal();
	virtual ~TFDTD_PBC_Boundary_normal();
protected:
	// Field pointer
    int aa1m1, aa1m2;
    double _pi;

public:
	void FDTD_PBC_Ex_cal();
	void FDTD_PBC_Ey_cal();
	void FDTD_PBC_Ez_cal();
	void FDTD_PBC_H_cal();

    //void set_media_Ex(int ***media_Ex){this->Medi_Ex = media_Ex;};
    //void set_media_Ey(int ***media_Ey){this->Medi_Ey = media_Ey;};
    //void set_media_Ez(int ***media_Ez){this->Medi_Ez = media_Ez;};

private:
	double temp_save;
	int mx, my, mz;
    

};

#endif
