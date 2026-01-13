//---------------------------------------------------------------------------
#include "TFDTD_Hybrid_mod_with_MOM.h"
//---------------------------------------------------------------------------
TFDTD_Hybrid_mod_with_MOM::TFDTD_Hybrid_mod_with_MOM(double*** Ex, double*** Ey, double*** Ez, double*** Hx, double*** Hy, double*** Hz)
{
    this->Ex = Ex;
    this->Ey = Ey;
    this->Ez = Ez;

    this->Hx = Hx;
    this->Hy = Hy;
    this->Hz = Hz;

    pi = acos(-1);
    mu0 = 4 * pi * 1e-7;
    eps0 = 1e-9 / 36.0 / pi;
    alpha = 0.5;
    
}
void TFDTD_Hybrid_mod_with_MOM::insert_parameters(struct_Hybrid_source_info *pSource)
{
    
    // Json update 방법으로 교체필요.
    Datalen  = pSource[0].Datalen;
    f0 = pSource[0].f0; // 800e6;
    ds = pSource[0].ds; // 5e-3;

    PML = 10;
    SPX = pSource[0].SPX;
    SPY = pSource[0].SPY;
    SPZ = pSource[0].SPZ;
    
    px_plus    = Tptr1<struct_Surface_Current> (Datalen);
    px_minus   = Tptr1<struct_Surface_Current> (Datalen);
    py_plus    = Tptr1<struct_Surface_Current> (Datalen);
    py_minus   = Tptr1<struct_Surface_Current> (Datalen);
    pz_plus    = Tptr1<struct_Surface_Current> (Datalen);
    pz_minus   = Tptr1<struct_Surface_Current> (Datalen);

	pfXminus = fopen(pSource[0].Fname_Xminus,"r");
	pfXplus  = fopen(pSource[0].Fname_Xplus ,"r");
	pfYminus = fopen(pSource[0].Fname_Yminus,"r");
	pfYplus  = fopen(pSource[0].Fname_Yplus ,"r");
	pfZminus = fopen(pSource[0].Fname_Zminus,"r");
    pfZplus  = fopen(pSource[0].Fname_Zplus ,"r");
}

TFDTD_Hybrid_mod_with_MOM::~TFDTD_Hybrid_mod_with_MOM()
{
    fclose(pfXminus);
    fclose(pfXplus);
    fclose(pfYminus);
    fclose(pfYplus);
    fclose(pfZminus);
    fclose(pfZplus);

    memfree(&px_plus );
    memfree(&px_minus);
    memfree(&py_plus );
    memfree(&py_minus);
    memfree(&pz_plus );
    memfree(&pz_minus);
}
void TFDTD_Hybrid_mod_with_MOM::Initial_condition()
{
    int inc_plane_x, inc_plane_y, inc_plane_z;
	inc_plane_x = B_Tx_ter - B_Tx_ini + 1;
	inc_plane_y = B_Ty_ter - B_Ty_ini + 1;
	inc_plane_z = B_Tz_ter - B_Tz_ini + 1;

	Ex_xyN_correction = Tptr2<double >(inc_plane_x, inc_plane_y);
	Ex_xyP_correction = Tptr2<double >(inc_plane_x, inc_plane_y);
	Ex_xzN_correction = Tptr2<double >(inc_plane_x, inc_plane_z);
	Ex_xzP_correction = Tptr2<double >(inc_plane_x, inc_plane_z);

	Ey_xyN_correction = Tptr2<double >(inc_plane_x, inc_plane_y);
	Ey_xyP_correction = Tptr2<double >(inc_plane_x, inc_plane_y);
	Ey_yzN_correction = Tptr2<double >(inc_plane_y, inc_plane_z);
	Ey_yzP_correction = Tptr2<double >(inc_plane_y, inc_plane_z);

	Ez_xzN_correction = Tptr2<double >(inc_plane_x, inc_plane_z);
	Ez_xzP_correction = Tptr2<double >(inc_plane_x, inc_plane_z);
	Ez_yzN_correction = Tptr2<double >(inc_plane_y, inc_plane_z);
	Ez_yzP_correction = Tptr2<double >(inc_plane_y, inc_plane_z);

	Hx_xyN_correction = Tptr2<double >(inc_plane_x, inc_plane_y);
	Hx_xyP_correction = Tptr2<double >(inc_plane_x, inc_plane_y);
	Hx_xzN_correction = Tptr2<double >(inc_plane_x, inc_plane_z);
	Hx_xzP_correction = Tptr2<double >(inc_plane_x, inc_plane_z);

	Hy_xyN_correction = Tptr2<double >(inc_plane_x, inc_plane_y);
	Hy_xyP_correction = Tptr2<double >(inc_plane_x, inc_plane_y);
	Hy_yzN_correction = Tptr2<double >(inc_plane_y, inc_plane_z);
	Hy_yzP_correction = Tptr2<double >(inc_plane_y, inc_plane_z);

	Hz_xzN_correction = Tptr2<double >(inc_plane_x, inc_plane_z);
	Hz_xzP_correction = Tptr2<double >(inc_plane_x, inc_plane_z);
	Hz_yzN_correction = Tptr2<double >(inc_plane_y, inc_plane_z);
	Hz_yzP_correction = Tptr2<double >(inc_plane_y, inc_plane_z);
}
void TFDTD_Hybrid_mod_with_MOM::set_surface_current_data()
{
    double PX, PY, PZ;
    double Jxr, Jxi, Jyr, Jyi, Jzr, Jzi;
    double Mxr, Mxi, Myr, Myi, Mzr, Mzi;
    double tx, ty, tz;

	for(int ii = 0 ; ii < Datalen ; ii++)
	{
		fscanf(pfXminus, "  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le\n",
							&PX,  &PY,  &PZ, &Jxr, &Jxi, &Jyr, &Jyi, &Jzr, &Jzi, &Mxr, &Mxi, &Myr, &Myi, &Mzr, &Mzi);

        px_minus[ii].PX = round(PX,ds) + 19 + SPX + PML;
        px_minus[ii].PY = round(PY,ds) + 19 + SPY + PML-1;
        px_minus[ii].PZ = round(PZ,ds) + 19 + SPZ + PML-1;

        px_minus[ii].abs_Jx = get_abs(Jxr, Jxi);
        px_minus[ii].abs_Jy = get_abs(Jyr, Jyi);
        px_minus[ii].abs_Jz = get_abs(Jzr, Jzi);

        px_minus[ii].abs_Mx = get_abs(Mxr, Mxi);
        px_minus[ii].abs_My = get_abs(Myr, Myi);
        px_minus[ii].abs_Mz = get_abs(Mzr, Mzi);

        px_minus[ii].phase_Jx = get_phase(Jxr, Jxi);
        px_minus[ii].phase_Jy = get_phase(Jyr, Jyi);
        px_minus[ii].phase_Jz = get_phase(Jzr, Jzi);

        px_minus[ii].phase_Mx = get_phase(Mxr, Mxi);
        px_minus[ii].phase_My = get_phase(Myr, Myi);
        px_minus[ii].phase_Mz = get_phase(Mzr, Mzi);

		fscanf(pfXplus , "  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le\n",
							&PX,  &PY,  &PZ, &Jxr, &Jxi, &Jyr, &Jyi, &Jzr, &Jzi, &Mxr, &Mxi, &Myr, &Myi, &Mzr, &Mzi);

        px_plus[ii].PX = round(PX,ds) + 19 + SPX + PML-2;
        px_plus[ii].PY = round(PY,ds) + 19 + SPY + PML-1;
        px_plus[ii].PZ = round(PZ,ds) + 19 + SPZ + PML-1;

        px_plus[ii].abs_Jx = get_abs(Jxr, Jxi);
        px_plus[ii].abs_Jy = get_abs(Jyr, Jyi);
        px_plus[ii].abs_Jz = get_abs(Jzr, Jzi);

        px_plus[ii].abs_Mx = get_abs(Mxr, Mxi);
        px_plus[ii].abs_My = get_abs(Myr, Myi);
        px_plus[ii].abs_Mz = get_abs(Mzr, Mzi);

        px_plus[ii].phase_Jx = get_phase(Jxr, Jxi);
        px_plus[ii].phase_Jy = get_phase(Jyr, Jyi);
        px_plus[ii].phase_Jz = get_phase(Jzr, Jzi);

        px_plus[ii].phase_Mx = get_phase(Mxr, Mxi);
        px_plus[ii].phase_My = get_phase(Myr, Myi);
        px_plus[ii].phase_Mz = get_phase(Mzr, Mzi);

		fscanf(pfYminus, "  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le\n",
							&PX,  &PY,  &PZ, &Jxr, &Jxi, &Jyr, &Jyi, &Jzr, &Jzi, &Mxr, &Mxi, &Myr, &Myi, &Mzr, &Mzi);

        py_minus[ii].PX = round(PX,ds) + 19 + SPX + PML-1;
        py_minus[ii].PY = round(PY,ds) + 19 + SPY + PML;
        py_minus[ii].PZ = round(PZ,ds) + 19 + SPZ + PML-1;

        py_minus[ii].abs_Jx = get_abs(Jxr, Jxi);
        py_minus[ii].abs_Jy = get_abs(Jyr, Jyi);
        py_minus[ii].abs_Jz = get_abs(Jzr, Jzi);

        py_minus[ii].abs_Mx = get_abs(Mxr, Mxi);
        py_minus[ii].abs_My = get_abs(Myr, Myi);
        py_minus[ii].abs_Mz = get_abs(Mzr, Mzi);

        py_minus[ii].phase_Jx = get_phase(Jxr, Jxi);
        py_minus[ii].phase_Jy = get_phase(Jyr, Jyi);
        py_minus[ii].phase_Jz = get_phase(Jzr, Jzi);

        py_minus[ii].phase_Mx = get_phase(Mxr, Mxi);
        py_minus[ii].phase_My = get_phase(Myr, Myi);
        py_minus[ii].phase_Mz = get_phase(Mzr, Mzi);

		fscanf(pfYplus , "  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le\n",
							&PX,  &PY,  &PZ, &Jxr, &Jxi, &Jyr, &Jyi, &Jzr, &Jzi, &Mxr, &Mxi, &Myr, &Myi, &Mzr, &Mzi);

        py_plus[ii].PX = round(PX,ds) + 19 + SPX + PML-1;
        py_plus[ii].PY = round(PY,ds) + 19 + SPY + PML-2;
        py_plus[ii].PZ = round(PZ,ds) + 19 + SPZ + PML-1;

        py_plus[ii].abs_Jx = get_abs(Jxr, Jxi);
        py_plus[ii].abs_Jy = get_abs(Jyr, Jyi);
        py_plus[ii].abs_Jz = get_abs(Jzr, Jzi);

        py_plus[ii].abs_Mx = get_abs(Mxr, Mxi);
        py_plus[ii].abs_My = get_abs(Myr, Myi);
        py_plus[ii].abs_Mz = get_abs(Mzr, Mzi);

        py_plus[ii].phase_Jx = get_phase(Jxr, Jxi);
        py_plus[ii].phase_Jy = get_phase(Jyr, Jyi);
        py_plus[ii].phase_Jz = get_phase(Jzr, Jzi);

        py_plus[ii].phase_Mx = get_phase(Mxr, Mxi);
        py_plus[ii].phase_My = get_phase(Myr, Myi);
        py_plus[ii].phase_Mz = get_phase(Mzr, Mzi);

		fscanf(pfZminus, "  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le\n",
							&PX,  &PY,  &PZ, &Jxr, &Jxi, &Jyr, &Jyi, &Jzr, &Jzi, &Mxr, &Mxi, &Myr, &Myi, &Mzr, &Mzi);

        pz_minus[ii].PX = round(PX,ds) + 19 + SPX + PML-1;
        pz_minus[ii].PY = round(PY,ds) + 19 + SPY + PML-1;
        pz_minus[ii].PZ = round(PZ,ds) + 19 + SPZ + PML;

        pz_minus[ii].abs_Jx = get_abs(Jxr, Jxi);
        pz_minus[ii].abs_Jy = get_abs(Jyr, Jyi);
        pz_minus[ii].abs_Jz = get_abs(Jzr, Jzi);

        pz_minus[ii].abs_Mx = get_abs(Mxr, Mxi);
        pz_minus[ii].abs_My = get_abs(Myr, Myi);
        pz_minus[ii].abs_Mz = get_abs(Mzr, Mzi);

        pz_minus[ii].phase_Jx = get_phase(Jxr, Jxi);
        pz_minus[ii].phase_Jy = get_phase(Jyr, Jyi);
        pz_minus[ii].phase_Jz = get_phase(Jzr, Jzi);

        pz_minus[ii].phase_Mx = get_phase(Mxr, Mxi);
        pz_minus[ii].phase_My = get_phase(Myr, Myi);
        pz_minus[ii].phase_Mz = get_phase(Mzr, Mzi);

		fscanf(pfZplus , "  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le\n",
                            &PX,  &PY,  &PZ, &Jxr, &Jxi, &Jyr, &Jyi, &Jzr, &Jzi, &Mxr, &Mxi, &Myr, &Myi, &Mzr, &Mzi);

        pz_plus[ii].PX = round(PX,ds) + 19 + SPX + PML-1;
        pz_plus[ii].PY = round(PY,ds) + 19 + SPY + PML-1;
        pz_plus[ii].PZ = round(PZ,ds) + 19 + SPZ + PML-2;

        pz_plus[ii].abs_Jx = get_abs(Jxr, Jxi);
        pz_plus[ii].abs_Jy = get_abs(Jyr, Jyi);
        pz_plus[ii].abs_Jz = get_abs(Jzr, Jzi);

        pz_plus[ii].abs_Mx = get_abs(Mxr, Mxi);
        pz_plus[ii].abs_My = get_abs(Myr, Myi);
        pz_plus[ii].abs_Mz = get_abs(Mzr, Mzi);

        pz_plus[ii].phase_Jx = get_phase(Jxr, Jxi);
        pz_plus[ii].phase_Jy = get_phase(Jyr, Jyi);
        pz_plus[ii].phase_Jz = get_phase(Jzr, Jzi);

        pz_plus[ii].phase_Mx = get_phase(Mxr, Mxi);
        pz_plus[ii].phase_My = get_phase(Myr, Myi);
        pz_plus[ii].phase_Mz = get_phase(Mzr, Mzi);
	}

    B_Tx_ini = px_minus[0].PX;
    B_Tx_ter = px_plus[0].PX;

    B_Ty_ini = py_minus[0].PY;
    B_Ty_ter = py_plus[0].PY;

    B_Tz_ini = pz_minus[0].PZ;
    B_Tz_ter = pz_plus[0].PZ;

    ii_min = B_Tx_ini+1;
    ii_max = B_Tx_ter-1;
    jj_min = B_Ty_ini+1;
    jj_max = B_Ty_ter-1;
    kk_min = B_Tz_ini+1;
    kk_max = B_Tz_ter-1;

    Initial_condition();
}

void TFDTD_Hybrid_mod_with_MOM::Apply_sourceE()
{
    Apply_xE();
    Apply_yE();
    Apply_zE();
}

void TFDTD_Hybrid_mod_with_MOM::Apply_sourceH()
{
    Apply_xH();
    Apply_yH();
    Apply_zH();
}

void TFDTD_Hybrid_mod_with_MOM::Set_In_zeroE()
{
	for(int ii = ii_min ; ii < ii_max; ii++)
	for(int jj = jj_min ; jj < jj_max; jj++)
	for(int kk = kk_min ; kk < kk_max; kk++)
    {
        Ex[ii][jj][kk] = 0;
        Ey[ii][jj][kk] = 0;
        Ez[ii][jj][kk] = 0;
    }
}

void TFDTD_Hybrid_mod_with_MOM::Set_In_zeroH()
{
	for(ii = ii_min ; ii < ii_max-1; ii++)
	for(jj = jj_min ; jj < jj_max-1; jj++)
	for(kk = kk_min ; kk < kk_max-1; kk++)
    {
        Hx[ii][jj][kk] = 0;
        Hy[ii][jj][kk] = 0;
        Hz[ii][jj][kk] = 0;
    }
}

void TFDTD_Hybrid_mod_with_MOM::Apply_xE()
{
    for(int ll = 0 ; ll < Datalen ; ll++)
    {
        ii = px_minus[ll].PX;
        jj = px_minus[ll].PY;
        kk = px_minus[ll].PZ;

        if (jj <= B_Ty_ter-1)
            Ey[ii][jj][kk] = Ey[ii][jj][kk] - Ey_yzN_correction[jj - B_Ty_ini][kk - B_Tz_ini];
        if (kk <= B_Tz_ter)
            Ez[ii][jj][kk] = Ez[ii][jj][kk] - Ez_yzN_correction[jj - B_Ty_ini][kk - B_Tz_ini];

        ii = px_plus[ll].PX;
        if (jj <= B_Ty_ter-1)
            Ey[ii][jj][kk] = Ey[ii][jj][kk]  - Ey_yzP_correction[jj - B_Ty_ini][kk - B_Tz_ini];
        if (kk <= B_Tz_ter-1)
            Ez[ii][jj][kk] = Ez[ii][jj][kk]  - Ez_yzP_correction[jj - B_Ty_ini][kk - B_Tz_ini];
    }
}

void TFDTD_Hybrid_mod_with_MOM::Apply_xH()
{
    double  abs, phase;
    for(int ll = 0 ; ll < Datalen ; ll++)
    {
        ii = px_minus[ll].PX;
        jj = px_minus[ll].PY;
        kk = px_minus[ll].PZ;

        if (kk <= B_Tz_ter-1)
            Hy[ii][jj][kk] = Hy[ii][jj][kk]  - Hy_yzN_correction[jj - B_Ty_ini][kk - B_Tz_ini];
        if (jj <= B_Ty_ter-1)
            Hz[ii][jj][kk] = Hz[ii][jj][kk]  - Hz_yzN_correction[jj - B_Ty_ini][kk - B_Tz_ini];

        ii = px_plus[ll].PX;
        if (kk <= B_Tz_ter-1)
            Hy[ii][jj][kk] = Hy[ii][jj][kk] - Hy_yzP_correction[jj - B_Ty_ini][kk - B_Tz_ini];
        if (jj <= B_Ty_ter-1)
            Hz[ii][jj][kk] = Hz[ii][jj][kk] - Hz_yzP_correction[jj - B_Ty_ini][kk - B_Tz_ini];
    }
}

void TFDTD_Hybrid_mod_with_MOM::Apply_yE()
{
    int ii, jj, kk;
    double  abs, phase;
    for(int ll = 0 ; ll < Datalen ; ll++)
    {
        ii = py_minus[ll].PX;
        jj = py_minus[ll].PY;
        kk = py_minus[ll].PZ;

        if (ii <= B_Tx_ter-1)
            Ex[ii][jj][kk] = Ex[ii][jj][kk] - Ex_xzN_correction[ii - B_Tx_ini][kk - B_Tz_ini];
        if (kk <= B_Tz_ter-1)
            Ez[ii][jj][kk] = Ez[ii][jj][kk] - Ez_xzN_correction[ii - B_Tx_ini][kk - B_Tz_ini];

        jj = py_plus[ll].PY;
        if (ii <= B_Tx_ter-1)
            Ex[ii][jj][kk] = Ex[ii][jj][kk] - Ex_xzP_correction[ii - B_Tx_ini][kk - B_Tz_ini];
        if (kk <= B_Tz_ter-1)
            Ez[ii][jj][kk] = Ez[ii][jj][kk] - Ez_xzP_correction[ii - B_Tx_ini][kk - B_Tz_ini];
    }
}

void TFDTD_Hybrid_mod_with_MOM::Apply_yH()
{
    for(int ll = 0 ; ll < Datalen ; ll++)
    {
        ii = py_minus[ll].PX;
        jj = py_minus[ll].PY;
        kk = py_minus[ll].PZ;

        if (kk <= B_Tz_ter-1)
            Hx[ii][jj][kk] = Hx[ii][jj][kk] - Hx_xzN_correction[ii - B_Tx_ini][kk - B_Tz_ini];
        if (ii <= B_Tx_ter-1)
            Hz[ii][jj][kk] = Hz[ii][jj][kk] - Hz_xzN_correction[ii - B_Tx_ini][kk - B_Tz_ini];

        jj = py_plus[ll].PY;
        if (kk <= B_Tz_ter-1)
            Hx[ii][jj][kk] = Hx[ii][jj][kk] - Hx_xzP_correction[ii - B_Tx_ini][kk - B_Tz_ini];
        if (ii <= B_Tx_ter-1)
            Hz[ii][jj][kk] = Hz[ii][jj][kk] - Hz_xzP_correction[ii - B_Tx_ini][kk - B_Tz_ini];
     }
}

void TFDTD_Hybrid_mod_with_MOM::Apply_zE()
{
    for(int ll = 0 ; ll < Datalen ; ll++)
    {
        ii = pz_minus[ll].PX;
        jj = pz_minus[ll].PY;

        kk = pz_minus[ll].PZ;
        if (ii <= B_Tx_ter-1)
            Ex[ii][jj][kk] = Ex[ii][jj][kk] - Ex_xyN_correction[ii - B_Tx_ini][jj - B_Ty_ini];
        if (jj <= B_Ty_ter-1)
            Ey[ii][jj][kk] = Ey[ii][jj][kk] - Ey_xyN_correction[ii - B_Tx_ini][jj - B_Ty_ini];

        kk = pz_plus[ll].PZ;
        if (ii <= B_Tx_ter-1)
            Ex[ii][jj][kk] = Ex[ii][jj][kk] - Ex_xyP_correction[ii - B_Tx_ini][jj - B_Ty_ini];
        if (jj <= B_Ty_ter-1)
            Ey[ii][jj][kk] = Ey[ii][jj][kk] - Ey_xyP_correction[ii - B_Tx_ini][jj - B_Ty_ini];
     }
}

void TFDTD_Hybrid_mod_with_MOM::Apply_zH()
{
    for(int ll = 0 ; ll < Datalen ; ll++)
    {
        ii = pz_minus[ll].PX;
        jj = pz_minus[ll].PY;

        kk = pz_minus[ll].PZ;
        if (jj <= B_Ty_ter-1)
            Hx[ii][jj][kk] = Hx[ii][jj][kk] - Hx_xyN_correction[ii - B_Tx_ini][jj - B_Ty_ini];
        if (ii <= B_Tx_ter-1)
            Hy[ii][jj][kk] = Hy[ii][jj][kk] - Hy_xyN_correction[ii - B_Tx_ini][jj - B_Ty_ini];

        kk = pz_plus[ll].PZ;
        if (jj <= B_Ty_ter-1)
            Hx[ii][jj][kk] = Hx[ii][jj][kk] - Hx_xyP_correction[ii - B_Tx_ini][jj - B_Ty_ini];
        if (ii <= B_Tx_ter-1)
            Hy[ii][jj][kk] = Hy[ii][jj][kk] - Hy_xyP_correction[ii - B_Tx_ini][jj - B_Ty_ini];
    }
}

void TFDTD_Hybrid_mod_with_MOM::Calculate_Incident_Fiend()
{
    // x+- plane
    for(int ll = 0 ; ll < Datalen ; ll++)
    {
        jj = px_minus[ll].PY;
        kk = px_minus[ll].PZ;
        if (jj <= B_Ty_ter-1)
        {
            Ey_yzN_correction[jj - B_Ty_ini][kk - B_Tz_ini] = px_minus[ll].abs_Jy*dt/((double )ds)/eps0 * sinrE(px_minus[ll].phase_Jy);
            Ey_yzP_correction[jj - B_Ty_ini][kk - B_Tz_ini] = px_plus[ll].abs_Jy*dt/((double )ds)/eps0 * sinrE(px_plus[ll].phase_Jy);
        }
        if (kk <= B_Tz_ter-1)
        {
            Ez_yzN_correction[jj - B_Ty_ini][kk - B_Tz_ini] = px_minus[ll].abs_Jz*dt/((double )ds)/eps0 * sinrE(px_minus[ll].phase_Jz);
            Ez_yzP_correction[jj - B_Ty_ini][kk - B_Tz_ini] = px_plus[ll].abs_Jz*dt/((double )ds)/eps0 * sinrE(px_plus[ll].phase_Jz);
        }
    }

    // y+- plane
    for(int ll = 0 ; ll < Datalen ; ll++)
    {
        ii = py_minus[ll].PX;
        kk = py_minus[ll].PZ;
        if (ii <= B_Tx_ter-1)
        {
        	Ex_xzN_correction[ii - B_Tx_ini][kk - B_Tz_ini] = py_minus[ll].abs_Jx*dt/((double )ds)/eps0 * sinrE(py_minus[ll].phase_Jx);
        	Ex_xzP_correction[ii - B_Tx_ini][kk - B_Tz_ini] = py_plus[ll].abs_Jx*dt/((double )ds)/eps0 * sinrE(py_plus[ll].phase_Jx);
        }
        if (kk <= B_Tz_ter-1)
        {
        	Ez_xzN_correction[ii - B_Tx_ini][kk - B_Tz_ini] = py_minus[ll].abs_Jz*dt/((double )ds)/eps0 * sinrE(py_minus[ll].phase_Jz);
        	Ez_xzP_correction[ii - B_Tx_ini][kk - B_Tz_ini] = py_plus[ll].abs_Jz*dt/((double )ds)/eps0 * sinrE(py_plus[ll].phase_Jz);
        }
    }
    // z+- plane
    for(int ll = 0 ; ll < Datalen ; ll++)
    {
        ii = pz_minus[ll].PX;
        jj = pz_minus[ll].PY;

        if (ii <= B_Tx_ter-1)
        {
            Ex_xyN_correction[ii - B_Tx_ini][jj - B_Ty_ini] = pz_minus[ll].abs_Jx*dt/((double )ds)/eps0 * sinrE(pz_minus[ll].phase_Jx);
            Ex_xyP_correction[ii - B_Tx_ini][jj - B_Ty_ini] = pz_plus[ll].abs_Jx*dt/((double )ds)/eps0 * sinrE(pz_plus[ll].phase_Jx);
        }
        if (jj <= B_Ty_ter-1)
        {
            Ey_xyN_correction[ii - B_Tx_ini][jj - B_Ty_ini] = pz_minus[ll].abs_Jy*dt/((double )ds)/eps0 * sinrE(pz_minus[ll].phase_Jy);
            Ey_xyP_correction[ii - B_Tx_ini][jj - B_Ty_ini] = pz_plus[ll].abs_Jy*dt/((double )ds)/eps0 * sinrE(pz_plus[ll].phase_Jy);
        }
    }

    for(int ll = 0 ; ll < Datalen ; ll++)
    {
        jj = px_minus[ll].PY;
        kk = px_minus[ll].PZ;

        if (jj <= B_Ty_ter-1)
        {
        	Hz_yzN_correction[jj - B_Ty_ini][kk - B_Tz_ini] = px_minus[ll].abs_Mz*dt/((double )ds)/mu0 * sinrH(px_minus[ll].phase_Mz);
        	Hz_yzP_correction[jj - B_Ty_ini][kk - B_Tz_ini] = px_plus[ll].abs_Mz*dt/((double )ds)/mu0 * sinrH(px_plus[ll].phase_Mz);;
        }
        if (kk <= B_Tz_ter-1)
        {
        	Hy_yzN_correction[jj - B_Ty_ini][kk - B_Tz_ini] = px_minus[ll].abs_My*dt/((double )ds)/mu0 * sinrH(px_minus[ll].phase_My);
        	Hy_yzP_correction[jj - B_Ty_ini][kk - B_Tz_ini] = px_plus[ll].abs_My*dt/((double )ds)/mu0 * sinrH(px_plus[ll].phase_My);
        }
    }

    for(int ll = 0 ; ll < Datalen ; ll++)
    {
        ii = py_minus[ll].PX;
        kk = py_minus[ll].PZ;

        if (kk <= B_Tz_ter-1)
        {
            Hx_xzN_correction[ii - B_Tx_ini][kk - B_Tz_ini] = py_minus[ll].abs_Mx*dt/((double )ds)/mu0 * sinrH(py_minus[ll].phase_Mx);
            Hx_xzP_correction[ii - B_Tx_ini][kk - B_Tz_ini] = py_plus[ll].abs_Mx*dt/((double )ds)/mu0 * sinrH(py_plus[ll].phase_Mx);
        }
        if (ii <= B_Tx_ter-1)
        {
            Hz_xzN_correction[ii - B_Tx_ini][kk - B_Tz_ini] = py_minus[ll].abs_Mz*dt/((double )ds)/mu0 * sinrH(py_minus[ll].phase_Mz);
            Hz_xzP_correction[ii - B_Tx_ini][kk - B_Tz_ini] = py_plus[ll].abs_Mz*dt/((double )ds)/mu0 * sinrH(py_plus[ll].phase_Mz);
        }
    }

    for(int ll = 0 ; ll < Datalen ; ll++)
    {
        ii = pz_minus[ll].PX;
        jj = pz_minus[ll].PY;

        if (jj <= B_Ty_ter-1)
        {
            Hx_xyN_correction[ii - B_Tx_ini][jj - B_Ty_ini] = pz_minus[ll].abs_Mx*dt/((double )ds)/mu0 * sinrH(pz_minus[ll].phase_Mx);
            Hx_xyP_correction[ii - B_Tx_ini][jj - B_Ty_ini] = pz_plus[ll].abs_Mx*dt/((double )ds)/mu0 * sinrH( pz_plus[ll].phase_Mx);
        }
        if (ii <= B_Tx_ter-1)
        {
            Hy_xyN_correction[ii - B_Tx_ini][jj - B_Ty_ini] = pz_minus[ll].abs_My*dt/((double )ds)/mu0 * sinrH(pz_minus[ll].phase_My);
            Hy_xyP_correction[ii - B_Tx_ini][jj - B_Ty_ini] = pz_plus[ll].abs_My*dt/((double )ds)/mu0 * sinrH(pz_plus[ll].phase_My);
        }
    }
}


double  TFDTD_Hybrid_mod_with_MOM::sinrE(double  phase)
{
    double  Pulse;
    if (ss*dt < alpha*1.0/f0)
        Pulse = 0.5*(1.0-cos(pi*f0*dt*ss/alpha))*sin(2.0*pi*f0*dt*ss + phase);
    else
        Pulse = sin(2.0*pi*f0*dt*ss + phase);
    return Pulse;
}

double  TFDTD_Hybrid_mod_with_MOM::sinrH(double  phase)
{
    double  Pulse;
    if (ss*dt < alpha*1.0/f0)
        Pulse = 0.5*(1.0-cos(pi*f0*dt*(ss+0.5)/alpha))*sin(2.0*pi*f0*dt*(ss+0.5) + phase);
    else
        Pulse = sin(2.0*pi*f0*dt*(ss+0.5) + phase);
    return Pulse;
}

double  TFDTD_Hybrid_mod_with_MOM::get_abs(double re, double im)
{
    double abs;
    abs = sqrt(re*re+im*im);
    return (double )abs;
}

double  TFDTD_Hybrid_mod_with_MOM::get_phase(double y, double x)
{
    double phase;
    if ( x > 0)
    {
        phase = 2*atan(y/(sqrt(x*x+y*y)+x));
    }
    else if (x <= 0 && y != 0)
    {
        phase = 2*atan( (sqrt(x*x+y*y) - x) / y);
    }
    else if(x < 0 && y == 0)
    {
        phase = pi;
    }
    else if (x == 0 && y == 0)
    {
        phase = 0;
    }

    return (double )phase;
}

int TFDTD_Hybrid_mod_with_MOM::round(double nominator, double denominator)
{
    int temp;

    if (nominator >= 0)
    {
        temp = (int)(nominator/denominator + 0.5);
    }
    else if(nominator <0 )
    {
        temp = (int)(nominator/denominator - 0.5);
    }
    return temp;
}