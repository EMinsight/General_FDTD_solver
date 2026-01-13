#include "TFDTD_Hybrid_TFSF.h"
#include "TSource_Excitation.h"



TFDTD_Hybrid_TFSF::TFDTD_Hybrid_TFSF(struct_Field_carrier* pFc , coordinate coord)
{
    this->pFc = pFc;
    m_pi = acos(-1);
    eps_0 = 1e-9 / 36 / m_pi;
    mu_0 = 4 * m_pi * 1e-7;
    c_0 = 1.0 / sqrt(eps_0 * mu_0);
	eta = sqrt(mu_0 / eps_0);
	alpha = 3.0;   // 원래는 0.5였음. 강한 펄스 신호에 대응하기 위해 3.0으로 변경
    n_T = 3; // 이것도 넘겨받나...? 없을텐데?

    this->Ex = pFc->Ex_address; this->Hx = pFc->Hx_address;
    this->Ey = pFc->Ey_address; this->Hy = pFc->Hy_address;
    this->Ez = pFc->Ez_address; this->Hz = pFc->Hz_address;

    this->dx = pFc->dx;
    this->dy = pFc->dy;
    this->dz = pFc->dz;
    this->dt = pFc->dt;
    this->f0 = pFc->f0;
    this->txsize = pFc->txsize;
    this->tysize = pFc->tysize;
    this->tzsize = pFc->tzsize;

    T_0 = 1 / f0;
    omega = 2 * m_pi * f0;
    beta = omega * sqrt(eps_0 * mu_0);

    // 인덱스 및 거리 계산을 위해 Hybrid class가 좌표값을 직접 알고 있어야 함.
	this->xcoordinate = coord.xcoordinate;
	this->ycoordinate = coord.ycoordinate;
	this->zcoordinate = coord.zcoordinate;

	this->xStaggered = coord.xStaggered;
	this->yStaggered = coord.yStaggered;
	this->zStaggered = coord.zStaggered;


    // 일단 전계와 mesh는 여기에서 넘겨받았다. 
    // 이제 source 파일에서 box 설정을 해야 한다. 
}

TFDTD_Hybrid_TFSF::~TFDTD_Hybrid_TFSF() {}

void TFDTD_Hybrid_TFSF::set_Hybrid_TFSF_box(struct_Hybrid_TFSF_source_info* pSource)
{
    insert_parameters(pSource); // 이 함수를 실행하면?
}



double TFDTD_Hybrid_TFSF::get_plane_constant_x(FILE* fp)
{
    double PX, PY, PZ, RR, II;
    long pos = ftell(fp);           // 현재 위치 저장
    fseek(fp, 0L, SEEK_SET);        // 맨 앞으로
    fscanf(fp, " %le %le %le %le %le", &PX, &PY, &PZ, &RR, &II);
    fseek(fp, pos, SEEK_SET);       // 원래 위치 복구 (안전)
    return PX;
}



double TFDTD_Hybrid_TFSF::circular_mean_phase(struct_Hybrid_TFSF** A, int Ny, int Nz)
{
    printf("Ny: %d , Nz: %d\n", Ny, Nz);
    double s = 0.0, c = 0.0;
    for (int j = 0; j < Ny; ++j)
        for (int k = 0; k < Nz; ++k)
        {
            s += sin(A[j][k].Phase);
            c += cos(A[j][k].Phase);
        }
    return atan2(s, c);
}


void TFDTD_Hybrid_TFSF::set_plane_phase_constant(struct_Hybrid_TFSF** A, int Ny, int Nz, double phi)
{
    for (int j = 0; j < Ny; ++j)
        for (int k = 0; k < Nz; ++k)
            A[j][k].Phase = phi;
}

void TFDTD_Hybrid_TFSF::add_plane_phase_offset(struct_Hybrid_TFSF** A, int Ny, int Nz, double dphi)
{
    for (int j = 0; j < Ny; ++j)
        for (int k = 0; k < Nz; ++k)
            A[j][k].Phase = wrap_to_pi(A[j][k].Phase + dphi);
}



void set_phase_ramp_x_for_E(struct_Hybrid_TFSF** A, int Nx, int N2, double phi_ref,
    double beta, double x_minus, double dx, double x_offset_cell)
{
    // A[x][*] 형태 (두 번째 인덱스는 y 또는 z)
    for (int i = 0; i < Nx; ++i)
    {
        double x = x_minus + (i + x_offset_cell) * dx; // x_offset_cell: 0.0 or 0.5
        double phi = wrap_to_pi(phi_ref + beta * (x - x_minus));
        for (int j = 0; j < N2; ++j)
            A[i][j].Phase = phi;
    }
}

void set_phase_ramp_x_for_H(struct_Hybrid_TFSF** A, int Nx, int N2, double phi_ref,
    double beta, double x_minus, double dx, double x_offset_cell)
{
    // H도 동일하게, 단 x_offset_cell을 H의 Yee 위치에 맞춰 다르게 줄 것
    for (int i = 0; i < Nx; ++i)
    {
        double x = x_minus + (i + x_offset_cell) * dx;
        double phi = wrap_to_pi(phi_ref + beta * (x - x_minus));
        for (int j = 0; j < N2; ++j)
            A[i][j].Phase = phi;
    }
}


void TFDTD_Hybrid_TFSF::insert_parameters(struct_Hybrid_TFSF_source_info* pSource)
{
    pf_Ey_myz.fp = fopen(pSource[0].plane1_Ey, "r");
    pf_Ez_myz.fp = fopen(pSource[0].plane1_Ez, "r");
    pf_Ex_xmz.fp = fopen(pSource[0].plane2_Ex, "r");
    pf_Ez_xmz.fp = fopen(pSource[0].plane2_Ez, "r");
    pf_Ey_pyz.fp = fopen(pSource[0].plane3_Ey, "r");
    pf_Ez_pyz.fp = fopen(pSource[0].plane3_Ez, "r");
    pf_Ex_xpz.fp = fopen(pSource[0].plane4_Ex, "r");
    pf_Ez_xpz.fp = fopen(pSource[0].plane4_Ez, "r");
    pf_Ex_xyp.fp = fopen(pSource[0].plane5_Ex, "r");
    pf_Ey_xyp.fp = fopen(pSource[0].plane5_Ey, "r");
    pf_Ey_xym.fp = fopen(pSource[0].plane6_Ey, "r");
    pf_Ex_xym.fp = fopen(pSource[0].plane6_Ex, "r");

    pf_Hy_myz.fp = fopen(pSource[0].plane1_Hy, "r");
    pf_Hz_myz.fp = fopen(pSource[0].plane1_Hz, "r");
    pf_Hx_xmz.fp = fopen(pSource[0].plane2_Hx, "r");
    pf_Hz_xmz.fp = fopen(pSource[0].plane2_Hz, "r");
    pf_Hy_pyz.fp = fopen(pSource[0].plane3_Hy, "r");
    pf_Hz_pyz.fp = fopen(pSource[0].plane3_Hz, "r");
    pf_Hx_xpz.fp = fopen(pSource[0].plane4_Hx, "r");
    pf_Hz_xpz.fp = fopen(pSource[0].plane4_Hz, "r");
    pf_Hx_xyp.fp = fopen(pSource[0].plane5_Hx, "r");
    pf_Hy_xyp.fp = fopen(pSource[0].plane5_Hy, "r");
    pf_Hy_xym.fp = fopen(pSource[0].plane6_Hy, "r");
    pf_Hx_xym.fp = fopen(pSource[0].plane6_Hx, "r");

    pf_Ey_myz.Datalen = count_line(pf_Ey_myz.fp);
    pf_Ez_myz.Datalen = count_line(pf_Ez_myz.fp); //1640 나옴
    pf_Ex_xmz.Datalen = count_line(pf_Ex_xmz.fp);
    pf_Ez_xmz.Datalen = count_line(pf_Ez_xmz.fp);
    pf_Ey_pyz.Datalen = count_line(pf_Ey_pyz.fp);
    pf_Ez_pyz.Datalen = count_line(pf_Ez_pyz.fp);
    pf_Ex_xpz.Datalen = count_line(pf_Ex_xpz.fp);
    pf_Ez_xpz.Datalen = count_line(pf_Ez_xpz.fp);
    pf_Ex_xyp.Datalen = count_line(pf_Ex_xyp.fp);
    pf_Ey_xyp.Datalen = count_line(pf_Ey_xyp.fp);
    pf_Ey_xym.Datalen = count_line(pf_Ey_xym.fp);
    pf_Ex_xym.Datalen = count_line(pf_Ex_xym.fp);

    pf_Hy_myz.Datalen = count_line(pf_Hy_myz.fp);
    pf_Hz_myz.Datalen = count_line(pf_Hz_myz.fp);
    pf_Hx_xmz.Datalen = count_line(pf_Hx_xmz.fp);
    pf_Hz_xmz.Datalen = count_line(pf_Hz_xmz.fp);
    pf_Hy_pyz.Datalen = count_line(pf_Hy_pyz.fp);
    pf_Hz_pyz.Datalen = count_line(pf_Hz_pyz.fp);
    pf_Hx_xpz.Datalen = count_line(pf_Hx_xpz.fp);
    pf_Hz_xpz.Datalen = count_line(pf_Hz_xpz.fp);
    pf_Hx_xyp.Datalen = count_line(pf_Hx_xyp.fp);
    pf_Hy_xyp.Datalen = count_line(pf_Hy_xyp.fp);
    pf_Hy_xym.Datalen = count_line(pf_Hy_xym.fp);
    pf_Hx_xym.Datalen = count_line(pf_Hx_xym.fp);

    //일단 임시로 설정. 나중에 불균일 격자에 맞게 바꿔야 함.
    double dxx = 0.005;

    Cedx1 = dt / eps_0 / dxx;
    Cedy1 = dt / eps_0 / dxx;
    Cedz1 = dt / eps_0 / dxx;
    Chdx1 = dt / mu_0 / dxx;
    Chdy1 = dt / mu_0 / dxx;
    Chdz1 = dt / mu_0 / dxx;

    Cedx2 = dt / eps_0 / dxx;
    Cedy2 = dt / eps_0 / dxx;
    Cedz2 = dt / eps_0 / dxx;
    Chdx2 = dt / mu_0 / dxx;
    Chdy2 = dt / mu_0 / dxx;
    Chdz2 = dt / mu_0 / dxx;



    Load_tangential_TFSF_fields();
}


void TFDTD_Hybrid_TFSF::Load_tangential_TFSF_fields()
{
    // ---------------------------
    // 0) Allocate & Load 6-face tangential fields from files (phasor Abs/Phase)
    // ---------------------------
    // (이미 너 코드처럼 Calculate_2D_size()가 파일에서 PXmin/PXmax 잡아서 배열 할당+채움)
    Calculate_2D_size_EyYZ(pf_Ey_myz, &Ey_myz, "Ey_myz"); Calculate_2D_size_HyYZ(pf_Hy_myz, &Hy_myz, "Hy_myz");
    Calculate_2D_size_EzYZ(pf_Ez_myz, &Ez_myz, "Ez_myz"); Calculate_2D_size_HzYZ(pf_Hz_myz, &Hz_myz, "Hz_myz");

    Calculate_2D_size_ExXZ(pf_Ex_xmz, &Ex_xmz, "Ex_xmz"); Calculate_2D_size_HxXZ(pf_Hx_xmz, &Hx_xmz, "Hx_xmz");
    Calculate_2D_size_EzXZ(pf_Ez_xmz, &Ez_xmz, "Ez_xmz"); Calculate_2D_size_HzXZ(pf_Hz_xmz, &Hz_xmz, "Hz_xmz");

    Calculate_2D_size_EyYZ(pf_Ey_pyz, &Ey_pyz, "Ey_pyz"); Calculate_2D_size_HyYZ(pf_Hy_pyz, &Hy_pyz, "Hy_pyz");
    Calculate_2D_size_EzYZ(pf_Ez_pyz, &Ez_pyz, "Ez_pyz"); Calculate_2D_size_HzYZ(pf_Hz_pyz, &Hz_pyz, "Hz_pyz");

    Calculate_2D_size_ExXZ(pf_Ex_xpz, &Ex_xpz, "Ex_xpz"); Calculate_2D_size_HxXZ(pf_Hx_xpz, &Hx_xpz, "Hx_xpz");
    Calculate_2D_size_EzXZ(pf_Ez_xpz, &Ez_xpz, "Ez_xpz"); Calculate_2D_size_HzXZ(pf_Hz_xpz, &Hz_xpz, "Hz_xpz");

    Calculate_2D_size_ExXY(pf_Ex_xyp, &Ex_xyp, "Ex_xyp"); Calculate_2D_size_HxXY(pf_Hx_xyp, &Hx_xyp, "Hx_xyp");
    Calculate_2D_size_EyXY(pf_Ey_xyp, &Ey_xyp, "Ey_xyp"); Calculate_2D_size_HyXY(pf_Hy_xyp, &Hy_xyp, "Hy_xyp");

    Calculate_2D_size_EyXY(pf_Ey_xym, &Ey_xym, "Ey_xym"); Calculate_2D_size_HyXY(pf_Hy_xym, &Hy_xym, "Hy_xym");
    Calculate_2D_size_ExXY(pf_Ex_xym, &Ex_xym, "Ex_xym"); Calculate_2D_size_HxXY(pf_Hx_xym, &Hx_xym, "Hx_xym");

    // ---------------------------
    // 1) Compute TF/SF box offsets (node 기준 박스 경계)
    // ---------------------------
    Calculate_start_pointE();
    Calculate_end_pointE();
    
    B_Tx_ini = (int)PXoffsetE;
    B_Ty_ini = (int)PYoffsetE;
    B_Tz_ini = (int)PZoffsetE;

    B_Tx_ter = (int)PXoffsetE2;
    B_Ty_ter = (int)PYoffsetE2;
    B_Tz_ter = (int)PZoffsetE2;

    const int Nx_node = (B_Tx_ter - B_Tx_ini + 1);
    const int Ny_node = (B_Ty_ter - B_Ty_ini + 1);
    const int Nz_node = (B_Tz_ter - B_Tz_ini + 1);

    // analytic 6faces (Yee) 격자 규칙:
    const int Nx_ctr = Nx_node - 1;
    const int Ny_ctr = Ny_node - 1;
    const int Nz_ctr = Nz_node - 1;

    const int Nx_out = Nx_node + 1;
    const int Ny_out = Ny_node + 1;
    const int Nz_out = Nz_node + 1;

    printf("[TFSF box node] Nx=%d Ny=%d Nz=%d\n", Nx_node, Ny_node, Nz_node);
    printf("[ctr] Nx=%d Ny=%d Nz=%d | [out] Nx=%d Ny=%d Nz=%d\n",Nx_ctr, Ny_ctr, Nz_ctr, Nx_out, Ny_out, Nz_out);

    // ---------------------------
    // 2) +x 진행파: yz(-x) 면을 기준 위상(phi_ref)으로 잡고
    //    yz(+x) 면은 phi_ref + beta*(x_plus-x_minus) 관계를 강제
    //    ★ 핵심: Ey/Ez/Hy/Hz는 각각 격자 크기가 다르므로 "성분별 Ny,Nz"를 써야 함
    // ---------------------------
    // yz면 실제 격자 크기 (파일 정의 그대로)
    const int Ny_Ey_yz = Ny_ctr;     // Ey: y_ctr × z_node  -> (Ny-1) × Nz
    const int Nz_Ey_yz = Nz_node;

    const int Ny_Ez_yz = Ny_node;    // Ez: y_node × z_ctr  -> Ny × (Nz-1)
    const int Nz_Ez_yz = Nz_ctr;

    const int Ny_Hy_yz = Ny_out;     // Hy: y_out × z_node  -> (Ny+1) × Nz
    const int Nz_Hy_yz = Nz_node;

    const int Ny_Hz_yz = Ny_node;    // Hz: y_node × z_out  -> Ny × (Nz+1)
    const int Nz_Hz_yz = Nz_out;

    // yz(-x) 기준면 평균 위상 (E만 섞어서 phi_ref 잡기)
    double phi_Ey_m = circular_mean_phase(Ey_myz, Ny_Ey_yz, Nz_Ey_yz);
    double phi_Ez_m = circular_mean_phase(Ez_myz, Ny_Ez_yz, Nz_Ez_yz);

    double phi_ref = atan2(
        sin(phi_Ey_m) + sin(phi_Ez_m),
        cos(phi_Ey_m) + cos(phi_Ez_m)
    );

    

    // yz(+x) 면도 우선 자기 평균으로 평탄화한 뒤, 진행파 위상관계로 dphi를 강제
    double phi_Ey_p = circular_mean_phase(Ey_pyz, Ny_Ey_yz, Nz_Ey_yz);
    double phi_Ez_p = circular_mean_phase(Ez_pyz, Ny_Ez_yz, Nz_Ez_yz);

    double phi_p = atan2(
        sin(phi_Ey_p) + sin(phi_Ez_p),
        cos(phi_Ey_p) + cos(phi_Ez_p)
    );

    

    // yz planes physical x 위치 (파일 첫 줄 PX 기준)
    double x_minus = get_plane_constant_x(pf_Ey_myz.fp); // yz(-x)
    double x_plus = get_plane_constant_x(pf_Ey_pyz.fp); // yz(+x)

    double phi_target = wrap_to_pi(phi_ref + beta * (x_plus - x_minus));
    double dphi = wrap_to_pi(phi_target - phi_p);

    

    // ---------------------------
    // 4) (기존) 시간 지연/속도 파라미터 세팅 (너 코드 유지)
    // ---------------------------
    double dx_plane = x_plus - x_minus;
    printf("yz planes: x_minus=%g, x_plus=%g, dx_plane=%g\n", x_minus, x_plus, dx_plane);

    printf("xini:%d yini:%d zini:%d xter:%d yter:%d zter:%d\n",
        B_Tx_ini, B_Ty_ini, B_Tz_ini, B_Tx_ter, B_Ty_ter, B_Tz_ter);

    // 임시: 기존 값 유지
    V_tilt = (0.99 * min1D(dz)) / dt / sqrt(3);
    printf("Vtilt: %le\n", V_tilt);

    Distance = (B_Tx_ter - B_Tx_ini) * min1D(dx);
    Time0 = 0;
    Time1E = 0; // Distance / V_tilt;  (원래 코드가 이렇게 돼 있던데, 필요하면 복구)
    Time1H = 0; // (Distance + dx) / V_tilt;

    Initial_condition();
}




void  TFDTD_Hybrid_TFSF::set_current_steps(int ss)
{
    this->nsteps = ss;
    Calcuate_inc_Field();
}

void TFDTD_Hybrid_TFSF::Calcuate_inc_Field()
{
    fopen_s(&pf_debugHybrid, "debug_Hybrid_TFSF_Calcuate_inc_Field.txt", "a");

    const int Nx = (B_Tx_ter - B_Tx_ini + 1);

    // x 위치 기반 distance (x_minus 기준)
    auto dist_int = [&](int i) { return (i + 0.0) * min1D(dx); };
    auto dist_half = [&](int i) { return (i + 0.5) * min1D(dx); };

    // yz planes는 x 고정이므로 상수 distance
    const double dist_yzN_int = 0.0;                 // x = x_minus
    const double dist_yzN_half = 0.0;                 // (필요시)
    const double dist_yzP_int = (Nx - 1.0) * min1D(dx);     // x = x_plus (integer 기준)
    const double dist_yzP_half = (Nx - 0.5) * min1D(dx);     // half 기준 (필요시)


    // ============================
    // (2) XY planes (z- / z+): Ex,Ey and Hx,Hy
    // ============================
    for (int ii = B_Tx_ini; ii <= B_Tx_ter; ii++)
    {
        for (int jj = B_Ty_ini; jj <= B_Ty_ter; jj++)
        {
            const int i = ii - B_Tx_ini;
            const int j = jj - B_Ty_ini;

            if (ii <= B_Tx_ter - 1)
            {
                // Ex_inc (i+1/2, j, z- / z+)
                Ex_xyN_correction[i][j] = Get_PulseE(nsteps, Ex_xym[i][j], dist_half(i));
                Ex_xyP_correction[i][j] = Get_PulseE(nsteps, Ex_xyp[i][j], dist_half(i));
            }

            if (jj <= B_Ty_ter - 1)
            {
                // Ey_inc (i, j+1/2, z- / z+)
                Ey_xyN_correction[i][j] = Get_PulseE(nsteps, Ey_xym[i][j], dist_int(i));
                Ey_xyP_correction[i][j] = Get_PulseE(nsteps, Ey_xyp[i][j], dist_int(i));
            }
        }
    }

    // ============================
    // (3) YZ planes (x- / x+): Ey,Ez and Hy,Hz  (E)
    // ============================
    for (int jj = B_Ty_ini; jj <= B_Ty_ter; jj++)
    {
        for (int kk = B_Tz_ini; kk <= B_Tz_ter; kk++)
        {
            const int j = jj - B_Ty_ini;
            const int k = kk - B_Tz_ini;

            if (jj <= B_Ty_ter - 1)
            {
                // Ey_inc (x-, j+1/2, k) / (x+, j+1/2, k)
                Ey_yzN_correction[j][k] = Get_PulseE(nsteps, Ey_myz[j][k], dist_yzN_int);
                Ey_yzP_correction[j][k] = Get_PulseE(nsteps, Ey_pyz[j][k], dist_yzP_int);
            }

            if (kk <= B_Tz_ter - 1)
            {
                // Ez_inc (x-, j, k+1/2) / (x+, j, k+1/2)
                Ez_yzN_correction[j][k] = Get_PulseE(nsteps, Ez_myz[j][k], dist_yzN_int);
                Ez_yzP_correction[j][k] = Get_PulseE(nsteps, Ez_pyz[j][k], dist_yzP_int);
            }
        }
    }

    // ============================
    // (4) XZ planes (y- / y+): Ex,Ez  (E)
    // ============================
    for (int ii = B_Tx_ini; ii <= B_Tx_ter; ii++)
    {
        for (int kk = B_Tz_ini; kk <= B_Tz_ter; kk++)
        {
            const int i = ii - B_Tx_ini;
            const int k = kk - B_Tz_ini;

            if (ii <= B_Tx_ter - 1)
            {
                // Ex_inc (i+1/2, y-, k) / (i+1/2, y+, k)
                Ex_xzN_correction[i][k] = Get_PulseE(nsteps, Ex_xmz[i][k], dist_half(i));
                Ex_xzP_correction[i][k] = Get_PulseE(nsteps, Ex_xpz[i][k], dist_half(i));
            }

            if (kk <= B_Tz_ter - 1)
            {
                // Ez_inc (i, y-, k+1/2) / (i, y+, k+1/2)
                Ez_xzN_correction[i][k] = Get_PulseE(nsteps, Ez_xmz[i][k], dist_int(i));
                Ez_xzP_correction[i][k] = Get_PulseE(nsteps, Ez_xpz[i][k], dist_int(i));
            }
        }
    }

    // ============================
    // (5) YZ planes H (x- / x+): Hz,Hy
    // ============================
    for (int jj = B_Ty_ini; jj <= B_Ty_ter; jj++)
    {
        for (int kk = B_Tz_ini; kk <= B_Tz_ter; kk++)
        {
            const int j = jj - B_Ty_ini;
            const int k = kk - B_Tz_ini;

            if (jj <= B_Ty_ter - 1)
            {
                // Hz_inc (x-, j+1/2, k) / (x+, j+1/2, k)
                Hz_yzN_correction[j][k] = Get_PulseH(nsteps, Hz_myz[j][k], dist_yzN_half);
                Hz_yzP_correction[j][k] = Get_PulseH(nsteps, Hz_pyz[j][k], dist_yzP_half);
            }

            if (kk <= B_Tz_ter - 1)
            {
                // Hy_inc (x-, j, k+1/2) / (x+, j, k+1/2)
                Hy_yzN_correction[j][k] = Get_PulseH(nsteps, Hy_myz[j][k], dist_yzN_half);
                Hy_yzP_correction[j][k] = Get_PulseH(nsteps, Hy_pyz[j][k], dist_yzP_half);
            }
        }
    }

    // ============================
    // (6) XZ planes H (y- / y+): Hx,Hz
    // ============================
    for (int ii = B_Tx_ini; ii <= B_Tx_ter; ii++)
    {
        for (int kk = B_Tz_ini; kk <= B_Tz_ter; kk++)
        {
            const int i = ii - B_Tx_ini;
            const int k = kk - B_Tz_ini;

            if (kk <= B_Tz_ter - 1)
            {
                // Hx_inc (i, y-, k+1/2) / (i, y+, k+1/2)
                Hx_xzN_correction[i][k] = Get_PulseH(nsteps, Hx_xmz[i][k], dist_int(i));
                Hx_xzP_correction[i][k] = Get_PulseH(nsteps, Hx_xpz[i][k], dist_int(i));
            }

            if (ii <= B_Tx_ter - 1)
            {
                // Hz_inc (i+1/2, y-, k) / (i+1/2, y+, k)
                Hz_xzN_correction[i][k] = Get_PulseH(nsteps, Hz_xmz[i][k], dist_half(i));
                Hz_xzP_correction[i][k] = Get_PulseH(nsteps, Hz_xpz[i][k], dist_half(i));
            }
        }
    }

    // ============================
    // (7) XY planes H (z- / z+): Hx,Hy
    // ============================
    for (int ii = B_Tx_ini; ii <= B_Tx_ter; ii++)
    {
        for (int jj = B_Ty_ini; jj <= B_Ty_ter; jj++)
        {
            const int i = ii - B_Tx_ini;
            const int j = jj - B_Ty_ini;

            if (jj <= B_Ty_ter - 1)
            {
                // Hx_inc (i, j+1/2, z-) / (i, j+1/2, z+)
                Hx_xyN_correction[i][j] = Get_PulseH(nsteps, Hx_xym[i][j], dist_int(i));
                Hx_xyP_correction[i][j] = Get_PulseH(nsteps, Hx_xyp[i][j], dist_int(i));
            }

            if (ii <= B_Tx_ter - 1)
            {
                // Hy_inc (i+1/2, j, z-) / (i+1/2, j, z+)
                Hy_xyN_correction[i][j] = Get_PulseH(nsteps, Hy_xym[i][j], dist_half(i));
                Hy_xyP_correction[i][j] = Get_PulseH(nsteps, Hy_xyp[i][j], dist_half(i));
            }
        }
    }

    // ============================
    // Debug 출력 (유지)
    // ============================
    //double val = Ey_xyP_correction[(B_Tx_ter - B_Tx_ini) / 2][(B_Ty_ter - B_Ty_ini) / 2];
    //fprintf(pf_debugHybrid, "%.6e \n", val);

    fclose(pf_debugHybrid);
}



int TFDTD_Hybrid_TFSF::round(double nominator, double denominator)
{
    int temp;

    if (nominator >= 0)
    {
        temp = (int)(nominator / denominator + 0.5);
    }
    else if (nominator < 0)
    {
        temp = (int)(nominator / denominator - 0.5);
    }
    return temp;
}


void TFDTD_Hybrid_TFSF::Calculate_start_pointE()
{
    double PXminE = 1e100, PYminE = 1e100, PZminE = 1e100;
    // TFSF 박스의 시작점 (최솟값) 을 구하기 위한 함수
    Calculate_start_point(pf_Ey_myz, &PXminE, &PYminE, &PZminE);
    cout << "PXminE = " << PXminE << " PYminE = " << PYminE << " PZminE = " << PZminE << endl;
    Calculate_start_point(pf_Ez_myz, &PXminE, &PYminE, &PZminE);
    cout << "PXminE = " << PXminE << " PYminE = " << PYminE << " PZminE = " << PZminE << endl;
    Calculate_start_point(pf_Ex_xmz, &PXminE, &PYminE, &PZminE);
    cout << "PXminE = " << PXminE << " PYminE = " << PYminE << " PZminE = " << PZminE << endl;
    Calculate_start_point(pf_Ez_xmz, &PXminE, &PYminE, &PZminE);
    cout << "PXminE = " << PXminE << " PYminE = " << PYminE << " PZminE = " << PZminE << endl;
    Calculate_start_point(pf_Ey_xym, &PXminE, &PYminE, &PZminE);
    cout << "PXminE = " << PXminE << " PYminE = " << PYminE << " PZminE = " << PZminE << endl;
    Calculate_start_point(pf_Ex_xym, &PXminE, &PYminE, &PZminE);
    cout << "PXminE = " << PXminE << " PYminE = " << PYminE << " PZminE = " << PZminE << endl;

    PXoffsetE = FindIndex(PXminE, 'X');
    PYoffsetE = FindIndex(PYminE, 'Y');
    PZoffsetE = FindIndex(PZminE, 'Z');

    printf("PXoffsetE = %d PYoffsetE = %d PZoffsetE = %d\n", PXoffsetE, PYoffsetE, PZoffsetE);
}


void TFDTD_Hybrid_TFSF::Calculate_end_pointE()
{
    double PXmaxE = -1e100, PYmaxE = -1e100, PZmaxE = -1e100;
    // TFSF 박스의 끝점 (최댓값) 을 구하기 위한 함수
    Calculate_end_point(pf_Ey_myz, &PXmaxE, &PYmaxE, &PZmaxE);
    cout << "PXmaxE = " << PXmaxE << " PYmaxE = " << PYmaxE << " PZmaxE = " << PZmaxE << endl;
    Calculate_end_point(pf_Ez_myz, &PXmaxE, &PYmaxE, &PZmaxE);
    cout << "PXmaxE = " << PXmaxE << " PYmaxE = " << PYmaxE << " PZmaxE = " << PZmaxE << endl;
    Calculate_end_point(pf_Ex_xmz, &PXmaxE, &PYmaxE, &PZmaxE);
    cout << "PXmaxE = " << PXmaxE << " PYmaxE = " << PYmaxE << " PZmaxE = " << PZmaxE << endl;
    Calculate_end_point(pf_Ez_xmz, &PXmaxE, &PYmaxE, &PZmaxE);
    cout << "PXmaxE = " << PXmaxE << " PYmaxE = " << PYmaxE << " PZmaxE = " << PZmaxE << endl;
    Calculate_end_point(pf_Ey_xym, &PXmaxE, &PYmaxE, &PZmaxE);
    cout << "PXmaxE = " << PXmaxE << " PYmaxE = " << PYmaxE << " PZmaxE = " << PZmaxE << endl;
    Calculate_end_point(pf_Ex_xym, &PXmaxE, &PYmaxE, &PZmaxE);
    cout << "PXmaxE = " << PXmaxE << " PYmaxE = " << PYmaxE << " PZmaxE = " << PZmaxE << endl;

    PXoffsetE2 = FindIndex(PXmaxE, 'X');
    PYoffsetE2 = FindIndex(PYmaxE, 'Y');
    PZoffsetE2 = FindIndex(PZmaxE, 'Z');

	printf("PXoffsetE2 = %d PYoffsetE2 = %d PZoffsetE2 = %d\n", PXoffsetE2, PYoffsetE2, PZoffsetE2);
}



void TFDTD_Hybrid_TFSF::Calculate_start_pointH()
{
    PXminH = 1e100, PYminH = 1e100, PZminH = 1e100;
    Calculate_start_point(pf_Hy_myz, &PXminH, &PYminH, &PZminH);
    Calculate_start_point(pf_Hz_myz, &PXminH, &PYminH, &PZminH);
    Calculate_start_point(pf_Hx_xmz, &PXminH, &PYminH, &PZminH);
    Calculate_start_point(pf_Hz_xmz, &PXminH, &PYminH, &PZminH);
    Calculate_start_point(pf_Hy_xym, &PXminH, &PYminH, &PZminH);
    Calculate_start_point(pf_Hx_xym, &PXminH, &PYminH, &PZminH);

    
	PXoffsetH1 = FindIndex(PXminH, 'X');
	PYoffsetH1 = FindIndex(PYminH, 'Y');
	PZoffsetH1 = FindIndex(PZminH, 'Z');

}

void TFDTD_Hybrid_TFSF::Calculate_end_pointH()
{
    PXminH = -1e100, PYminH = -1e100, PZminH = -1e100;
    Calculate_end_point(pf_Hy_myz, &PXminH, &PYminH, &PZminH);
    Calculate_end_point(pf_Hz_myz, &PXminH, &PYminH, &PZminH);
    Calculate_end_point(pf_Hx_xmz, &PXminH, &PYminH, &PZminH);
    Calculate_end_point(pf_Hz_xmz, &PXminH, &PYminH, &PZminH);
    Calculate_end_point(pf_Hy_xym, &PXminH, &PYminH, &PZminH);
    Calculate_end_point(pf_Hx_xym, &PXminH, &PYminH, &PZminH);


    PXoffsetH2 = FindIndex(PXminH, 'X');
    PYoffsetH2 = FindIndex(PYminH, 'Y');
    PZoffsetH2 = FindIndex(PZminH, 'Z');

}



//void TFDTD_Hybrid_TFSF::Calculate_2D_size(struct_Hybrid_TFSF_File psHTFSFFile, struct_Hybrid_TFSF*** psHTFSF, string HTFSFname)
//{
//    double PX, PY, PZ;
//    double RR, II;
//    double PXmin = 1e100, PYmin = 1e100, PZmin = 1e100;
//    double PXmax = -1e100, PYmax = -1e100, PZmax = -1e100;
//
//    cout << HTFSFname << endl;
//
//    Calculate_start_point(psHTFSFFile, &PXmin, &PYmin, &PZmin);
//    Calculate_end_point(psHTFSFFile, &PXmax, &PYmax, &PZmax);
//
//	int xstart, ystart, zstart;
//	int xend, yend, zend;
//
//	xstart = FindIndex(PXmin, 'X');
//	ystart = FindIndex(PYmin, 'Y');
//	zstart = FindIndex(PZmin, 'Z');
//	xend = FindIndex(PXmax, 'X');
//	yend = FindIndex(PYmax, 'Y');
//	zend = FindIndex(PZmax, 'Z');
//
//	int xlen = xend - xstart + 1;
//    int ylen = yend - ystart + 1;
//    int zlen = zend - zstart + 1;
//    int ii = 0, jj = 0, kk = 0;
//
//	printf("xlen=%d ylen=%d zlen=%d\n", xlen, ylen, zlen);
//
//    if (xlen == 0)
//    {
//        *psHTFSF = Tptr2< struct_Hybrid_TFSF>(ylen + 1, zlen + 1);
//        for (int ll = 0; ll < psHTFSFFile.Datalen; ll++)
//        {
//            fscanf(psHTFSFFile.fp, "  %le  %le  %le  %le  %le\n", &PX, &PY, &PZ, &RR, &II);
//            int pY = FindIndex(PY, 'Y');
//			int pZ = FindIndex(PZ, 'Z');
//			jj = pY - ystart;
//            kk = pZ - zstart;
//            //cout << PY << " " << PYmin << " " <<  jj << "  " << PZ << " " << PZmin << " " << kk << endl;
//            (*(*psHTFSF + jj) + kk)->Abs = get_abs(RR, II);
//            (*(*psHTFSF + jj) + kk)->Phase = get_Phase(RR, II);
//            /// xplane, yini, yter, zini, zter;
//
//
//
//        }
//        fseek(psHTFSFFile.fp, 0L, SEEK_SET);
//    }
//    else if (ylen == 0)
//    {
//        *psHTFSF = Tptr2< struct_Hybrid_TFSF>(xlen + 1, zlen + 1);
//        for (int ll = 0; ll < psHTFSFFile.Datalen; ll++)
//        {
//            fscanf(psHTFSFFile.fp, "  %le  %le  %le  %le  %le\n", &PX, &PY, &PZ, &RR, &II);
//            int pX = FindIndex(PX, 'X');
//            int pZ = FindIndex(PZ, 'Z');
//            ii = pX - xstart;
//            kk = pZ - zstart;
//            (*(*psHTFSF + ii) + kk)->Abs = get_abs(RR, II);
//            (*(*psHTFSF + ii) + kk)->Phase = get_Phase(RR, II);
//        }
//        fseek(psHTFSFFile.fp, 0L, SEEK_SET);
//    }
//    else if (zlen == 0)
//    {
//        *psHTFSF = Tptr2< struct_Hybrid_TFSF>(xlen + 1, ylen + 1); //여기에서 차원 변환이 일어남. 처음 생성할 때는 3차원 포인터 벡터였음.
//        for (int ll = 0; ll < psHTFSFFile.Datalen; ll++)
//        {
//            fscanf(psHTFSFFile.fp, "  %le  %le  %le  %le  %le\n", &PX, &PY, &PZ, &RR, &II);
//            int pX = FindIndex(PX, 'X');
//            int pY = FindIndex(PY, 'Y');
//            ii = pX - xstart;
//            jj = pY - ystart;
//            (*(*psHTFSF + ii) + jj)->Abs = get_abs(RR, II);
//            (*(*psHTFSF + ii) + jj)->Phase = get_Phase(RR, II);
//        }
//        fseek(psHTFSFFile.fp, 0L, SEEK_SET);
//    }
//}


void TFDTD_Hybrid_TFSF::Calculate_start_point(struct_Hybrid_TFSF_File psHTFSFFile, double* PXmin, double* PYmin, double* PZmin)
{
    double PX, PY, PZ;
    double RR, II;
    for (int ii = 0; ii < psHTFSFFile.Datalen; ii++)
    {
        fscanf(psHTFSFFile.fp, "  %le  %le  %le  %le  %le\n", &PX, &PY, &PZ, &RR, &II);
        if (PX < *PXmin)
            *PXmin = PX;
        if (PY < *PYmin)
            *PYmin = PY;
        if (PZ < *PZmin)
            *PZmin = PZ;
    }
    fseek(psHTFSFFile.fp, 0L, SEEK_SET);
}



void TFDTD_Hybrid_TFSF::Calculate_end_point(struct_Hybrid_TFSF_File psHTFSFFile, double* PXmax, double* PYmax, double* PZmax)
{
    double PX, PY, PZ;
    double RR, II;
    for (int ii = 0; ii < psHTFSFFile.Datalen; ii++)
    {
        fscanf(psHTFSFFile.fp, "  %le  %le  %le  %le  %le\n", &PX, &PY, &PZ, &RR, &II);
        if (PX > *PXmax)
            *PXmax = PX;
        if (PY > *PYmax)
            *PYmax = PY;
        if (PZ > *PZmax)
            *PZmax = PZ;
    }
    fseek(psHTFSFFile.fp, 0L, SEEK_SET);
}



double  TFDTD_Hybrid_TFSF::get_abs(double re, double im)
{
    double abs;
    abs = sqrt(re * re + im * im);
    return (double)abs;
}

double  TFDTD_Hybrid_TFSF::get_Phase(double x, double y)
{
    double temp = 1e-30;
    double phase = atan2(y, x + temp);
    //    if (x > 0)
    //    {
    //        phase = 2 * atan(y / (sqrt(x * x + y * y) + x));
    //    }
    //    else if (x <= 0 && y != 0)
    //    {
    //        phase = 2 * atan((sqrt(x * x + y * y) - x) / y);
    //    }
    //    else if (x < 0 && y == 0)
    //    {
    //        phase = m_pi;
    //    }
    //    else if (x == 0 && y == 0)
    //    {
    //        phase = 0;
    //    }
    return (double)phase;
}

int TFDTD_Hybrid_TFSF::count_line(FILE* pfCount)
{
    int length = 0;
    char temp[256];
    while (!feof(pfCount))
    {
        fgets(temp, 256, pfCount);
        length++; // 줄 수 증가
    }
    fseek(pfCount, 0L, SEEK_SET);
    return length;
}


double TFDTD_Hybrid_TFSF::Get_PulseE(double ss, struct_Hybrid_TFSF psHTFSF, double distance) //psHTFSF[ii][jj].Abs, psHTFSF[ii][jj].Phase
{
    double t = ss * dt;
    double tau = t - distance / V_tilt;
    double window = sin(m_pi * f0 * tau / (2 * alpha)) * sin(m_pi * f0 * tau / (2 * alpha));
    if (tau < 0)
        Pulse = 0;
    else if (tau < alpha * T_0)
        Pulse = psHTFSF.Abs * window * sin(2. * m_pi * f0 * t + psHTFSF.Phase);
    else
        Pulse = psHTFSF.Abs * sin(2. * m_pi * f0 * t + psHTFSF.Phase);

    return Pulse; // time delay 부분은 여기 봐야됨. 
}



double TFDTD_Hybrid_TFSF::Get_PulseH(double ss, struct_Hybrid_TFSF psHTFSF, double distance)
{
    double t = (ss + 0.5) * dt;            // global time
    double tau = t - distance / V_tilt;                 // 이 셀에서의 local time
    double window = sin(m_pi * f0 * tau / (2 * alpha)) * sin(m_pi * f0 * tau / (2 * alpha));
    if (tau < 0)
        Pulse = 0;
    else if (tau < alpha * T_0)
    {
        //temp = psHTFSF.Abs * 0.5 * (1. - cos((m_pi * f0 * dt * ss) / alpha)) * sin(2. * m_pi * f0 * dt * (ss + 0) + psHTFSF.Phase - omega * time);
        Pulse = psHTFSF.Abs * window * sin(2.0 * m_pi * f0 * t + psHTFSF.Phase);
    }
    else
    {
        //temp = psHTFSF.Abs * sin(2. * m_pi * f0 * dt * (ss + 0) + psHTFSF.Phase - omega * time);
        Pulse = psHTFSF.Abs * sin(2.0 * m_pi * f0 * t + psHTFSF.Phase);
    }
    return Pulse;
}


void TFDTD_Hybrid_TFSF::Calculate_2D_size_EyYZ(struct_Hybrid_TFSF_File psHTFSFFile, struct_Hybrid_TFSF*** psHTFSF, string HTFSFname)
{
    double PX, PY, PZ;
    double RR, II;
    double PXmin = 1e100, PYmin = 1e100, PZmin = 1e100;
    double PXmax = -1e100, PYmax = -1e100, PZmax = -1e100;

    cout << HTFSFname << endl;

    Calculate_start_point(psHTFSFFile, &PXmin, &PYmin, &PZmin);
    Calculate_end_point(psHTFSFFile, &PXmax, &PYmax, &PZmax);

    int xstart, ystart, zstart;
    int xend, yend, zend;

    
    ystart = findStaggeredIndex(PYmin, 'Y');
    zstart = findIndex(PZmin, 'Z');
    
    yend = findStaggeredIndex(PYmax, 'Y');
    zend = findIndex(PZmax, 'Z');

    
    int ylen = yend - ystart + 1;
    int zlen = zend - zstart + 1;
    int ii = 0, jj = 0, kk = 0;

    printf("ylen=%d zlen=%d\n", ylen, zlen);

    *psHTFSF = Tptr2< struct_Hybrid_TFSF>(ylen + 1, zlen + 1);
    for (int ll = 0; ll < psHTFSFFile.Datalen; ll++)
    {
        fscanf(psHTFSFFile.fp, "  %le  %le  %le  %le  %le\n", &PX, &PY, &PZ, &RR, &II);
        int pY = findStaggeredIndex(PY, 'Y');
        int pZ = findIndex(PZ, 'Z');
        jj = pY - ystart;
        kk = pZ - zstart;
        //cout << PY << " " << PYmin << " " <<  jj << "  " << PZ << " " << PZmin << " " << kk << endl;
        (*(*psHTFSF + jj) + kk)->Abs = get_abs(RR, II);
        (*(*psHTFSF + jj) + kk)->Phase = get_Phase(RR, II);
        /// xplane, yini, yter, zini, zter;



    }
    fseek(psHTFSFFile.fp, 0L, SEEK_SET); 
}

void TFDTD_Hybrid_TFSF::Calculate_2D_size_EzYZ(struct_Hybrid_TFSF_File psHTFSFFile, struct_Hybrid_TFSF*** psHTFSF, string HTFSFname)
{
    double PX, PY, PZ;
    double RR, II;
    double PXmin = 1e100, PYmin = 1e100, PZmin = 1e100;
    double PXmax = -1e100, PYmax = -1e100, PZmax = -1e100;

    cout << HTFSFname << endl;

    Calculate_start_point(psHTFSFFile, &PXmin, &PYmin, &PZmin);
    Calculate_end_point(psHTFSFFile, &PXmax, &PYmax, &PZmax);

    int xstart, ystart, zstart;
    int xend, yend, zend;


    ystart = findIndex(PYmin, 'Y');
    zstart = findStaggeredIndex(PZmin, 'Z');

    yend = findIndex(PYmax, 'Y');
    zend = findStaggeredIndex(PZmax, 'Z');


    int ylen = yend - ystart + 1;
    int zlen = zend - zstart + 1;
    int ii = 0, jj = 0, kk = 0;

    printf("ylen=%d zlen=%d\n", ylen, zlen);

    *psHTFSF = Tptr2< struct_Hybrid_TFSF>(ylen + 1, zlen + 1);
    for (int ll = 0; ll < psHTFSFFile.Datalen; ll++)
    {
        fscanf(psHTFSFFile.fp, "  %le  %le  %le  %le  %le\n", &PX, &PY, &PZ, &RR, &II);
        int pY = findIndex(PY, 'Y');
        int pZ = findStaggeredIndex(PZ, 'Z');
        jj = pY - ystart;
        kk = pZ - zstart;
        //cout << PY << " " << PYmin << " " <<  jj << "  " << PZ << " " << PZmin << " " << kk << endl;
        (*(*psHTFSF + jj) + kk)->Abs = get_abs(RR, II);
        (*(*psHTFSF + jj) + kk)->Phase = get_Phase(RR, II);
        /// xplane, yini, yter, zini, zter;



    }
    fseek(psHTFSFFile.fp, 0L, SEEK_SET);
}

void TFDTD_Hybrid_TFSF::Calculate_2D_size_HyYZ(struct_Hybrid_TFSF_File psHTFSFFile, struct_Hybrid_TFSF*** psHTFSF, string HTFSFname)
{
    double PX, PY, PZ;
    double RR, II;
    double PXmin = 1e100, PYmin = 1e100, PZmin = 1e100;
    double PXmax = -1e100, PYmax = -1e100, PZmax = -1e100;

    cout << HTFSFname << endl;

    Calculate_start_point(psHTFSFFile, &PXmin, &PYmin, &PZmin);
    Calculate_end_point(psHTFSFFile, &PXmax, &PYmax, &PZmax);

    int xstart, ystart, zstart;
    int xend, yend, zend;


    ystart = findStaggeredIndex(PYmin, 'Y');
    zstart = findIndex(PZmin, 'Z');

    yend = findStaggeredIndex(PYmax, 'Y');
    zend = findIndex(PZmax, 'Z');


    int ylen = yend - ystart + 1;
    int zlen = zend - zstart + 1;
    int ii = 0, jj = 0, kk = 0;

    printf("ylen=%d zlen=%d\n", ylen, zlen);

    *psHTFSF = Tptr2< struct_Hybrid_TFSF>(ylen + 1, zlen + 1);
    for (int ll = 0; ll < psHTFSFFile.Datalen; ll++)
    {
        fscanf(psHTFSFFile.fp, "  %le  %le  %le  %le  %le\n", &PX, &PY, &PZ, &RR, &II);
        int pY = findStaggeredIndex(PY, 'Y');
        int pZ = findIndex(PZ, 'Z');
        jj = pY - ystart;
        kk = pZ - zstart;
        //cout << PY << " " << PYmin << " " <<  jj << "  " << PZ << " " << PZmin << " " << kk << endl;
        (*(*psHTFSF + jj) + kk)->Abs = get_abs(RR, II);
        (*(*psHTFSF + jj) + kk)->Phase = get_Phase(RR, II);
        /// xplane, yini, yter, zini, zter;



    }
    fseek(psHTFSFFile.fp, 0L, SEEK_SET);
}

void TFDTD_Hybrid_TFSF::Calculate_2D_size_HzYZ(struct_Hybrid_TFSF_File psHTFSFFile, struct_Hybrid_TFSF*** psHTFSF, string HTFSFname)
{
    double PX, PY, PZ;
    double RR, II;
    double PXmin = 1e100, PYmin = 1e100, PZmin = 1e100;
    double PXmax = -1e100, PYmax = -1e100, PZmax = -1e100;

    cout << HTFSFname << endl;

    Calculate_start_point(psHTFSFFile, &PXmin, &PYmin, &PZmin);
    Calculate_end_point(psHTFSFFile, &PXmax, &PYmax, &PZmax);

    int xstart, ystart, zstart;
    int xend, yend, zend;


    ystart = findIndex(PYmin, 'Y');
    zstart = findStaggeredIndex(PZmin, 'Z');

    yend = findIndex(PYmax, 'Y');
    zend = findStaggeredIndex(PZmax, 'Z');


    int ylen = yend - ystart + 1;
    int zlen = zend - zstart + 1;
    int ii = 0, jj = 0, kk = 0;

    printf("ylen=%d zlen=%d\n", ylen, zlen);

    *psHTFSF = Tptr2< struct_Hybrid_TFSF>(ylen + 1, zlen + 1);
    for (int ll = 0; ll < psHTFSFFile.Datalen; ll++)
    {
        fscanf(psHTFSFFile.fp, "  %le  %le  %le  %le  %le\n", &PX, &PY, &PZ, &RR, &II);
        int pY = findIndex(PY, 'Y');
        int pZ = findStaggeredIndex(PZ, 'Z');
        jj = pY - ystart;
        kk = pZ - zstart;
        //cout << PY << " " << PYmin << " " <<  jj << "  " << PZ << " " << PZmin << " " << kk << endl;
        (*(*psHTFSF + jj) + kk)->Abs = get_abs(RR, II);
        (*(*psHTFSF + jj) + kk)->Phase = get_Phase(RR, II);
        /// xplane, yini, yter, zini, zter;



    }
    fseek(psHTFSFFile.fp, 0L, SEEK_SET);
}

void TFDTD_Hybrid_TFSF::Calculate_2D_size_ExXZ(struct_Hybrid_TFSF_File psHTFSFFile, struct_Hybrid_TFSF*** psHTFSF, string HTFSFname)
{
    double PX, PY, PZ;
    double RR, II;
    double PXmin = 1e100, PYmin = 1e100, PZmin = 1e100;
    double PXmax = -1e100, PYmax = -1e100, PZmax = -1e100;

    cout << HTFSFname << endl;

    Calculate_start_point(psHTFSFFile, &PXmin, &PYmin, &PZmin);
    Calculate_end_point(psHTFSFFile, &PXmax, &PYmax, &PZmax);

    int xstart, ystart, zstart;
    int xend, yend, zend;

    xstart = findStaggeredIndex(PXmin, 'X');
    zstart = findIndex(PZmin, 'Z');
    xend = findStaggeredIndex(PXmax, 'X');
    zend = findIndex(PZmax, 'Z');

    int xlen = xend - xstart + 1;
    int zlen = zend - zstart + 1;
    int ii = 0, jj = 0, kk = 0;

    printf("xlen=%d zlen=%d\n", xlen, zlen);

    
    
    
    *psHTFSF = Tptr2< struct_Hybrid_TFSF>(xlen + 1, zlen + 1);
    for (int ll = 0; ll < psHTFSFFile.Datalen; ll++)
    {
        fscanf(psHTFSFFile.fp, "  %le  %le  %le  %le  %le\n", &PX, &PY, &PZ, &RR, &II);
        int pX = findStaggeredIndex(PX, 'X');
        int pZ = findIndex(PZ, 'Z');
        ii = pX - xstart;
        kk = pZ - zstart;
        (*(*psHTFSF + ii) + kk)->Abs = get_abs(RR, II);
        (*(*psHTFSF + ii) + kk)->Phase = get_Phase(RR, II);
    }
    fseek(psHTFSFFile.fp, 0L, SEEK_SET);
     
}

void TFDTD_Hybrid_TFSF::Calculate_2D_size_EzXZ(struct_Hybrid_TFSF_File psHTFSFFile, struct_Hybrid_TFSF*** psHTFSF, string HTFSFname)
{
    double PX, PY, PZ;
    double RR, II;
    double PXmin = 1e100, PYmin = 1e100, PZmin = 1e100;
    double PXmax = -1e100, PYmax = -1e100, PZmax = -1e100;

    cout << HTFSFname << endl;

    Calculate_start_point(psHTFSFFile, &PXmin, &PYmin, &PZmin);
    Calculate_end_point(psHTFSFFile, &PXmax, &PYmax, &PZmax);

    int xstart, ystart, zstart;
    int xend, yend, zend;

    xstart = findIndex(PXmin, 'X');
    zstart = findStaggeredIndex(PZmin, 'Z');
    xend = findIndex(PXmax, 'X');
    zend = findStaggeredIndex(PZmax, 'Z');

    int xlen = xend - xstart + 1;
    int zlen = zend - zstart + 1;
    int ii = 0, jj = 0, kk = 0;

    printf("xlen=%d zlen=%d\n", xlen, zlen);




    *psHTFSF = Tptr2< struct_Hybrid_TFSF>(xlen + 1, zlen + 1);
    for (int ll = 0; ll < psHTFSFFile.Datalen; ll++)
    {
        fscanf(psHTFSFFile.fp, "  %le  %le  %le  %le  %le\n", &PX, &PY, &PZ, &RR, &II);
        int pX = findIndex(PX, 'X');
        int pZ = findStaggeredIndex(PZ, 'Z');
        ii = pX - xstart;
        kk = pZ - zstart;
        (*(*psHTFSF + ii) + kk)->Abs = get_abs(RR, II);
        (*(*psHTFSF + ii) + kk)->Phase = get_Phase(RR, II);
    }
    fseek(psHTFSFFile.fp, 0L, SEEK_SET);

}

void TFDTD_Hybrid_TFSF::Calculate_2D_size_HxXZ(struct_Hybrid_TFSF_File psHTFSFFile, struct_Hybrid_TFSF*** psHTFSF, string HTFSFname)
{
    double PX, PY, PZ;
    double RR, II;
    double PXmin = 1e100, PYmin = 1e100, PZmin = 1e100;
    double PXmax = -1e100, PYmax = -1e100, PZmax = -1e100;

    cout << HTFSFname << endl;

    Calculate_start_point(psHTFSFFile, &PXmin, &PYmin, &PZmin);
    Calculate_end_point(psHTFSFFile, &PXmax, &PYmax, &PZmax);

    int xstart, ystart, zstart;
    int xend, yend, zend;

    xstart = findStaggeredIndex(PXmin, 'X');
    zstart = findIndex(PZmin, 'Z');
    xend = findStaggeredIndex(PXmax, 'X');
    zend = findIndex(PZmax, 'Z');

    int xlen = xend - xstart + 1;
    int zlen = zend - zstart + 1;
    int ii = 0, jj = 0, kk = 0;

    printf("xlen=%d zlen=%d\n", xlen, zlen);




    *psHTFSF = Tptr2< struct_Hybrid_TFSF>(xlen + 1, zlen + 1);
    for (int ll = 0; ll < psHTFSFFile.Datalen; ll++)
    {
        fscanf(psHTFSFFile.fp, "  %le  %le  %le  %le  %le\n", &PX, &PY, &PZ, &RR, &II);
        int pX = findStaggeredIndex(PX, 'X');
        int pZ = findIndex(PZ, 'Z');
        ii = pX - xstart;
        kk = pZ - zstart;
        (*(*psHTFSF + ii) + kk)->Abs = get_abs(RR, II);
        (*(*psHTFSF + ii) + kk)->Phase = get_Phase(RR, II);
    }
    fseek(psHTFSFFile.fp, 0L, SEEK_SET);

}

void TFDTD_Hybrid_TFSF::Calculate_2D_size_HzXZ(struct_Hybrid_TFSF_File psHTFSFFile, struct_Hybrid_TFSF*** psHTFSF, string HTFSFname)
{
    double PX, PY, PZ;
    double RR, II;
    double PXmin = 1e100, PYmin = 1e100, PZmin = 1e100;
    double PXmax = -1e100, PYmax = -1e100, PZmax = -1e100;

    cout << HTFSFname << endl;

    Calculate_start_point(psHTFSFFile, &PXmin, &PYmin, &PZmin);
    Calculate_end_point(psHTFSFFile, &PXmax, &PYmax, &PZmax);

    int xstart, ystart, zstart;
    int xend, yend, zend;

    xstart = findIndex(PXmin, 'X');
    zstart = findStaggeredIndex(PZmin, 'Z');
    xend = findIndex(PXmax, 'X');
    zend = findStaggeredIndex(PZmax, 'Z');

    int xlen = xend - xstart + 1;
    int zlen = zend - zstart + 1;
    int ii = 0, jj = 0, kk = 0;

    printf("xlen=%d zlen=%d\n", xlen, zlen);




    *psHTFSF = Tptr2< struct_Hybrid_TFSF>(xlen + 1, zlen + 1);
    for (int ll = 0; ll < psHTFSFFile.Datalen; ll++)
    {
        fscanf(psHTFSFFile.fp, "  %le  %le  %le  %le  %le\n", &PX, &PY, &PZ, &RR, &II);
        int pX = findIndex(PX, 'X');
        int pZ = findStaggeredIndex(PZ, 'Z');
        ii = pX - xstart;
        kk = pZ - zstart;
        (*(*psHTFSF + ii) + kk)->Abs = get_abs(RR, II);
        (*(*psHTFSF + ii) + kk)->Phase = get_Phase(RR, II);
    }
    fseek(psHTFSFFile.fp, 0L, SEEK_SET);

}


void TFDTD_Hybrid_TFSF::Calculate_2D_size_ExXY(struct_Hybrid_TFSF_File psHTFSFFile, struct_Hybrid_TFSF*** psHTFSF, string HTFSFname)
{
    double PX, PY, PZ;
    double RR, II;
    double PXmin = 1e100, PYmin = 1e100, PZmin = 1e100;
    double PXmax = -1e100, PYmax = -1e100, PZmax = -1e100;

    cout << HTFSFname << endl;

    Calculate_start_point(psHTFSFFile, &PXmin, &PYmin, &PZmin);
    Calculate_end_point(psHTFSFFile, &PXmax, &PYmax, &PZmax);

    int xstart, ystart, zstart;
    int xend, yend, zend;

    xstart = findStaggeredIndex(PXmin, 'X');
    ystart = findIndex(PYmin, 'Y');
    xend = findStaggeredIndex(PXmax, 'X');
    yend = findIndex(PYmax, 'Y');

    int xlen = xend - xstart + 1;
    int ylen = yend - ystart + 1;
    int ii = 0, jj = 0, kk = 0;

    printf("xlen=%d ylen=%d\n", xlen, ylen);

    
    
    *psHTFSF = Tptr2< struct_Hybrid_TFSF>(xlen + 1, ylen + 1); //여기에서 차원 변환이 일어남. 처음 생성할 때는 3차원 포인터 벡터였음.
    for (int ll = 0; ll < psHTFSFFile.Datalen; ll++)
    {
        fscanf(psHTFSFFile.fp, "  %le  %le  %le  %le  %le\n", &PX, &PY, &PZ, &RR, &II);
        int pX = findStaggeredIndex(PX, 'X');
        int pY = findIndex(PY, 'Y');
        ii = pX - xstart;
        jj = pY - ystart;
        (*(*psHTFSF + ii) + jj)->Abs = get_abs(RR, II);
        (*(*psHTFSF + ii) + jj)->Phase = get_Phase(RR, II);
    }
    fseek(psHTFSFFile.fp, 0L, SEEK_SET);
    
}

void TFDTD_Hybrid_TFSF::Calculate_2D_size_EyXY(struct_Hybrid_TFSF_File psHTFSFFile, struct_Hybrid_TFSF*** psHTFSF, string HTFSFname)
{
    double PX, PY, PZ;
    double RR, II;
    double PXmin = 1e100, PYmin = 1e100, PZmin = 1e100;
    double PXmax = -1e100, PYmax = -1e100, PZmax = -1e100;

    cout << HTFSFname << endl;

    Calculate_start_point(psHTFSFFile, &PXmin, &PYmin, &PZmin);
    Calculate_end_point(psHTFSFFile, &PXmax, &PYmax, &PZmax);

    int xstart, ystart, zstart;
    int xend, yend, zend;

    xstart = findIndex(PXmin, 'X');
    ystart = findStaggeredIndex(PYmin, 'Y');
    xend = findIndex(PXmax, 'X');
    yend = findStaggeredIndex(PYmax, 'Y');

    int xlen = xend - xstart + 1;
    int ylen = yend - ystart + 1;
    int ii = 0, jj = 0, kk = 0;

    printf("xlen=%d ylen=%d\n", xlen, ylen);



    *psHTFSF = Tptr2< struct_Hybrid_TFSF>(xlen + 1, ylen + 1); //여기에서 차원 변환이 일어남. 처음 생성할 때는 3차원 포인터 벡터였음.
    for (int ll = 0; ll < psHTFSFFile.Datalen; ll++)
    {
        fscanf(psHTFSFFile.fp, "  %le  %le  %le  %le  %le\n", &PX, &PY, &PZ, &RR, &II);
        int pX = findIndex(PX, 'X');
        int pY = findStaggeredIndex(PY, 'Y');
        ii = pX - xstart;
        jj = pY - ystart;
        (*(*psHTFSF + ii) + jj)->Abs = get_abs(RR, II);
        (*(*psHTFSF + ii) + jj)->Phase = get_Phase(RR, II);
    }
    fseek(psHTFSFFile.fp, 0L, SEEK_SET);

}

void TFDTD_Hybrid_TFSF::Calculate_2D_size_HxXY(struct_Hybrid_TFSF_File psHTFSFFile, struct_Hybrid_TFSF*** psHTFSF, string HTFSFname)
{
    double PX, PY, PZ;
    double RR, II;
    double PXmin = 1e100, PYmin = 1e100, PZmin = 1e100;
    double PXmax = -1e100, PYmax = -1e100, PZmax = -1e100;

    cout << HTFSFname << endl;

    Calculate_start_point(psHTFSFFile, &PXmin, &PYmin, &PZmin);
    Calculate_end_point(psHTFSFFile, &PXmax, &PYmax, &PZmax);

    int xstart, ystart, zstart;
    int xend, yend, zend;

    xstart = findStaggeredIndex(PXmin, 'X');
    ystart = findIndex(PYmin, 'Y');
    xend = findStaggeredIndex(PXmax, 'X');
    yend = findIndex(PYmax, 'Y');

    int xlen = xend - xstart + 1;
    int ylen = yend - ystart + 1;
    int ii = 0, jj = 0, kk = 0;

    printf("xlen=%d ylen=%d\n", xlen, ylen);



    *psHTFSF = Tptr2< struct_Hybrid_TFSF>(xlen + 1, ylen + 1); //여기에서 차원 변환이 일어남. 처음 생성할 때는 3차원 포인터 벡터였음.
    for (int ll = 0; ll < psHTFSFFile.Datalen; ll++)
    {
        fscanf(psHTFSFFile.fp, "  %le  %le  %le  %le  %le\n", &PX, &PY, &PZ, &RR, &II);
        int pX = findStaggeredIndex(PX, 'X');
        int pY = findIndex(PY, 'Y');
        ii = pX - xstart;
        jj = pY - ystart;
        (*(*psHTFSF + ii) + jj)->Abs = get_abs(RR, II);
        (*(*psHTFSF + ii) + jj)->Phase = get_Phase(RR, II);
    }
    fseek(psHTFSFFile.fp, 0L, SEEK_SET);

}

void TFDTD_Hybrid_TFSF::Calculate_2D_size_HyXY(struct_Hybrid_TFSF_File psHTFSFFile, struct_Hybrid_TFSF*** psHTFSF, string HTFSFname)
{
    double PX, PY, PZ;
    double RR, II;
    double PXmin = 1e100, PYmin = 1e100, PZmin = 1e100;
    double PXmax = -1e100, PYmax = -1e100, PZmax = -1e100;

    cout << HTFSFname << endl;

    Calculate_start_point(psHTFSFFile, &PXmin, &PYmin, &PZmin);
    Calculate_end_point(psHTFSFFile, &PXmax, &PYmax, &PZmax);

    int xstart, ystart, zstart;
    int xend, yend, zend;

    xstart = findIndex(PXmin, 'X');
    ystart = findStaggeredIndex(PYmin, 'Y');
    xend = findIndex(PXmax, 'X');
    yend = findStaggeredIndex(PYmax, 'Y');

    int xlen = xend - xstart + 1;
    int ylen = yend - ystart + 1;
    int ii = 0, jj = 0, kk = 0;

    printf("xlen=%d ylen=%d\n", xlen, ylen);



    *psHTFSF = Tptr2<struct_Hybrid_TFSF>(xlen + 1, ylen + 1); //여기에서 차원 변환이 일어남. 처음 생성할 때는 3차원 포인터 벡터였음.
    for (int ll = 0; ll < psHTFSFFile.Datalen; ll++)
    {
        fscanf(psHTFSFFile.fp, "  %le  %le  %le  %le  %le\n", &PX, &PY, &PZ, &RR, &II);
        int pX = findIndex(PX, 'X');
        int pY = findStaggeredIndex(PY, 'Y');
        ii = pX - xstart;
        jj = pY - ystart;
        (*(*psHTFSF + ii) + jj)->Abs = get_abs(RR, II);
        (*(*psHTFSF + ii) + jj)->Phase = get_Phase(RR, II);
    }
    fseek(psHTFSFFile.fp, 0L, SEEK_SET);

}

int TFDTD_Hybrid_TFSF::findIndex(double val, char coordinate)
{
    int index = -1;
    double eps = 1e-12;
    if (coordinate == 'X')
    {
        for (int ii = 0; ii < txsize + 1; ii++)
        {
            if (fabs(val - xcoordinate[ii]) < eps)
            {
                index = ii;
            }
        }
    }
    else if (coordinate == 'Y')
    {
        for (int jj = 0; jj < tysize + 1; jj++)
        {
            if (fabs(val - ycoordinate[jj]) < eps)
            {
                index = jj;
			}
        }
    }
    else if (coordinate == 'Z')
    {
        for (int kk = 0; kk < tzsize + 1; kk++)
        {
            if (fabs(val - zcoordinate[kk]) < eps)
            {
                index = kk;
            }
        }
	}

    return index;
}

int TFDTD_Hybrid_TFSF::findStaggeredIndex(double val, char coordinate)
{
    int index = -1;
    double eps = 1e-12;
    if (coordinate == 'X')
    {
        for (int ii = 0; ii < txsize; ii++)
        {
            if (fabs(val - xStaggered[ii]) < eps)
            {
                index = ii;
            }
        }
    }
    else if (coordinate == 'Y')
    {
        for (int jj = 0; jj < tysize; jj++)
        {
            if (fabs(val - yStaggered[jj]) < eps)
            {
                index = jj;
            }
        }
    }
    else if (coordinate == 'Z')
    {
        for (int kk = 0; kk < tzsize; kk++)
        {
            if (fabs(val - zStaggered[kk]) < eps)
            {
                index = kk;
            }
        }
    }
    return index;
}