#pragma once
#include "TFDTD_Boundary_TFSF.h"
#include <unordered_map>
#include <complex>
#include <algorithm>
#include <cmath>
#include <array>
#include <vector>

constexpr double PI = 3.141592653589793;
static constexpr double SCALE = 1e6; // 1e-6 m 해상도(원 코드 동일) :contentReference[oaicite:4]{index=4}

using cd = std::complex<double>;


typedef struct Vec3 {
    double x, y, z;
} Vec3;

typedef struct Key3 {
    int ix, iy, iz;
    bool operator==(const Key3& o) const { return ix == o.ix && iy == o.iy && iz == o.iz; }
} Key3;

typedef struct Key3Hash {
    size_t operator()(const Key3& k) const {
        // 간단 해시
        return (size_t)k.ix * 73856093u ^ (size_t)k.iy * 19349663u ^ (size_t)k.iz * 83492791u;
    }
} Key3Hash;

typedef struct Sample3D {
    double x, y, z;
    std::complex<double> f;
} Sample3D;

typedef struct NormalEq {
    double ATA[3][3];
    double ATb[3];
} NormalEq;

typedef struct KeyLL { long long ix, iy, iz; } KeyLL;

typedef struct KeyLLHash {
    size_t operator()(const KeyLL& k) const noexcept {
        size_t h = 1469598103934665603ull;
        auto mix = [&](long long v) {
            h ^= std::hash<long long>()(v) + 0x9e3779b97f4a7c15ull + (h << 6) + (h >> 2);
            };
        mix(k.ix); mix(k.iy); mix(k.iz);
        return h;
    }
} KeyLLHash;

typedef struct KeyLLEq {
    bool operator()(const KeyLL& a, const KeyLL& b) const noexcept {
        return a.ix == b.ix && a.iy == b.iy && a.iz == b.iz;
    }
} KeyLLEq;


inline double wrap_to_pi(double x)
{
    while (x > PI) x -= 2.0 * PI;
    while (x < -PI) x += 2.0 * PI;
    return x;
}

inline double SmoothRamp(double t, double Tr)
{
    if (t <= 0.0) return 0.0;
    if (t >= Tr)  return 1.0;

    double s = sin(0.5 * M_PI * t / Tr);
    return s * s;
}

inline double dot(const Vec3& a, const Vec3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

inline double norm(const Vec3& v) {
    return std::sqrt(dot(v, v));
}

inline cd phasor(double A, double phase) {
    return A * std::exp(cd(0.0, phase));
}



inline Key3 make_key(double x, double y, double z, double tol) {
    return Key3{
        (int)llround(x / tol),
        (int)llround(y / tol),
        (int)llround(z / tol)
    };
}

inline Vec3 normalize(const Vec3& v) {
    double n = norm(v);
    return (n > 0) ? Vec3{ v.x / n, v.y / n, v.z / n } : Vec3{ 0,0,0 };
}

static inline void init_normal(NormalEq& ne) {
    for (int i = 0;i < 3;++i) {
        ne.ATb[i] = 0.0;
        for (int j = 0;j < 3;++j) ne.ATA[i][j] = 0.0;
    }
};

static inline void accumulate(NormalEq& ne, double dx, double dy, double dz, double b) {
    double dr[3] = { dx,dy,dz };
    for (int i = 0;i < 3;++i)
        for (int j = 0;j < 3;++j)
            ne.ATA[i][j] += dr[i] * dr[j];
    for (int i = 0;i < 3;++i)
        ne.ATb[i] += dr[i] * b;
};

static inline KeyLL make_keyLL(double x, double y, double z) {
    return KeyLL{ (long long)llround(x * SCALE), (long long)llround(y * SCALE), (long long)llround(z * SCALE) };
}


typedef struct struct_Hybrid_TFSF_source_info
{
    char plane1_Ey[256], plane1_Ez[256], plane3_Ey[256], plane3_Ez[256];
    char plane2_Ex[256], plane2_Ez[256], plane4_Ex[256], plane4_Ez[256];
    char plane5_Ex[256], plane5_Ey[256], plane6_Ex[256], plane6_Ey[256];

    char plane1_Hy[256], plane1_Hz[256], plane3_Hy[256], plane3_Hz[256];
    char plane2_Hx[256], plane2_Hz[256], plane4_Hx[256], plane4_Hz[256];
    char plane5_Hx[256], plane5_Hy[256], plane6_Hx[256], plane6_Hy[256];

    double f0;

} struct_Hybrid_TFSF_source_info;

typedef struct struct_Hybrid_TFSF_File
{
    int Datalen;
    FILE* fp;
} struct_Hybrid_TFSF_File;

typedef struct struct_Hybrid_TFSF
{
    double Abs;
    double Phase;
    int ini1, ter1;
    int ini2, ter2;
    int planepos;
} struct_Hybrid_TFSF;

class TFDTD_Hybrid_TFSF : public TFDTD_Boundary_TFSF
{
public:
    __fastcall TFDTD_Hybrid_TFSF(struct_Field_carrier* pFc, coordinate coord);
    __fastcall ~TFDTD_Hybrid_TFSF();
public:
    void insert_parameters(struct_Hybrid_TFSF_source_info* pSource);
    void set_current_steps(int ss);
    void set_Hybrid_TFSF_box(struct_Hybrid_TFSF_source_info* pSource);
    //double get_plane_constant_x(FILE* fp);
protected:
    FILE* pf_debugHybrid;
private:
    double temp;
    double eta, beta, omega;
    double Distance, V_tilt, Time, Time0, Time1E, Time1H;
    int SPX, SPY, SPZ; // start points
    double x0, y0, z0; // 기준점 좌표

    int PXoffsetE, PYoffsetE, PZoffsetE;
    int PXoffsetE2, PYoffsetE2, PZoffsetE2;

    int PXoffsetH1, PYoffsetH1, PZoffsetH1;
	int PXoffsetH2, PYoffsetH2, PZoffsetH2;
    double PXminE, PYminE, PZminE;
    double PXminH, PYminH, PZminH;
    int nsteps;

    double k_hat_x;
    double k_hat_y;
    double k_hat_z;

	double* xcoordinate, * ycoordinate, * zcoordinate;
	double* xStaggered, * yStaggered, * zStaggered;

    Vec3 k_hat;     // 기본값: +x

    Vec3 r0;        // 기준점(phase/time reference)

    double alpha, T_0, n_T, Pulse;
    // m:minus, p: plus

	struct_Hybrid_TFSF_source_info* pSource;
    struct_Hybrid_TFSF_File pf_Ex_xmz, pf_Ex_xpz, pf_Ex_xym, pf_Ex_xyp;
    struct_Hybrid_TFSF_File pf_Ey_myz, pf_Ey_pyz, pf_Ey_xym, pf_Ey_xyp;
    struct_Hybrid_TFSF_File pf_Ez_myz, pf_Ez_pyz, pf_Ez_xmz, pf_Ez_xpz;

    struct_Hybrid_TFSF_File pf_Hx_xmz, pf_Hx_xpz, pf_Hx_xym, pf_Hx_xyp;
    struct_Hybrid_TFSF_File pf_Hy_myz, pf_Hy_pyz, pf_Hy_xym, pf_Hy_xyp;
    struct_Hybrid_TFSF_File pf_Hz_myz, pf_Hz_pyz, pf_Hz_xmz, pf_Hz_xpz;

    struct_Hybrid_TFSF** Ex_xmz, ** Ex_xpz, ** Ex_xym, ** Ex_xyp;
    struct_Hybrid_TFSF** Ey_myz, ** Ey_pyz, ** Ey_xym, ** Ey_xyp;
    struct_Hybrid_TFSF** Ez_myz, ** Ez_pyz, ** Ez_xmz, ** Ez_xpz;

    struct_Hybrid_TFSF** Hx_xmz, ** Hx_xpz, ** Hx_xym, ** Hx_xyp;
    struct_Hybrid_TFSF** Hy_myz, ** Hy_pyz, ** Hy_xym, ** Hy_xyp;
    struct_Hybrid_TFSF** Hz_myz, ** Hz_pyz, ** Hz_xmz, ** Hz_xpz;
    
    
    void Load_tangential_TFSF_fields();
    void Calculate_start_pointE();
    void Calculate_end_pointE();
    void Calculate_start_pointH();
    void Calculate_end_pointH();
    void Calculate_start_point(struct_Hybrid_TFSF_File psHTFSFFile, double* PXmin, double* PYmin, double* PZmin);
    void Calculate_end_point(struct_Hybrid_TFSF_File psHTFSFFile, double* PXmax, double* PYmax, double* PZmax);
    void Calculate_2D_size(struct_Hybrid_TFSF_File psHTFSFFile, struct_Hybrid_TFSF*** psHTFSF);
    void Calculate_2D_size(struct_Hybrid_TFSF_File psHTFSFFile, struct_Hybrid_TFSF*** psHTFSF, string HTFSFname);

    // calculation 2D size of EH fields on xy, yz, xz planes
    void Calculate_2D_size_EyYZ(struct_Hybrid_TFSF_File psHTFSFFile, struct_Hybrid_TFSF*** psHTFSF, string HTFSFname);
    void Calculate_2D_size_EzYZ(struct_Hybrid_TFSF_File psHTFSFFile, struct_Hybrid_TFSF*** psHTFSF, string HTFSFname);
    void Calculate_2D_size_HyYZ(struct_Hybrid_TFSF_File psHTFSFFile, struct_Hybrid_TFSF*** psHTFSF, string HTFSFname);
    void Calculate_2D_size_HzYZ(struct_Hybrid_TFSF_File psHTFSFFile, struct_Hybrid_TFSF*** psHTFSF, string HTFSFname);

	void Calculate_2D_size_ExXZ(struct_Hybrid_TFSF_File psHTFSFFile, struct_Hybrid_TFSF*** psHTFSF, string HTFSFname);
	void Calculate_2D_size_EzXZ(struct_Hybrid_TFSF_File psHTFSFFile, struct_Hybrid_TFSF*** psHTFSF, string HTFSFname);
	void Calculate_2D_size_HxXZ(struct_Hybrid_TFSF_File psHTFSFFile, struct_Hybrid_TFSF*** psHTFSF, string HTFSFname);
	void Calculate_2D_size_HzXZ(struct_Hybrid_TFSF_File psHTFSFFile, struct_Hybrid_TFSF*** psHTFSF, string HTFSFname);

	void Calculate_2D_size_ExXY(struct_Hybrid_TFSF_File psHTFSFFile, struct_Hybrid_TFSF*** psHTFSF, string HTFSFname);
	void Calculate_2D_size_EyXY(struct_Hybrid_TFSF_File psHTFSFFile, struct_Hybrid_TFSF*** psHTFSF, string HTFSFname);
	void Calculate_2D_size_HxXY(struct_Hybrid_TFSF_File psHTFSFFile, struct_Hybrid_TFSF*** psHTFSF, string HTFSFname);
	void Calculate_2D_size_HyXY(struct_Hybrid_TFSF_File psHTFSFFile, struct_Hybrid_TFSF*** psHTFSF, string HTFSFname);

    double circular_mean_phase(struct_Hybrid_TFSF** A, int Ny, int Nz);
    void set_plane_phase_constant(struct_Hybrid_TFSF** A, int Ny, int Nz, double phi);
    void add_plane_phase_offset(struct_Hybrid_TFSF** A, int Ny, int Nz, double dphi);
    double get_plane_constant_x(FILE* fp);
    double get_plane_constant_y(FILE* fp);
    double get_plane_constant_z(FILE* fp);

    void SetIncidentDirection_ThetaPhi(double theta_deg, double phi_deg);
    void SetBeta_FromFreq(double freq_hz);
    double proj_s(double x, double y, double z);
    double estimate_dphi_from_overlap(
        const std::unordered_map<Key3, cd, Key3Hash>& ref_vals,
        const std::unordered_map<Key3, cd, Key3Hash>& face_vals
    );
    double estimate_dphi_line(
        int N,                               // 교집합 선의 길이
        const std::complex<double>* cref,    // 기준 면의 phasor 배열 (선)
        const std::complex<double>* cface    // 맞출 면의 phasor 배열 (선)
    );
    bool Estimate_khat_from_faces(Vec3& out_khat);


    void append_samples_yz_Ez(
        std::vector<Sample3D>& out,
        struct_Hybrid_TFSF** Ez_yz,
        int Ny, int Nz,
        double x_const, double PYmin, double PZmin,
        double dy, double dz
    );
    void append_samples_xz_Ez(
        std::vector<Sample3D>& out,
        struct_Hybrid_TFSF** Ez_xz,   // [i][k]
        int Nx, int Nz,
        double PXmin, double PZmin,
        double y_const,
        double dx, double dz
    );
    void append_samples_xy_Ey(
        std::vector<Sample3D>& out,
        struct_Hybrid_TFSF** Ey_xy,   // [i][j]
        int Nx, int Ny,
        double PXmin, double PYmin,
        double z_const,
        double dx, double dy
    );

    bool Get_2D_dims_only(
        struct_Hybrid_TFSF_File f,
        int& outN1, int& outN2,
        int& outPlaneAxis // 0: yz(x const), 1: xz(y const), 2: xy(z const)
    );

    void process_samples(const std::vector<Sample3D>& data, NormalEq& ne, long long& total_eq);

    double Get_PulseE(double ss, struct_Hybrid_TFSF psHTFSF);
    double Get_PulseH(double ss, struct_Hybrid_TFSF psHTFSF);
    double Get_PulseE(double ss, struct_Hybrid_TFSF psHTFSF, double time);
    double Get_PulseH(double ss, struct_Hybrid_TFSF psHTFSF, double time);
    double Get_PulseE(int nsteps, const struct_Hybrid_TFSF& ps, double s, double t);
    double Get_PulseE(double ss);
    double Get_PulseH(double ss, double spatial);

    void Calcuate_inc_Field();
    void Calcuate_inc_Field_validate();
    int count_line(FILE* pfCount);

    double get_abs(double re, double im);
    double get_Phase(double y, double x);

    int round(double nominator, double denominator);

	int findIndex(double val, char coordinate);
	int findStaggeredIndex(double val, char coordinate);

};


