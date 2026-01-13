//---------------------------------------------------------------------------
#ifndef TSource_ExcitationH
#define TSource_ExcitationH
//---------------------------------------------------------------------------
#include <memory.h>
#include <math.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "TFunctions.h"
#include "TFDTD_Boundary_TFSF.h"
#include "TFDTDparameters.h"
#include "TFDTD_Hybrid_mod_with_MOM.h"
#include "TFDTD_Hybrid_TFSF.h"

typedef struct struct_Source_type
{
	int Sx1, Sy1, Sz1, Sx2, Sy2, Sz2;
	double Rsx1, Rsy1, Rsz1, Rsx2, Rsy2, Rsz2;
	int point_dir;
	double impedance;

}  struct_Source_type;


int FindIndex(double value, char Coordinate);

class TSource_Excitation
{
public:
	//__fastcall TSource_Excitation(double*** Ex, double*** Ey, double*** Ez, double*** Hx, double*** Hy, double*** Hz);
	TSource_Excitation(struct_Field_carrier* pFc);
	TSource_Excitation();
	~TSource_Excitation();
protected:
	double eps0;
	double pi;
private:
	double* dx, * dy, * dz;
	int PMLX, PMLY, PMLZ;
	int txsize, tysize, tzsize;
	double Period;
	int NSTEPS;

	double ***Ex, ***Ey, ***Ez;
	double ***Hx, ***Hy, ***Hz;
	// No. Source
	int no_source;
	char JsonFileName[256];
	double* Phase;
	double *PW_Theta;
	double *PW_Fhi  ;
	double *PW_Polar;
	int PbufferX, PbufferY, PbufferZ;
	TFDTD_Boundary_TFSF** planewave;
	// Hybrid
	TFDTD_Hybrid_mod_with_MOM* pHybrid;
	TFDTD_Hybrid_TFSF** pHybrid_TFSF;
	struct_Field_carrier* pFc;
public:
	// double Get_Pulse(int ss);  // Define Pulse for planewave
	double Get_Pulse(int ss, double phi);
	
	double Get_sourceType() { return SourceType; };
	double get_f0(){return f0;};
	double get_BW(){return bw;};
	//int Set_Pulse_infomation();		// Read dat file for Pulse
	void Set_Pulse_parameters(); // Basic parameters for Pulse
	//void Get_Pulse_subplementary_parameters();	// additional parameters for Pulse, (not available)
	
	void Set_dsdt(double *dx, double *dy, double *dz, double dt){
		this->dx = dx;
		this->dy = dy;
		this->dz = dz;
		this->dt = dt;};
	int getNSTEPS() { return NSTEPS; };
	
//	double Define_Source(int ss); // Defile value of Pulse

// �������� ���
	int Get_Number_of_Source(){return no_source;};
    struct_Source_type *Get_source_data(){return sst;};
	struct_Hybrid_source_info* Get_hybrid_data() {return pHsi;};

	// int Get_Sx_point() { return sst[0].Sx1; };
	// int Get_Sy_point() { return sst[0].Sy1; };
	// int Get_Sz_point() { return sst[0].Sz1; };
	void Inject_SourceE(int ss);
	void Inject_SourceH(int ss);
	void Inject_Source_add(int n);
protected:
	void Get_User_define_source();
	double Set_User_define_source(int ss);
	double *pPulse;
	double saved_Pulse_length;
	void Set_TFSF();
	void Set_Hybrid_19();
	void Set_Hybrid_TFSF();
	
	void Inject_Edge_Source();    // Apply Point souorce
	//void Inject_TFSF_Source();
	void Inject_Pointsource();
	//void Inject_Hybrid();
protected:
	int SourceType; //1:Voltage, 2: Current, 3: Point, 4: Plane
	struct_Source_type* sst;
	struct_Hybrid_source_info* pHsi;
	struct_Hybrid_TFSF_source_info* pHybTFSF;
	coordinate coordinate;
	int Opt_Pulse;

	double* xcoordinate, * ycoordinate, * zcoordinate;
	double* xStaggered, * yStaggered, * zStaggered;

	double E0 ;
	double f0 ;
	double bw ;
	double T  ;
	double n_T;
	double T_0;
	double alpha;
	double Pulse;
	double x0, y0, z0;
private:
	double dt;
private:
};

#endif
