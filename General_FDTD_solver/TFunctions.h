#ifndef TFUNCTIONSH
#define TFUNCTIONSH
#include <memory>
#include <complex>
#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <direct.h>
#include <omp.h>
#include <vector>

#include "../rapidjsonLibrary/document.h"
#include "../rapidjsonLibrary/writer.h"
#include "../rapidjsonLibrary/stringbuffer.h"
#include "../rapidjsonLibrary/filereadstream.h"
#include "../rapidjsonLibrary/filewritestream.h"
#include "../rapidjsonLibrary/istreamwrapper.h"
#include "../rapidjsonLibrary/ostreamwrapper.h"

using namespace rapidjson;
using namespace std;

// 모든 클래스가 참조하는 공용 구조체 / 함수. Mesh/Material/input/output 구조체도 여기서 만든다.
typedef struct Material_properties
{
	char TypeOfMaterials;
	// dispersion materials
	double* Ca, * Cb, * Cc, * Cd, * Ce;
	// Non-dispersion materials
	double Epsr, Sigma, loss_tan;
	char Name[512];
	//int Material_index, Material_type;
} Material_properties;

typedef struct struct_Field_carrier
{
	double*** Ex_address, *** Ey_address, *** Ez_address;
	double*** Ex_p_address, *** Ey_p_address, *** Ez_p_address;
	double*** Ex_pp_address, *** Ey_pp_address, *** Ez_pp_address;
	double*** Hx_address, *** Hy_address, *** Hz_address;
	double**** Jx_address, **** Jy_address, **** Jz_address;
	double**** Jx_p_address, **** Jy_p_address, **** Jz_p_address;
	int*** media, *** Media_Ex, *** Media_Ey, *** Media_Ez;

	int txsize, tysize, tzsize;
	double* dx, * dy, * dz;
	double* E_dx, * E_dy, * E_dz;
	double* H_dx, * H_dy, * H_dz;
	int PML[3];
	double dt;
	double NSTEPS;
	double f0, BW;
	int total_pol;
	int aa1m1, aa2m2; // BC (PBC) option

	// thread
	int num_thread;
} struct_Field_carrier;

typedef struct coordinate
{
	double* xcoordinate, * ycoordinate, * zcoordinate;
	double* xStaggered, * yStaggered, * zStaggered;
} coordinate;

// 여기는 API가 출력하는 구조체들. 이거를 맞춰야 한다. 
struct MeshInput {
	int tx, ty, tz;
	std::vector<double> x, y, z;     // boundary coords: len=tx+1 ...
	std::vector<int> media;          // len=tx*ty*tz (material key)
	int PML[3];
	double dt;
	int nsteps;
	double unit_scale;               // UNIT 변환용
};

struct MaterialEntry {
	int key;
	int type;                        // nondispersive/dispersion...
	std::string name;
	double epsr, sigma, loss_tan;
	std::vector<double> Ca, Cb, Cc, Cd, Ce; // dispersion일 때만
};

struct MaterialTableInput {
	std::vector<MaterialEntry> entries;
};

struct Waveform {
	std::string name;
	int type;
	double f0, bw, phase;
};

struct Tx {
	std::string waveform_name;
	int type; // point/plane/...
	double x, y, z;  // 또는 box
	double theta, phi, polar;
	double amplitude;
};

struct Rx {
	int type; // point/plane
	// 위치 + 저장 요청(Ey, Ez 등)
};

struct TxRxInput {
	std::vector<Waveform> waveforms;
	std::vector<Tx> txs;
	std::vector<Rx> rxs;
};

// 이런 데이터를 출력한다. solver의 입력부를 여기에 맞춰야 한다. 

template <class T> T conv(T *Aresp, T *BFunc, int tt)
{
    T sumConvG = 0.0;
    for(int tau = 0 ; tau < tt ; tau++)
    {
        sumConvG += Aresp[tau]*BFunc[tt-tau];
    }
    return sumConvG;
}

template <class T> T min1D(T* Array)
{
	int ArraySize = 0;
	while (1)
	{
		if (Array[ArraySize] == -1)
			break;
		ArraySize++;
	}

	double min = 1e30;
	for (int ii = 0; ii < ArraySize; ii++)
	{
		if (min > Array[ii])
			min = Array[ii];
	}

	return min;
}

// Complex form used only to update J fields 
template <class T> T**** Tptr4(int x, int y, int z, int pol)
{
	T**** temp;
	temp = new T ***[x];
	for (int i = 0; i < x; i++)
	{
		temp[i] = new T **[y];
		for (int j = 0; j < y; j++)
		{
			temp[i][j] = new T *[z];
			for (int k = 0; k < z; k++)
			{
				temp[i][j][k] = new T[pol];
				memset(temp[i][j][k], 0, sizeof(T) * pol);
			}
		}
	}
	return temp;
}


template <class T> T ***Tptr3(int x, int y, int z)
{
	T ***temp;
	temp = new T **[x];
	for (int i = 0; i < x; i++)
	{
		temp[i] = new T *[y];
		for (int j = 0; j < y; j++)
		{
			temp[i][j] = new T[z];
			memset(temp[i][j], 0, sizeof(T)*z);
		}
	}
	return temp;
}

template <class T> T **Tptr2(int x, int y)
{
	T **temp;
	temp = new T*[x];
	for (int i = 0; i < x; i++)
	{
		temp[i] = new T[y];
		memset(temp[i], 0, sizeof(T)*y);
	}
	return temp;
}

template <class T> T** Identity_matrix_2D(int x, int y)
{
	T** temp;
	temp = new T * [x];
	for (int i = 0; i < x; i++)
	{
		temp[i] = new T[y];
		memset(temp[i], 0, sizeof(T) * y);
	}
	temp[0][0] = 1;
	temp[1][1] = 1;
	temp[2][2] = 1;
	return temp;
}

template <class T> T *Tptr11(int x)
{
	T *temp;
	temp = new T[x];
	memset(temp, 1, sizeof(T)*x);

	return temp;
}

template <class T> T *Tptr1(int x)
{
	T *temp;
	temp = new T[x];
	memset(temp, 0, sizeof(T)*x);

	return temp;
}

template <typename T>
void memfree(T **corr_field)
{
	T **Free_field;
	Free_field = corr_field;
	delete [] *Free_field;
}

template <typename T>
void memfree(T ***corr_field, int x)
{
	T ***Free_field;
	Free_field = corr_field;
	for (int i = 0 ; i < x ; i++)
		delete [] *(*Free_field+i);
	delete [] *Free_field;
}

template <typename T>
void memfree(T ****corr_field, int x, int y)
{
	T ****Free_field;
	Free_field = corr_field;

	for (int i = 0 ; i < x ; i++)
	{
		for (int j = 0 ; j < y ; j++)
			delete [] *(*(*Free_field+i)+j);
		delete [] *(*Free_field+i);
	}
	delete [] *Free_field;
}

template <typename T>
void memfree(T***** corr_field, int x, int y, int z)
{
	T***** Free_field;
	Free_field = corr_field;
	for (int i = 0; i < x; i++)
	{
		for (int j = 0; j < y; j++)
		{
			for (int k = 0; k < z; k++)
				delete[] *(*(*(*Free_field +i)+ j) + k);
			delete[] *(* (*Free_field + i) + j);
		}
		delete[] * (*Free_field + i);
	
	}
	delete[] * Free_field;
}

#endif