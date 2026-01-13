#pragma once

// main���� ȣ���ϴ� Ŭ������, postprocessing�� �����͸� �����ϰ� 
// �����Ϸ��� ������ json���κ��� �о�� �� 
// postprocessing�� �ʿ��� Ŭ�������� ȣ��


#include "TFunctions.h"
#include "TSaveObservation.h"
#include "TSaveDistribution.h"
#include "TSaveComposit.h"
#include "Toptical_chirality.h"
#include "TSaveCurrentDistribution.h"
#include "TChargeDensity.h"

//using namespace std;

class TPostprocessing 
{
public:
	TPostprocessing() {};
	TPostprocessing(struct_Field_carrier* pFc, Material_properties* pMp);
	~TPostprocessing() {};
	
	void Save(int ss);
	void Prepare2save();

private:
	Document Jdocument;
	char JsonFileName[1024];
	double m_pi, eps_0, mu_0, c_0;
	struct_Field_carrier* pFc;
	Material_properties* pMp;
	int no_points = 0, no_planes = 0, no_comp = 0, no_current = 0, no_opticalC = 0, no_chargeD = 0;
	
	TSaveObservation* pSaveObs;
	TSaveComposit* pSaveComp;
	TSaveDistribution* pSaveDist;
	TSaveCurrentDistribution* pSaveCur_tDist;
	Toptical_chirality* pOpticalC;
	TSaveChargeDensity* pChargeD;
	
	void Read_json();
	void Assign_output();

	void setup_saveObservation();
	void setup_saveComposit();
	void setup_saveDistribution();
	void setup_optical_chirality();
	void setup_SaveCurrentDistribution();
	void setup_Charge_Density();

	void Check_infomations_null();
};