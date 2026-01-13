#include "TPostprocessing.h"

TPostprocessing::TPostprocessing(struct_Field_carrier* pFc, Material_properties* pMp)
{
	m_pi = acos(-1.0);
	eps_0 = 1e-9 / 36 / m_pi;
	mu_0 = 4 * m_pi * 1e-7;
	c_0 = 1 / sqrt(eps_0 * mu_0);
	strcpy(JsonFileName, "./IEM_solver_input_format.json");
	this->pFc = pFc;
	this->pMp = pMp;
}

void TPostprocessing::Prepare2save()
{
	Read_json();
	Assign_output();
}

void TPostprocessing::Check_infomations_null()
{

}

void TPostprocessing::Read_json()
{
	// Document Jdocument;
	ifstream ifs(JsonFileName);
	IStreamWrapper isw(ifs);
	Jdocument.ParseStream(isw);

	if (Jdocument["SIMULATION_ENVIORONMENT"].HasMember("Observation") == true)
		no_points = Jdocument["SIMULATION_ENVIORONMENT"]["Observation"].Size();

	if (Jdocument["SIMULATION_ENVIORONMENT"].HasMember("Distribution") == true)
		no_planes = Jdocument["SIMULATION_ENVIORONMENT"]["Distribution"].Size();

	if (Jdocument["SIMULATION_ENVIORONMENT"].HasMember("Composition") == true)
		// 추후 Solver 에서 Location 삭제
		// Json에서 Position추가.
		no_comp = Jdocument["SIMULATION_ENVIORONMENT"]["Composition"].Size();

	if (Jdocument["SIMULATION_ENVIORONMENT"].HasMember("Current") == true)
		no_current = Jdocument["SIMULATION_ENVIORONMENT"]["Current"].Size();

	if (Jdocument["SIMULATION_ENVIORONMENT"].HasMember("OpticalC") == true)
		no_opticalC = Jdocument["SIMULATION_ENVIORONMENT"]["OpticalC"].Size();

	if (Jdocument["SIMULATION_ENVIORONMENT"].HasMember("Chargedensity") == true)
		no_chargeD = Jdocument["SIMULATION_ENVIORONMENT"]["Chargedensity"].Size(); ;

	// cout << " no_points :" << no_points << " no_planes :" << no_planes << " no_comp :" << no_comp << " no_current :" << no_current << " no_opticalC :" << no_opticalC << " no_chargeD :" << no_chargeD;
	ifs.close();
}

void TPostprocessing::Assign_output()
{
	if (no_points != 0)
		setup_saveObservation();
	if (no_comp != 0)
		setup_saveComposit();
	if (no_planes != 0)
		setup_saveDistribution();
	if (no_opticalC != 0)
		setup_optical_chirality(); 
	if (no_current != 0)
		setup_SaveCurrentDistribution(); 
	if (no_chargeD != 0)
		setup_Charge_Density();

}

void TPostprocessing::Save(int ss)
{
	if (no_points != 0)
		pSaveObs->Save(ss);
	if (no_comp != 0)
		pSaveComp->Save(ss);
	// 전체 저장 (pulse) or 특정영역 저장 (harmonic)
	if (no_planes != 0) 
		pSaveDist->Save(ss);
	// harmonic
	if (no_opticalC != 0)
		pOpticalC->Save(ss);
	// 전체 저장 (pulse) or 특정영역 저장 (harmonic)
	if (no_current != 0)
		pSaveCur_tDist->Save(ss);
	if (no_chargeD != 0)
		pChargeD->Save(ss);
}

void TPostprocessing::setup_saveObservation() 
{
	pSaveObs = new TSaveObservation(); 
	pSaveObs->set_postprocessing_info(pFc, pMp);
	pSaveObs->Prepare2save();
}

void TPostprocessing::setup_saveComposit() 
{
	pSaveComp = new TSaveComposit();  
	pSaveComp->set_postprocessing_info(pFc, pMp);
	pSaveComp->Prepare2save();
}

void TPostprocessing::setup_saveDistribution() 
{
	// 전체 저장 (pulse) or 특정영역 저장 (harmonic)
	pSaveDist = new TSaveDistribution(); 
	pSaveDist->set_postprocessing_info(pFc, pMp);
	pSaveDist->Prepare2save();// 여기서 처리.
}

void TPostprocessing::setup_optical_chirality() 
{
	pOpticalC = new Toptical_chirality(); 
	pOpticalC->set_postprocessing_info(pFc, pMp);
	pOpticalC->Prepare2save();// 여기서 처리.
}

void TPostprocessing::setup_SaveCurrentDistribution() 
{
	// 전체 저장 (pulse) or 특정영역 저장 (harmonic)
	pSaveCur_tDist = new TSaveCurrentDistribution(); 
	pSaveCur_tDist->set_postprocessing_info(pFc, pMp);
}

void TPostprocessing::setup_Charge_Density()
{
	pChargeD = new TSaveChargeDensity();
	pChargeD->set_postprocessing_info(pFc, pMp);
	pChargeD->Prepare2save();
}