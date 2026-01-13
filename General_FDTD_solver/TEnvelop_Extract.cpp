#include "TEnvelop_Extract.h"

bool TEnvelop_Extract::Chiral_Envelop_Extract(struct_Field_carrier* pFc, Material_properties* pMp, int*** pSurface, string Mname)
{
	/*
	* Chiral_Envelop_Extract(struct_Field_carrier *pFc, Material_properties *pMp, double ***pSurface, char Mname[16]);
	* pFc : Field carrier
	* pMp : Material_properties
	* pSurface : envelop of the material
	* Mname : full name of mateiral
	*/
	int txsize = pFc->txsize; 
	int tysize = pFc->tysize;
	int tzsize = pFc->tzsize;
	int*** Media = pFc->media;

	//Get_Pulse_infomation();
	//Layer = no_composit - 2 + 3;
	// Layer = 7;
	// Num. Layer check;
	int ** pExtract_Optical_Chirality_plane = Tptr2<int>(txsize, tysize);
	int ** pExtract_chiral_envelope = Tptr2<int>(txsize, tysize);
	for (int nn = 0; nn < 2; nn++)
		for (int ii = 0; ii < txsize; ii++)
			memset(pSurface[nn][ii], 0, sizeof(int) * tysize);
	// ±¼°î °æ°è ÃßÃâ ¾Ë°í¸®Áò.
	bool flag = false;
	for (int ii = 0; ii < txsize; ii++)
	for (int jj = 0; jj < tysize; jj++)
	for (int kk = 1; kk < tzsize - 1; kk++)
	{
		int mat_ = Media[ii][jj][kk];
		int mat_1p = Media[ii][jj][kk+1];
		int mat_1n = Media[ii][jj][kk-1];
		if ((mat_ != mat_1p)&&
			strcmp(pMp[mat_].Name, Mname.c_str()) == 0)  // Gold_up
		{
			pSurface[0][ii][jj] = kk;
			flag = true;
		}
		else if ((mat_ != mat_1n) &&
			strcmp(pMp[mat_].Name, Mname.c_str()) == 0)  // Gold_down
		{
			pSurface[1][ii][jj] = kk;
			flag = true;
		}
	}

	return flag;
}
