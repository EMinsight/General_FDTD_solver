#pragma once
#include "TFunctions.h"

class TEnvelop_Extract
{
public:
	TEnvelop_Extract();
	~TEnvelop_Extract();
	static bool Chiral_Envelop_Extract(struct_Field_carrier* pFc, Material_properties* pMp, int*** pSurface, string mName);
	// Global (for current)

};

