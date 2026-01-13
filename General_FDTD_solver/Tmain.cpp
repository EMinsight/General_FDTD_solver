#include "FDTD_Terminal.h"

int main()
{
	// watermark 추가
	/*
	here.
	*/
	//
	// git 주석 추가.

	Tmain *pTmain;

	pTmain = new Tmain();
	pTmain->TFDTD_main_routine();

	

	delete pTmain;
	return 0;
}