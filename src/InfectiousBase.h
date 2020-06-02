#pragma once
#include "IInfectionState.h"
#include "Param.h"

class InfectiousBase 
{

protected:

	void DoDeath(int ai, int tn, int run);
	void DoRecover(int ai, int tn, int run);


	Param* P;

	InfectiousBase(Param* p)
	{
		P = p;
	}
	InfectiousBase()
	{

	}
};

