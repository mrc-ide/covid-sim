#pragma once
#include "Param.h"
class IInfectionState
{
public:
	virtual void GetsWorse(int ai, double t, int tn, int run) {}
	virtual void GetsWorse(int ai, unsigned short ts, int tn, int run) {}

	virtual void GetsBetter(int ai, double t, int tn, int run) {}

	Param* P;
public:
	IInfectionState(Param* p)
	{
		P = p;
	}
	IInfectionState()
	{
			
	}

};
