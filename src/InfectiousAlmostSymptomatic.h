#pragma once
#include "IInfectionState.h"
#include "Param.h"

class InfectiousAlmostSymptomatic : public IInfectionState
{
public:
	InfectiousAlmostSymptomatic(Param* p)
		: IInfectionState(p)
	{
		
	}
public:
	virtual void GetsWorse(int ai, double t, int tn, int run);
	virtual void GetsBetter(int ai, double t, int tn, int run);

private:
	void DoCase(int ai, double t, unsigned short int ts, int tn);
	void DoILI(int ai, int tn);

};
