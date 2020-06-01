#pragma once
#include "IInfectionState.h"
#include "Param.h"

class InfectiousAsymptomaticNotCase : public IInfectionState
{

public:
	InfectiousAsymptomaticNotCase(Param* p)
		: IInfectionState(p)
	{
	}
public:
	virtual void GetsWorse(int ai, double t, int tn, int run);
	virtual void GetsBetter(int ai, double t, int tn, int run);

private:

	void DoDeath(int ai, int tn, int run);
	void DoRecover(int ai, int tn, int run);

};

