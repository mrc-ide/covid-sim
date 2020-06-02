#pragma once
#include "IInfectionState.h"
#include "InfectiousBase.h"
#include "Param.h"

class InfectiousAsymptomaticNotCase : public IInfectionState,
	protected InfectiousBase
{

public:
	InfectiousAsymptomaticNotCase(Param* p)
		: IInfectionState(p), InfectiousBase(p)
	{
	}
public:
	virtual void GetsWorse(int ai, double t, int tn, int run);
	virtual void GetsBetter(int ai, double t, int tn, int run);
};

