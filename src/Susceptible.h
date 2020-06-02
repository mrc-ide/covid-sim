#pragma once
#include "Param.h"
#include "IInfectionState.h"
#include "Bitmap.h"
class Susceptible : public IInfectionState
{

public:

	Susceptible(Param* p)
		: IInfectionState(p)
	{
	}
	virtual void GetsWorse(int ai, double t, int tn, int run);
	virtual void GetsBetter(int ai, double t, int tn, int run);

private:

	void BecomesImmune(int ai);
	void BecomesInfected(int ai, double t, int tn, int run);
	void RecordEvent(double t, int ai, int run, int type, int tn);
	

};

