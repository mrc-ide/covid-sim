#pragma once
#include "IInfectionState.h"
#include "Param.h"
class Latent :
	public IInfectionState
{
	Param P;
public:
	Latent(Param p)
	{
		P = p;
	}
	virtual void GetsWorse(int ai, double t, int tn, int run);
	virtual void GetsBetter(int ai, double t, int tn, int run);

private:
	int DoIncub(int ai, unsigned short int ts, int tn, int run);
};

