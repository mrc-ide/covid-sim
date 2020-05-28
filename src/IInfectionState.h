#pragma once

class IInfectionState
{
public:
	virtual void GetsWorse(int ai, double t, int tn, int run) = 0;
	virtual void GetsBetter(int ai, double t, int tn, int run) = 0;
};
