#include "InfectiousAsymptomaticNotCase.h"
#include "Model.h"

void InfectiousAsymptomaticNotCase::GetsWorse(int ai, double t, int tn, int run)
{
	DoDeath(ai, tn, run);

	Person* a = Hosts + ai;
	a->infectionState = Hosts->stateHandlers[InfStatType_Dead_WasAsymp];
}

void InfectiousAsymptomaticNotCase::GetsBetter(int ai, double t, int tn, int run)
{
	DoRecover(ai, tn, run);

	Person* a = Hosts + ai;
	a->infectionState = Hosts->stateHandlers[InfStatType_RecoveredFromAsymp];
}
