#include "Case.h"
#include "Model.h"

void Case::GetsWorse(int ai, double t, int tn, int run)
{
	DoDeath(ai, tn, run);

	Person* a = Hosts + ai;

	a->infectionState = Hosts->stateHandlers[InfStatType_Dead_WasSymp];
}

void Case::GetsBetter(int ai, double t, int tn, int run)
{
	DoRecover(ai, tn, run);

	Person* a = Hosts + ai;

	a->infectionState = Hosts->stateHandlers[InfStatType_RecoveredFromSymp];
}
