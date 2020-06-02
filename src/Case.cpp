#include "Case.h"
#include "Model.h"
#include "Dead_WasSymp.h"
#include "RecoveredFromSymp.h"
#include <math.h>
#include "Bitmap.h"
#include "ModelMacros.h"

void Case::GetsWorse(int ai, double t, int tn, int run)
{
	DoDeath(ai, tn, run);

	// TODO - consolidate
	Person* a = Hosts + ai;

	a->infectionState = Hosts->stateHandlers[InfStatType_Dead_WasSymp];
}

void Case::GetsBetter(int ai, double t, int tn, int run)
{
	DoRecover(ai, tn, run);

	Person* a = Hosts + ai;

	a->infectionState = Hosts->stateHandlers[InfStatType_RecoveredFromSymp];
}
