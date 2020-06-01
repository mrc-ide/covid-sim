#include "Model.h"
#include "ModelMacros.h"
#include <math.h>
#include "Rand.h"
#include <stdio.h>
#include "Bitmap.h"
#include "InfectiousAlmostSymptomatic.h"
#include "Case.h"


void DoPlaceClose(int i, int j, unsigned short int ts, int tn, int DoAnyway);
void DoProph(int ai, unsigned short int ts, int tn);
void DoTreatCase(int ai, unsigned short int ts, int tn);
int DoVacc(int ai, unsigned short int ts);

void DoImmune(int ai)

{
	// This transfers a person straight from susceptible to immune. Used to start a run with a partially immune population.
	Person* a = Hosts + ai;

	a->infectionState->GetsBetter(ai, 0, 0, 0);
}

void DoInfect(int ai, double t, int tn, int run) // Change person from susceptible to latently infected.  added int as argument to DoInfect to record run number: ggilani - 15/10/14
{
	Person* a = Hosts + ai;

	a->infectionState->GetsWorse(ai, t, tn, run);

}

void DoIncub(int ai, unsigned short int ts, int tn, int run)
{
	Person* a = Hosts + ai;

	a->infectionState->GetsWorse(ai, ts, tn, run);
}

// stubs

void DoFalseCase(int ai, double t, unsigned short int ts, int tn)
{
	((InfectiousAlmostSymptomatic*)Hosts->stateHandlers[InfStatType_InfectiousAlmostSymptomatic])->DoFalseCase(ai, t, ts, tn);
}


void DoSARI(int ai, int tn)
{
	((InfectiousAlmostSymptomatic*)Hosts->stateHandlers[InfStatType_InfectiousAlmostSymptomatic])->DoSARI(ai, tn);
}



void DoCritical(int ai, int tn)
{
	((InfectiousAlmostSymptomatic*)Hosts->stateHandlers[InfStatType_InfectiousAlmostSymptomatic])->DoCritical(ai, tn);
}



void DoRecoveringFromCritical(int ai, int tn)
{
	((InfectiousAlmostSymptomatic*)Hosts->stateHandlers[InfStatType_InfectiousAlmostSymptomatic])->DoRecoveringFromCritical(ai, tn);
}

void DoDetectedCase(int ai, double t, unsigned short int ts, int tn)
{
	((InfectiousAlmostSymptomatic*)Hosts->stateHandlers[InfStatType_InfectiousAlmostSymptomatic])->DoDetectedCase(ai, t, ts, tn);
}

void DoProph(int ai, unsigned short int ts, int tn)
{
	((InfectiousAlmostSymptomatic*)Hosts->stateHandlers[InfStatType_InfectiousAlmostSymptomatic])->DoProph(ai, ts, tn);
}

void DoProphNoDelay(int ai, unsigned short int ts, int tn, int nc)
{
	((InfectiousAlmostSymptomatic*)Hosts->stateHandlers[InfStatType_InfectiousAlmostSymptomatic])->DoProphNoDelay(ai, ts, tn ,nc);
}

int DoVacc(int ai, unsigned short int ts)
{
	return ((InfectiousAlmostSymptomatic*)Hosts->stateHandlers[InfStatType_InfectiousAlmostSymptomatic])->DoVacc(State.mvacc_queue[ai], ts);
}


void DoVaccNoDelay(int ai, unsigned short int ts)
{
	((InfectiousAlmostSymptomatic*)Hosts->stateHandlers[InfStatType_InfectiousAlmostSymptomatic])->DoVaccNoDelay(State.mvacc_queue[ai], ts);
}

void DoPlaceOpen(int i, int j, unsigned short int ts, int tn)
{
	((InfectiousAlmostSymptomatic*)Hosts->stateHandlers[InfStatType_InfectiousAlmostSymptomatic])->DoPlaceOpen(i, j, ts, tn);
}

void DoPlaceClose(int i, int j, unsigned short int ts, int tn, int doAnyway)
{
	((InfectiousAlmostSymptomatic*)Hosts->stateHandlers[InfStatType_InfectiousAlmostSymptomatic])->DoPlaceClose(i, j, ts, tn, doAnyway);
}


void DoRecover(int ai, int tn, int run)
{
	((Case*)Hosts->stateHandlers[InfStatType_Case])->DoRecover(ai, tn, run);
}

void DoRecover_FromSeverity(int ai, int tn)
{
	((InfectiousAlmostSymptomatic*)Hosts->stateHandlers[InfStatType_InfectiousAlmostSymptomatic])->DoRecover_FromSeverity(ai, tn);
}


void DoDeath_FromCriticalorSARIorILI(int ai, int tn)
{
	((InfectiousAlmostSymptomatic*)Hosts->stateHandlers[InfStatType_InfectiousAlmostSymptomatic])->DoDeath_FromCriticalorSARIorILI(ai, tn);
}


