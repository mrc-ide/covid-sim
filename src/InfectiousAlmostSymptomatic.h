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

	void DoFalseCase(int ai, double t, unsigned short int ts, int tn);
	void DoDetectedCase(int ai, double t, unsigned short int ts, int tn);
	void DoCase(int ai, double t, unsigned short int ts, int tn);
	void DoPlaceClose(int i, int j, unsigned short int ts, int tn, int DoAnyway);
	void DoPlaceOpen(int i, int j, unsigned short int ts, int tn);
	void DoProphNoDelay(int ai, unsigned short int ts, int tn, int nc);
	int DoVacc(int ai, unsigned short int ts);
	void DoVaccNoDelay(int ai, unsigned short int ts);
	void DoTreatCase(int ai, unsigned short int ts, int tn);
	void DoProph(int ai, unsigned short int ts, int tn);

	void DoILI(int ai, int tn);
	void DoMild(int ai, int tn);
	void DoSARI(int ai, int tn);
	void DoCritical(int ai, int tn);
	void DoRecoveringFromCritical(int ai, int tn);
	void DoDeath_FromCriticalorSARIorILI(int ai, int tn);
	void DoRecover_FromSeverity(int ai, int tn);
	
	



};
