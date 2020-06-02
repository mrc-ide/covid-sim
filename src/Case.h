#include "IInfectionState.h"
#include "Param.h"
#include "InfectiousBase.h"

class Case :
	public IInfectionState,
	protected InfectiousBase
{
public:

	Case(Param* p)
		: IInfectionState(p), InfectiousBase(p)
	{
	}
	virtual void GetsWorse(int ai, double t, int tn, int run);
	virtual void GetsBetter(int ai, double t, int tn, int run);

};
