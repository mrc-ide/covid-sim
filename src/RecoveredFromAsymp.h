#include "IInfectionState.h"
#include "Param.h"
class RecoveredFromAsymp :
	public IInfectionState
{
public:
	RecoveredFromAsymp(Param* p)
		: IInfectionState(p)
	{
	}
	virtual void GetsWorse(int ai, double t, int tn, int run);
	virtual void GetsBetter(int ai, double t, int tn, int run);
};
