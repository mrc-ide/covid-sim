#include "IInfectionState.h"
#include "Param.h"
class RecoveredFromSymp :
	public IInfectionState
{
public:
	RecoveredFromSymp(Param* p)
		: IInfectionState(p)
	{
	}
	virtual void GetsWorse(int ai, double t, int tn, int run);
	virtual void GetsBetter(int ai, double t, int tn, int run);
};
