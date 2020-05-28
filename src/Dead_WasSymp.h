#include "IInfectionState.h"
#include "Param.h"
class Dead_WasSymp :
	public IInfectionState
{
public:
	Dead_WasSymp()
	{
	}
	virtual void GetsWorse(int ai, double t, int tn, int run);
	virtual void GetsBetter(int ai, double t, int tn, int run);
};

