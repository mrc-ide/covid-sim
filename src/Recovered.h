#include "IInfectionState.h"
#include "Param.h"
class Recovered :
	public IInfectionState
{
public:
	Recovered(Param* p)
		: IInfectionState(p)
	{
	}
	virtual void GetsWorse(int ai, double t, int tn, int run);
	virtual void GetsBetter(int ai, double t, int tn, int run);
};
