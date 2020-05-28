#include "IInfectionState.h"
#include "Param.h"
class ImmuneAtStart :
	public IInfectionState
{
public:
	ImmuneAtStart()
	{
	}
	virtual void GetsWorse(int ai, double t, int tn, int run);
	virtual void GetsBetter(int ai, double t, int tn, int run);
};
