#include "IInfectionState.h"
#include "Param.h"
class Dead_WasAsymp :
	public IInfectionState
{
public:
	Dead_WasAsymp()
	{
	}
	virtual void GetsWorse(int ai, double t, int tn, int run);
	virtual void GetsBetter(int ai, double t, int tn, int run);
};
