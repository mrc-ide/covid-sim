#include "IInfectionState.h"
#include "Param.h"
class Case :
	public IInfectionState
{
public:
	Case()
	{
	}
	virtual void GetsWorse(int ai, double t, int tn, int run);
	virtual void GetsBetter(int ai, double t, int tn, int run);

private:
	void DoRecover(int ai, int tn, int run);
	void DoDeath(int ai, int tn, int run);
};
