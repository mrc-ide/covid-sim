#include "InfectiousAsymptomaticNotCase.h"
#include "Model.h"
#include <math.h>
#include "ModelMacros.h"
#include "Bitmap.h"
#include "Dead_WasAsymp.h"
// get worse
/*
Dead
Dead_WasAsymp
*/

// get better
/*
Recovered
RecoveredFromAsymp
*/

void InfectiousAsymptomaticNotCase::GetsWorse(int ai, double t, int tn, int run)
{
	DoDeath(ai, tn, run);

	Person* a = Hosts + ai;
	a->infectionState = new Dead_WasAsymp();
}

void InfectiousAsymptomaticNotCase::GetsBetter(int ai, double t, int tn, int run)
{
}

void InfectiousAsymptomaticNotCase::DoDeath(int ai, int tn, int run)
{
	int i, x, y;
	Person* a = Hosts + ai;

	if ((a->inf == InfStat_InfectiousAsymptomaticNotCase || a->inf == InfStat_Case))
	{
		a->inf = (InfStat)(InfStat_Dead * a->inf / abs(a->inf));
		Cells[a->pcell].D++;
		Cells[a->pcell].I--;
		i = a->listpos;
		if (i < Cells[a->pcell].S + Cells[a->pcell].L + Cells[a->pcell].I)
		{
			Cells[a->pcell].susceptible[a->listpos] = Cells[a->pcell].infected[Cells[a->pcell].I];
			Hosts[Cells[a->pcell].susceptible[a->listpos]].listpos = i;
			a->listpos = Cells[a->pcell].S + Cells[a->pcell].L + Cells[a->pcell].I;
			Cells[a->pcell].susceptible[a->listpos] = ai;
		}

		/*		a->listpos=-1; */
		StateT[tn].cumDa[HOST_AGE_GROUP(ai)]++;
		if (P->DoAdUnits) StateT[tn].cumD_adunit[Mcells[a->mcell].adunit]++;
		if (P->OutputBitmap)
		{
			if ((P->OutputBitmapDetected == 0) || ((P->OutputBitmapDetected == 1) && (Hosts[ai].detected == 1)))
			{
				x = ((int)(Households[a->hh].loc_x * P->scalex)) - P->bminx;
				y = ((int)(Households[a->hh].loc_y * P->scaley)) - P->bminy;
				if ((x >= 0) && (x < P->bwidth) && (y >= 0) && (y < P->bheight))
				{
					unsigned j = y * bmh->width + x;
					if (j < bmh->imagesize)
					{
#pragma omp atomic
						bmRecovered[j]++;
#pragma omp atomic
						bmInfected[j]--;
					}
				}
			}
		}
	}
}
