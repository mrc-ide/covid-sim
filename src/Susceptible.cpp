#include "Susceptible.h"
#include <math.h>
#include "Model.h"
#include "Latent.h"
#include "InfectiousAlmostSymptomatic.h"
#include "ModelMacros.h"
#include "Rand.h"
#include "Case.h"


extern void RecordEvent(double t, int ai, int run, int type, int tn);

void Susceptible::GetsWorse(int ai, double t, int tn, int run)
{
	BecomesInfected(ai, t, tn, run);
	Person* a = Hosts + ai;
	a->infectionState = Hosts->stateHandlers[InfStatType_Latent];
}

void Susceptible::GetsBetter(int ai, double t, int tn, int run)
{
	BecomesImmune(tn);
	Person* a = Hosts + ai;
	a->infectionState = Hosts->stateHandlers[InfStatType_ImmuneAtStart];
}

void Susceptible::BecomesImmune(int ai)
{
	// This transfers a person straight from susceptible to immune. Used to start a run with a partially immune population.
	Person* a;
	int c;
	int x, y;

	a = Hosts + ai;

	c = a->pcell;

	a->inf = InfStat_ImmuneAtStart;


	Cells[c].S--;
	if (a->listpos < Cells[c].S)
	{
		Cells[c].susceptible[a->listpos] = Cells[c].susceptible[Cells[c].S];
		Hosts[Cells[c].susceptible[a->listpos]].listpos = a->listpos;
	}
	if (Cells[c].L > 0)
	{
		Cells[c].susceptible[Cells[c].S] = Cells[c].susceptible[Cells[c].S + Cells[c].L];
		Hosts[Cells[c].susceptible[Cells[c].S]].listpos = Cells[c].S;
	}
	if (Cells[c].I > 0)
	{
		Cells[c].susceptible[Cells[c].S + Cells[c].L] = Cells[c].susceptible[Cells[c].S + Cells[c].L + Cells[c].I];
		Hosts[Cells[c].susceptible[Cells[c].S + Cells[c].L]].listpos = Cells[c].S + Cells[c].L;
	}
	if (a->listpos < Cells[c].S + Cells[c].L + Cells[c].I)
	{
		Cells[c].susceptible[Cells[c].S + Cells[c].L + Cells[c].I] = ai;
		a->listpos = Cells[c].S + Cells[c].L + Cells[c].I;
	}
	Cells[c].latent--;
	Cells[c].infected--;
	Cells[c].R++;
	if (P->OutputBitmap)
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
			}
		}
	}
}

void Susceptible::BecomesInfected(int ai, double t, int tn, int run) // Change person from susceptible to latently infected.  added int as argument to DoInfect to record run number: ggilani - 15/10/14
{
	///// This updates a number of things concerning person ai (and their contacts/infectors/places etc.) at time t in thread tn for this run.
	int i;
	unsigned short int ts; //// time step
	double q, x, y; //// q radius squared, x and y coords. q later changed to be quantile of inverse CDF to choose latent period.
	Person* a;

	a = Hosts + ai; //// pointer arithmetic. a = pointer to person. ai = int person index.

	ts = (unsigned short int) (P->TimeStepsPerDay * t);
	a->inf = InfStat_Latent; //// set person a to be infected
	a->infection_time = (unsigned short int) ts; //// record their infection time
	///// Change threaded state variables to reflect new infection status of person a.
	StateT[tn].cumI++;
	StateT[tn].cumItype[a->infect_type % INFECT_TYPE_MASK]++;
	StateT[tn].cumIa[HOST_AGE_GROUP(ai)]++;
	//// calculate radius squared, and increment sum of radii squared.
	x = (Households[a->hh].loc_x - P->LocationInitialInfection[0][0]);
	y = (Households[a->hh].loc_y - P->LocationInitialInfection[0][1]);
	q = x * x + y * y;
	StateT[tn].sumRad2 += q;

	if (q > StateT[tn].maxRad2) StateT[tn].maxRad2 = q; //// update maximum radius squared from seeding infection
	{
		Cells[a->pcell].S--;
		Cells[a->pcell].L++;			//// number of latently infected people increases by one.
		Cells[a->pcell].latent--;		//// pointer to latent in that cell decreased.
		if (a->listpos < Cells[a->pcell].S)
		{
			Cells[a->pcell].susceptible[a->listpos] = Cells[a->pcell].susceptible[Cells[a->pcell].S];
			Hosts[Cells[a->pcell].susceptible[a->listpos]].listpos = a->listpos;
			a->listpos = Cells[a->pcell].S;	//// person a's position with cell.members now equal to number of susceptibles in cell.
			Cells[a->pcell].latent[0] = ai; //// person ai joins front of latent queue.
		}
	}
	StateT[tn].cumI_keyworker[a->keyworker]++;

	if (P->DoLatent)
	{
		i = (int)floor((q = ranf_mt(tn) * CDF_RES));
		q -= ((double)i);
		a->latent_time = (unsigned short int) floor(0.5 + (t - P->LatentPeriod * log(q * P->latent_icdf[i + 1] + (1.0 - q) * P->latent_icdf[i])) * P->TimeStepsPerDay);
	}
	else
		a->latent_time = (unsigned short int) (t * P->TimeStepsPerDay);

	//if (P->DoLatent)	a->latent_time = a->infection_time + ChooseFromICDF(P->latent_icdf, P->LatentPeriod, tn);
	//else			a->latent_time = (unsigned short int) (t * P->TimeStepsPerDay);

	if (P->DoAdUnits)		StateT[tn].cumI_adunit[Mcells[a->mcell].adunit]++;

	if (P->OutputBitmap)
	{
		if ((P->OutputBitmapDetected == 0) || ((P->OutputBitmapDetected == 1) && (Hosts[ai].detected == 1)))
		{
			int ix = ((int)(Households[a->hh].loc_x * P->scalex)) - P->bminx;
			int iy = ((int)(Households[a->hh].loc_y * P->scaley)) - P->bminy;
			if ((ix >= 0) && (ix < P->bwidth) && (iy >= 0) && (iy < P->bheight))
			{
				unsigned j = iy * bmh->width + ix;
				if (j < bmh->imagesize)
				{
#pragma omp atomic
					bmInfected[j]++;
				}
			}
		}
	}
	//added this to record event if flag is set to 1 : ggilani - 10/10/2014
	if (P->DoRecordInfEvents)
	{
		RecordEvent(t, ai, run, 0, tn); //added int as argument to RecordEvent to record run number: ggilani - 15/10/14
	}
	if ((t > 0) && (P->DoOneGen))
	{
		((Latent*)Hosts->stateHandlers[InfStatType_Latent])->GetsWorse(ai, ts, tn, run);
		((InfectiousAlmostSymptomatic*)Hosts->stateHandlers[InfStatType_InfectiousAlmostSymptomatic])->GetsWorse(ai, t, tn, run);
		((Case*)Hosts->stateHandlers[InfStatType_Case])->GetsBetter(ai, t, tn, run);
	}
}
