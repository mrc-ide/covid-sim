#include "InfectiousAlmostSymptomatic.h"
#include "Model.h"
#include "ModelMacros.h"
#include "Rand.h"
#include "Update.h"
#include <math.h>
#include "Case.h"
#include <stdio.h>
#include "Bitmap.h"

extern void DoDetectedCase(int ai, double t, unsigned short int ts, int tn);
extern void DoILI(int ai, int tn);
extern void DoMild(int ai, int tn);

void InfectiousAlmostSymptomatic::GetsWorse(int ai, double t, int tn, int run)
{
	DoCase(ai, t, tn, run);

	Person* a = Hosts + ai;
	a->infectionState = Hosts->stateHandlers[InfStatType_Case];

}

void InfectiousAlmostSymptomatic::DoCase(int ai, double t, unsigned short int ts, int tn) //// makes an infectious (but asymptomatic) person symptomatic. Called in IncubRecoverySweep (and DoInfect if P->DoOneGen)
{
	int j, k, f, j1, j2;
	Person* a;
	int age;

	age = HOST_AGE_GROUP(ai);
	if (age >= NUM_AGE_GROUPS) age = NUM_AGE_GROUPS - 1;
	a = Hosts + ai;
	a->inf = InfStat_Case; //// make person symptomatic and infectious (i.e. a case)
	if (HOST_ABSENT(ai))
	{
		if (a->absent_stop_time < ts + P->usCaseAbsenteeismDelay + P->usCaseAbsenteeismDuration)
			a->absent_stop_time = ts + P->usCaseAbsenteeismDelay + P->usCaseAbsenteeismDuration;
	}
	else if ((P->DoRealSymptWithdrawal) && (P->DoPlaces))
	{
		a->absent_start_time = USHRT_MAX - 1;
		for (j = 0; j < P->PlaceTypeNum; j++)
			if ((a->PlaceLinks[j] >= 0) && (j != P->HotelPlaceType) && (!HOST_ABSENT(ai)) && (P->SymptPlaceTypeWithdrawalProp[j] > 0))
			{
				if ((P->SymptPlaceTypeWithdrawalProp[j] == 1) || (ranf_mt(tn) < P->SymptPlaceTypeWithdrawalProp[j]))
				{
					a->absent_start_time = ts + P->usCaseAbsenteeismDelay;
					a->absent_stop_time = ts + P->usCaseAbsenteeismDelay + P->usCaseAbsenteeismDuration;
					if (P->AbsenteeismPlaceClosure)
					{
						if ((t >= P->PlaceCloseTimeStart) && (!P->DoAdminTriggers) && (!P->DoGlobalTriggers))
							for (j = 0; j < P->PlaceTypeNum; j++)
								if ((j != P->HotelPlaceType) && (a->PlaceLinks[j] >= 0))
									DoPlaceClose(j, a->PlaceLinks[j], ts, tn, 0);
					}
					if ((!HOST_QUARANTINED_PTR(ai)) && (Hosts[ai].PlaceLinks[P->PlaceTypeNoAirNum - 1] >= 0) && (HOST_AGE_YEAR(ai) >= P->CaseAbsentChildAgeCutoff))
						StateT[tn].cumAC++;
					/* This calculates adult absenteeism from work due to care of sick children. Note, children not at school not counted (really this should
					be fixed in population setup by having adult at home all the time for such kids. */
					if ((P->DoHouseholds) && (HOST_AGE_YEAR(ai) < P->CaseAbsentChildAgeCutoff))
					{
						if (!HOST_QUARANTINED_PTR(ai)) StateT[tn].cumACS++;
						if (Hosts[ai].ProbCare < P->CaseAbsentChildPropAdultCarers)
						{
							j1 = Households[Hosts[ai].hh].FirstPerson; j2 = j1 + Households[Hosts[ai].hh].nh;
							f = 0;
							for (int j = j1; (j < j2) && (!f); j++)
								f = ((abs(Hosts[j].inf) != InfStat_Dead) && (HOST_AGE_YEAR(j) >= P->CaseAbsentChildAgeCutoff) && ((Hosts[j].PlaceLinks[P->PlaceTypeNoAirNum - 1] < 0) || (HOST_ABSENT(j)) || (HOST_QUARANTINED_PTR(j))));
							if (!f)
							{
								for (int j = j1; (j < j2) && (!f); j++)
									if ((HOST_AGE_YEAR(j) >= P->CaseAbsentChildAgeCutoff) && (abs(Hosts[j].inf) != InfStat_Dead)) { k = j; f = 1; }
								if (f)
								{
									if (!HOST_ABSENT(k)) Hosts[k].absent_start_time = ts + P->usCaseIsolationDelay;
									Hosts[k].absent_stop_time = ts + P->usCaseIsolationDelay + P->usCaseIsolationDuration;
									StateT[tn].cumAA++;
								}
							}
						}
					}
				}
			}
	}

	//added some case detection code here: ggilani - 03/02/15
	if (Hosts[ai].detected == 1)
		//if ((P->ControlPropCasesId == 1) || (ranf_mt(tn) < P->ControlPropCasesId))
	{
		StateT[tn].cumDC++;
		StateT[tn].cumDC_adunit[Mcells[a->mcell].adunit]++;
		DoDetectedCase(ai, t, ts, tn);
		//add detection time

	}


	if (HOST_TREATED(ai)) Cells[Hosts[ai].pcell].cumTC++;
	StateT[tn].cumC++;
	StateT[tn].cumCa[age]++;
	StateT[tn].cumC_country[Mcells[Hosts[ai].mcell].country]++; //add to cumulative count of cases in that country: ggilani - 12/11/14
	StateT[tn].cumC_keyworker[a->keyworker]++;


	if (P->DoSeverity)
	{
		if (a->Severity_Final == Severity_Mild)
			DoMild(ai, tn);
		else
			DoILI(ai, tn); //// symptomatic cases either mild or ILI at symptom onset. SARI and Critical cases still onset with ILI.
	}
	if (P->DoAdUnits) StateT[tn].cumC_adunit[Mcells[a->mcell].adunit]++;


	

}
