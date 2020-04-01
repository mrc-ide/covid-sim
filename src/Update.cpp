#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "Update.h"
#include "Model.h"
#include "ModelMacros.h"
#include "Param.h"
#include "InfStat.h"
#include "Bitmap.h"
#include "Rand.h"

//adding function to record an event: ggilani - 10/10/2014
void RecordEvent(double, int, int, int, int); //added int as argument to InfectSweep to record run number: ggilani - 15/10/14

unsigned short int ChooseFromICDF(double *, double, int);
int ChooseFinalDiseaseSeverity(int, int);

void DoImmune(int ai)
{
	// This transfers a person straight from susceptible to immune. Used to start a run with a partially immune population.
	person* a;
	int c;
	int x, y;

	a = Hosts + ai;
	if (a->inf == InfStat_Susceptible)
	{
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
		if (P.OutputBitmap)
		{
			x = ((int)(Households[a->hh].loc_x * P.scalex)) - P.bminx;
			y = ((int)(Households[a->hh].loc_y * P.scaley)) - P.bminy;
			if ((x >= 0) && (x < P.bwidth) && (y >= 0) && (y < P.bheight))
			{
				unsigned j = y * bmh->width + x;
				if (j < bmh->imagesize)
				{
#pragma omp atomic
					bmi3[j]++;
				}
			}
		}
	}
}
void DoInfect(int ai, double t, int tn, int run) // Change person from susceptible to latently infected.  added int as argument to DoInfect to record run number: ggilani - 15/10/14
{
	///// This updates a number of things concerning person ai (and their contacts/infectors/places etc.) at time t in thread tn for this run. 
	int i;
	unsigned short int ts; //// time step
	double q, x, y; //// q radius squared, x and y coords. q later changed to be quantile of inverse CDF (I think) to choose latent period.
	person* a;
	a = Hosts + ai; //// pointer arithmetic. a = pointer to person. ai = int person index.

	if (a->inf == InfStat_Susceptible) //// Only change anything if person a/ai uninfected at start of this function.
	{
		ts = (unsigned short int) (P.TimeStepsPerDay * t);
		a->inf = InfStat_Latent; //// set person a to be infected
		a->infection_time = (unsigned short int) ts; //// record their infection time
		///// Change threaded state variables to reflect new infection status of person a.
		StateT[tn].cumI++;
		StateT[tn].cumItype[a->infect_type % INFECT_TYPE_MASK]++;
		StateT[tn].cumIa[HOST_AGE_GROUP(ai)]++;
		//// calculate radius squared, and increment sum of radii squared. 
		x = (Households[a->hh].loc_x - P.LocationInitialInfection[0][0]);
		y = (Households[a->hh].loc_y - P.LocationInitialInfection[0][1]);
		q = x * x + y * y;
		StateT[tn].sumRad2 += q;

		if (q > StateT[tn].maxRad2) StateT[tn].maxRad2 = q; //// update maxium radius squared from seeding infection
		{
			Cells[a->pcell].S--;
			Cells[a->pcell].L++;			//// number of latently infected people increases by one. 
			Cells[a->pcell].latent--;		//// pointer to latent in that cell decreased.
			if (Cells[a->pcell].S > 0)
			{
				Cells[a->pcell].susceptible[a->listpos] = Cells[a->pcell].susceptible[Cells[a->pcell].S];
				Hosts[Cells[a->pcell].susceptible[a->listpos]].listpos = a->listpos;
				a->listpos = Cells[a->pcell].S;	//// person a's position with cell.members now equal to number of susceptibles in cell. 
				Cells[a->pcell].latent[0] = ai; //// person ai joins front of latent queue. 
			}
		}
		StateT[tn].cumI_keyworker[a->keyworker]++;

		if (P.DoLatent)
		{
			i = (int)floor((q = ranf_mt(tn) * CDF_RES));
			q -= ((double)i);
			a->latent_time = (unsigned short int) floor(0.5 + (t - P.LatentPeriod * log(q * P.latent_icdf[i + 1] + (1.0 - q) * P.latent_icdf[i])) * P.TimeStepsPerDay);
		}
		else
			a->latent_time = (unsigned short int) (t * P.TimeStepsPerDay);

		//if (P.DoLatent)	a->latent_time = a->infection_time + ChooseFromICDF(P.latent_icdf, P.LatentPeriod, tn);
		//else			a->latent_time = (unsigned short int) (t * P.TimeStepsPerDay);

		if (P.DoAdUnits)		StateT[tn].cumI_adunit[Mcells[a->mcell].adunit]++;

		if (P.OutputBitmap)
		{
			if ((P.OutputBitmapDetected == 0) || ((P.OutputBitmapDetected == 1) && (Hosts[ai].detected == 1)))
			{
				int ix = ((int)(Households[a->hh].loc_x * P.scalex)) - P.bminx;
				int iy = ((int)(Households[a->hh].loc_y * P.scaley)) - P.bminy;
				if ((x >= 0) && (x < P.bwidth) && (y >= 0) && (y < P.bheight))
				{
					unsigned j = iy * bmh->width + ix;
					if (j < bmh->imagesize)
					{
#pragma omp atomic
						bmi2[j]++;
					}
				}
			}
		}
		//added this to record event if flag is set to 1 : ggilani - 10/10/2014
		if (P.DoRecordInfEvents)
		{
			if (*nEvents < P.MaxInfEvents)
			{
				RecordEvent(t, ai, run, 0, tn); //added int as argument to RecordEvent to record run number: ggilani - 15/10/14
			}
		}
		if ((t > 0) && (P.DoOneGen))
		{
			DoIncub(ai, ts, tn, run);
			DoCase(ai, t, ts, tn);
			DoRecover(ai, run, tn);
		}
	}
}
void RecordEvent(double t, int ai, int run, int type, int tn) //added int as argument to RecordEvent to record run number: ggilani - 15/10/14
{
	/* Function: RecordEvent(t, ai)
	 * Records an infection event in the event log
	 *
	 * Parameters:
	 *	t: time of infection event
	 *	ai: index of infectee
	 *	nEventsPoint: pointer to number of events
	 *
	 * Returns: void
	 *
	 * Author: ggilani, Date: 10/10/2014
	 */
	 //Declare int to store infector's index
	int bi;

	bi = Hosts[ai].infector;

	//Save information to event
#pragma omp critical (inf_event)
	{
		InfEventLog[*nEvents].run = run;
		InfEventLog[*nEvents].type = type;
		InfEventLog[*nEvents].t = t;
		InfEventLog[*nEvents].infectee_ind = ai;
		InfEventLog[*nEvents].infectee_adunit = Mcells[Hosts[ai].mcell].adunit;
		InfEventLog[*nEvents].infectee_x = Households[Hosts[ai].hh].loc_x + P.SpatialBoundingBox[0];
		InfEventLog[*nEvents].infectee_y = Households[Hosts[ai].hh].loc_y + P.SpatialBoundingBox[1];
		InfEventLog[*nEvents].listpos = Hosts[ai].listpos;
		InfEventLog[*nEvents].infectee_cell = Hosts[ai].pcell;
		InfEventLog[*nEvents].thread = tn;
		if (type == 0) //infection event - record time of onset of infector and infector
		{
			InfEventLog[*nEvents].infector_ind = bi;
			if (bi < 0)
			{
				InfEventLog[*nEvents].t_infector = -1;
				InfEventLog[*nEvents].infector_cell = -1;
			}
			else
			{
				InfEventLog[*nEvents].t_infector = (int)(Hosts[bi].infection_time / P.TimeStepsPerDay);
				InfEventLog[*nEvents].infector_cell = Hosts[bi].pcell;
			}
		}
		else if (type == 1) //onset event - record infectee's onset time
		{
			InfEventLog[*nEvents].t_infector = (int)(Hosts[ai].infection_time / P.TimeStepsPerDay);
		}
		else if ((type == 2) || (type == 3)) //recovery or death event - record infectee's onset time
		{
			InfEventLog[*nEvents].t_infector = (int)(Hosts[ai].latent_time / P.TimeStepsPerDay);
		}

		//increment the index of the infection event
		(*nEvents)++;
	}

}

void DoMild(int ai, int tn)
{
	if (P.DoSeverity) //// shouldn't need this but best be careful. 
	{
		person* a = Hosts + ai;
		if (a->Severity_Current == Severity_Asymptomatic)
		{
			a->Severity_Current = Severity_Mild;
			StateT[tn].Mild++;
			StateT[tn].cumMild++;
			if (P.DoAdUnits)
			{
				StateT[tn].Mild_adunit[Mcells[a->mcell].adunit]++;
				StateT[tn].cumMild_adunit[Mcells[a->mcell].adunit]++;
			}
		}
	}
}
void DoILI(int ai, int tn)
{
	if (P.DoSeverity) //// shouldn't need this but best be careful. 
	{
		person* a = Hosts + ai;
		if (a->Severity_Current == Severity_Asymptomatic)
		{
			a->Severity_Current = Severity_ILI;
			StateT[tn].ILI++;
			StateT[tn].cumILI++;
			if (P.DoAdUnits)
			{
				StateT[tn].ILI_adunit[Mcells[a->mcell].adunit]++;
				StateT[tn].cumILI_adunit[Mcells[a->mcell].adunit]++;
			}
		}
	}
}
void DoSARI(int ai, int tn)
{
	if (P.DoSeverity) //// shouldn't need this but best be careful. 
	{
		person* a = Hosts + ai;
		if (a->Severity_Current == Severity_ILI)
		{
			a->Severity_Current = Severity_SARI;
			StateT[tn].ILI--;
			StateT[tn].SARI++;
			StateT[tn].cumSARI++;

			if (P.DoAdUnits)
			{
				StateT[tn].ILI_adunit[Mcells[a->mcell].adunit]--;
				StateT[tn].SARI_adunit[Mcells[a->mcell].adunit]++;
				StateT[tn].cumSARI_adunit[Mcells[a->mcell].adunit]++;
			}
		}
	}
}
void DoCritical(int ai, int tn)
{
	if (P.DoSeverity) //// shouldn't need this but best be careful. 
	{
		person* a = Hosts + ai;
		if (a->Severity_Current == Severity_SARI)
		{
			a->Severity_Current = Severity_Critical;
			StateT[tn].SARI--;
			StateT[tn].Critical++;
			StateT[tn].cumCritical++;

			if (P.DoAdUnits)
			{
				StateT[tn].SARI_adunit[Mcells[a->mcell].adunit]--;
				StateT[tn].Critical_adunit[Mcells[a->mcell].adunit]++;
				StateT[tn].cumCritical_adunit[Mcells[a->mcell].adunit]++;
			}
		}
	}
}
void DoRecoveringFromCritical(int ai, int tn)
{
	//// note function different from DoRecover_FromSeverity. 
	//// DoRecover_FromSeverity assigns people to state Recovered (and bookkeeps accordingly). 
	//// DoRecoveringFromCritical assigns people to intermediate state "recovering from critical condition" (and bookkeeps accordingly). 
	if (P.DoSeverity) //// shouldn't need this but best be careful. 
	{
		person* a = Hosts + ai;
		if (a->Severity_Current == Severity_Critical && (!a->to_die)) //// second condition should be unnecessary but leave in for now. 
		{
			a->Severity_Current = Severity_RecoveringFromCritical;
			StateT[tn].Critical--;
			StateT[tn].CritRecov++;
			StateT[tn].cumCritRecov++;

			if (P.DoAdUnits)
			{
				StateT[tn].Critical_adunit[Mcells[a->mcell].adunit]--;
				StateT[tn].CritRecov_adunit[Mcells[a->mcell].adunit]++;
				StateT[tn].cumCritRecov_adunit[Mcells[a->mcell].adunit]++;
			}
		}
	}
}
void DoDeath_FromCriticalorSARI(int ai, int tn)
{
	//// moved this from DoDeath as I think threading where DoRecover called from IncubRecoverySweep a little weird. May have been something else.

	person* a = Hosts + ai;
	if (P.DoSeverity)
	{

		if (a->Severity_Current == Severity_Critical)  //// second condition should be unnecessary but leave in for now. 
		{
			StateT[tn].Critical--;
			if (P.DoAdUnits)	StateT[tn].Critical_adunit[Mcells[a->mcell].adunit]--;
			//// change current status so that flags work. 
		}
		else if (a->Severity_Current == Severity_SARI)
		{
			StateT[tn].SARI--;
			if (P.DoAdUnits)	StateT[tn].SARI_adunit[Mcells[a->mcell].adunit]--;
			//// change current status so that flags work. 
		}
		a->Severity_Current = Severity_Dead;
	}
}
void DoRecover_FromSeverity(int ai, int tn)
{
	//// note function different from DoRecoveringFromCritical. 
	//// DoRecover_FromSeverity assigns people to state Recovered (and bookkeeps accordingly). 
	//// DoRecoveringFromCritical assigns people to intermediate state "recovering from critical condition" (and bookkeeps accordingly). 

	//// moved this from DoRecover as I think threading where DoRecover called from IncubRecoverySweep a little weird. Talk to Gemma/Neil.
	person* a = Hosts + ai;

	if (P.DoSeverity)
		if (a->inf == InfStat_InfectiousAsymptomaticNotCase || a->inf == InfStat_Case) ///// i.e same condition in DoRecover (make sure you don't recover people twice). 
		{
			//StateT[tn].R++; ///// Don't think you need this. variables .S, .I, .L, .R, .D (in popvar/StateT) aren't used during runtime (threads), these states used in unthreaded State, and kept track of during runtime through .S, .I etc. in cell. Same thing with death .D value. 

			if (a->Severity_Current == Severity_Mild)
			{
				StateT[tn].Mild--;
				if (P.DoAdUnits) StateT[tn].Mild_adunit[Mcells[a->mcell].adunit]--;
			}
			else if (a->Severity_Current == Severity_ILI)
			{
				StateT[tn].ILI--;
				if (P.DoAdUnits) StateT[tn].ILI_adunit[Mcells[a->mcell].adunit]--;
			}
			else if (a->Severity_Current == Severity_SARI)
			{
				StateT[tn].SARI--;
				if (P.DoAdUnits) StateT[tn].SARI_adunit[Mcells[a->mcell].adunit]--;
			}
			else if (a->Severity_Current == Severity_RecoveringFromCritical)
			{
				StateT[tn].CritRecov--; //// decrement CritRecov, not critical. 
				if (P.DoAdUnits) StateT[tn].CritRecov_adunit[Mcells[a->mcell].adunit]--;
			}
			//// change current status so that flags work. 
			a->Severity_Current = Severity_Recovered;
		}
}


void DoIncub(int ai, unsigned short int ts, int tn, int run)
{
	person* a;
	double q;
	int age;

	age = HOST_AGE_GROUP(ai);
	if (age >= NUM_AGE_GROUPS) age = NUM_AGE_GROUPS - 1;

	a = Hosts + ai;
	if (a->inf == InfStat_Latent)
	{
		if (P.InfectiousnessSD == 0)	a->infectiousness = (float)P.AgeInfectiousness[age];
		else							a->infectiousness = (float)(P.AgeInfectiousness[age] * gen_gamma_mt(P.InfectiousnessGamA, P.InfectiousnessGamR, tn));

		q = P.ProportionSymptomatic[age]
			* (HOST_TREATED(ai) ? (1 - P.TreatSympDrop) : 1)
			* (HOST_VACCED(ai) ? (1 - P.VaccSympDrop) : 1);

		if (ranf_mt(tn) < q)
		{
			a->inf = InfStat_InfectiousAlmostSymptomatic;
			a->infectiousness = (float)(-P.SymptInfectiousness * a->infectiousness);
		}
		else
			a->inf = InfStat_InfectiousAsymptomaticNotCase;

		if (!P.DoSeverity || a->inf == InfStat_InfectiousAsymptomaticNotCase) //// if not doing severity or if person asymptomatic. 
		{
			if (P.DoInfectiousnessProfile)	a->recovery_time = a->latent_time + (unsigned short int) (P.InfectiousPeriod * P.TimeStepsPerDay);
			else							a->recovery_time = a->latent_time + ChooseFromICDF(P.infectious_icdf, P.InfectiousPeriod, tn);
		}
		else
		{
			//// choose final disease severity (either mild, ILI, SARI, Critical, not asymptomatic as covered above) by age
			a->Severity_Final = ChooseFinalDiseaseSeverity(age, tn);

			/// choose outcome recovery or death
			if (((a->Severity_Final == Severity_Critical) && (ranf_mt(tn) < P.CFR_Critical_ByAge[age])) || ((a->Severity_Final == Severity_SARI) && (ranf_mt(tn) < P.CFR_SARI_ByAge[age])))
				a->to_die = 1;

			//// re-define Final severity to Severity_Critical if person dies. 
			// #### Neil - don't want this ///
//			if (a->to_die) a->Severity_Final = Severity_Critical;

			//// choose events and event times
			if (a->Severity_Final == Severity_Mild)
				a->recovery_time = a->latent_time + ChooseFromICDF(P.MildToRecovery_icdf, P.Mean_MildToRecovery, tn);
			else if (a->Severity_Final == Severity_Critical)
			{
				a->SARI_time = a->latent_time + ChooseFromICDF(P.ILIToSARI_icdf, P.Mean_ILIToSARI, tn);
				a->Critical_time = a->SARI_time + ChooseFromICDF(P.SARIToCritical_icdf, P.Mean_SARIToCritical, tn);
				if (a->to_die)
					a->recovery_time = a->Critical_time + ChooseFromICDF(P.CriticalToDeath_icdf, P.Mean_CriticalToDeath, tn);
				else
				{
					a->RecoveringFromCritical_time = a->Critical_time + ChooseFromICDF(P.CriticalToCritRecov_icdf, P.Mean_CriticalToCritRecov, tn);
					a->recovery_time = a->RecoveringFromCritical_time + ChooseFromICDF(P.CritRecovToRecov_icdf, P.Mean_CritRecovToRecov, tn);
				}
			}
			else if (a->Severity_Final == Severity_SARI)
			{
				a->SARI_time = a->latent_time + ChooseFromICDF(P.ILIToSARI_icdf, P.Mean_ILIToSARI, tn);
				if (a->to_die)
					a->recovery_time = a->SARI_time + ChooseFromICDF(P.CriticalToDeath_icdf, P.Mean_CriticalToDeath, tn);  // update with its own time
				else
					a->recovery_time = a->SARI_time + ChooseFromICDF(P.SARIToRecovery_icdf, P.Mean_SARIToRecovery, tn);
			}
			else /*i.e. if Severity_Final == Severity_ILI*/
				a->recovery_time = a->latent_time + ChooseFromICDF(P.ILIToRecovery_icdf, P.Mean_ILIToRecovery, tn);
		}

		if ((P.ControlPropCasesId == 1) || (ranf_mt(tn) < P.ControlPropCasesId))
		{
			Hosts[ai].detected = 1;
		}


		//// update pointers
		Cells[a->pcell].L--; //// one fewer person latently infected. 
		Cells[a->pcell].infected--; //// first infected person is now one index earlier in array. 
		Cells[a->pcell].I++; //// one more infectious person. 
		if (Cells[a->pcell].L > 0)
		{
			Cells[a->pcell].susceptible[a->listpos] = Cells[a->pcell].latent[Cells[a->pcell].L]; //// reset pointers.
			Hosts[Cells[a->pcell].susceptible[a->listpos]].listpos = a->listpos;
			a->listpos = Cells[a->pcell].S + Cells[a->pcell].L; //// change person a's listpos, which will now refer to their position among infectious people, not latent. 
			Cells[a->pcell].infected[0] = ai; //// this person is now first infectious person in the array? I think because the pointer was moved back one so now that bit of memorty needs to refer to person ai. Alternative would be to move everyone back one which would take longer. 
		}
		////added this to record event if flag is set to 1 and if host isn't initial seed, i.e. if Hosts[ai].infector>=0: ggilani - 10/10/2014
		//if(P.DoRecordInfEvents)
		//	{
		//		if(*nEvents<P.MaxInfEvents)
		//		{
		//			RecordEvent(((double)a->latent_time)/P.TimeStepsPerDay,ai,run,1); //added int as argument to RecordEvent to record run number: ggilani - 15/10/14
		//		}
		//	}
	}
}
void DoDetectedCase(int ai, double t, unsigned short int ts, int tn)
{
	int i, j, k, f, j1, j2;
	person* a = Hosts + ai;

	//// Increment triggers (Based on numbers of detected cases) for interventions. Used in TreatSweep function when not doing Global or Admin triggers. And not when doing ICU triggers.
	if (Mcells[a->mcell].treat_trig				< USHRT_MAX - 1) Mcells[a->mcell].treat_trig++;
	if (Mcells[a->mcell].vacc_trig				< USHRT_MAX - 1) Mcells[a->mcell].vacc_trig++;
	if (Mcells[a->mcell].move_trig				< USHRT_MAX - 1) Mcells[a->mcell].move_trig++;
	if (Mcells[a->mcell].socdist_trig			< USHRT_MAX - 1) Mcells[a->mcell].socdist_trig++;
	if (Mcells[a->mcell].keyworkerproph_trig	< USHRT_MAX - 1) Mcells[a->mcell].keyworkerproph_trig++;
#ifndef ABSENTEEISM_PLACE_CLOSURE
#ifdef PLACE_CLOSE_ROUND_HOUSEHOLD
	if (Mcells[a->mcell].place_trig < USHRT_MAX - 1) Mcells[a->mcell].place_trig++;
#endif
	if ((t >= P.PlaceCloseTimeStart) && (!P.DoAdminTriggers) && (!P.DoGlobalTriggers))
		for (j = 0; j < P.PlaceTypeNum; j++)
			if ((j != HOTEL_PLACE_TYPE) && (a->PlaceLinks[j] >= 0))
			{
				DoPlaceClose(j, a->PlaceLinks[j], ts, tn, 0);
#ifndef PLACE_CLOSE_ROUND_HOUSEHOLD
				if (Mcells[Places[j][a->PlaceLinks[j]].mcell].place_trig < USHRT_MAX - 1)
				{
#pragma omp critical (place_trig)
					Mcells[Places[j][a->PlaceLinks[j]].mcell].place_trig++;
				}
#endif
			}
#endif
	if (t >= P.TreatTimeStart)
		if ((P.TreatPropCases == 1) || (ranf_mt(tn) < P.TreatPropCases))
		{
			DoTreatCase(ai, ts, tn);
			if (P.DoHouseholds)
			{
				if ((t < P.TreatTimeStart + P.TreatHouseholdsDuration) && ((P.TreatPropCaseHouseholds == 1) || (ranf_mt(tn) < P.TreatPropCaseHouseholds)))
				{
					j1 = Households[Hosts[ai].hh].FirstPerson; j2 = j1 + Households[Hosts[ai].hh].nh;
					for (j = j1; j < j2; j++)
						if (!HOST_TO_BE_TREATED(j)) DoProph(j, ts, tn);
				}
			}
			if (P.DoPlaces)
			{
				if (t < P.TreatTimeStart + P.TreatPlaceGeogDuration)
					for (j = 0; j < P.PlaceTypeNum; j++)
						if (a->PlaceLinks[j] >= 0)
						{
							if (P.DoPlaceGroupTreat)
							{
								if ((P.TreatPlaceProbCaseId[j] == 1) || (ranf_mt(tn) < P.TreatPlaceProbCaseId[j]))
								{
									StateT[tn].p_queue[j][StateT[tn].np_queue[j]] = a->PlaceLinks[j];
									StateT[tn].pg_queue[j][StateT[tn].np_queue[j]++] = a->PlaceGroupLinks[j];
								}
							}
							else
							{
								f = 0;
#pragma omp critical (starttreat)
								if (!Places[j][a->PlaceLinks[j]].treat) f = Places[j][a->PlaceLinks[j]].treat = 1;
								if (f)
								{
									if ((P.TreatPlaceProbCaseId[j] == 1) || (ranf_mt(tn) < P.TreatPlaceProbCaseId[j]))
										StateT[tn].p_queue[j][StateT[tn].np_queue[j]++] = a->PlaceLinks[j];
									else
										Places[j][a->PlaceLinks[j]].treat = 0;
								}
							}
						}
			}
		}
	if (P.DoHouseholds)
	{
		if ((!P.DoMassVacc) && (t >= P.VaccTimeStart) && (State.cumV < P.VaccMaxCourses))
			if ((t < P.VaccTimeStart + P.VaccHouseholdsDuration) && ((P.VaccPropCaseHouseholds == 1) || (ranf_mt(tn) < P.VaccPropCaseHouseholds)))
			{
				j1 = Households[Hosts[ai].hh].FirstPerson; j2 = j1 + Households[Hosts[ai].hh].nh;
				for (j = j1; j < j2; j++) DoVacc(j, ts);
			}

		//// Giant compound if statement. If doing delays by admin unit, then window of HQuarantine dependent on admin unit-specific duration. This if statement ensures that this timepoint within window, regardless of how window defined.
		if ((P.DoInterventionDelaysByAdUnit &&
			(t >= AdUnits[Mcells[a->mcell].adunit].HQuarantineTimeStart		&&	(t < AdUnits[Mcells[a->mcell].adunit].HQuarantineTimeStart + AdUnits[Mcells[a->mcell].adunit].HQuarantineDuration)))		||
			(t >= AdUnits[Mcells[a->mcell].adunit].HQuarantineTimeStart		&&	(t < AdUnits[Mcells[a->mcell].adunit].HQuarantineTimeStart + P.HQuarantinePolicyDuration))									)
		{
			j1 = Households[Hosts[ai].hh].FirstPerson; j2 = j1 + Households[Hosts[ai].hh].nh;
			if ((!HOST_TO_BE_QUARANTINED(j1)) || (P.DoHQretrigger))
			{
				Hosts[j1].quar_start_time = ts + ((unsigned short int) (P.TimeStepsPerDay * P.HQuarantineHouseDelay));
				k = (ranf_mt(tn) < P.HQuarantinePropHouseCompliant) ? 1 : 0;
				if (k) StateT[tn].cumHQ++;

				Hosts[j1].quar_comply = ((k == 0) ? 0 : ((ranf_mt(tn) < P.HQuarantinePropIndivCompliant) ? 1 : 0));
				if ((Hosts[j1].quar_comply) && (!HOST_ABSENT(j1)))
				{
					if (HOST_AGE_YEAR(j1) >= P.CaseAbsentChildAgeCutoff)
					{
						if (Hosts[j1].PlaceLinks[NUM_PLACE_TYPES_NOAIR - 1] >= 0) StateT[tn].cumAH++;
					}
					else		StateT[tn].cumACS++;
				}
				for (j = j1 + 1; j < j2; j++)
				{
					Hosts[j].quar_start_time = Hosts[j1].quar_start_time;
					Hosts[j].quar_comply = ((k == 0) ? 0 : ((ranf_mt(tn) < P.HQuarantinePropIndivCompliant) ? 1 : 0));
					if ((Hosts[j].quar_comply) && (!HOST_ABSENT(j)))
					{
						if (HOST_AGE_YEAR(j) >= P.CaseAbsentChildAgeCutoff)
						{
							if (Hosts[j].PlaceLinks[NUM_PLACE_TYPES_NOAIR - 1] >= 0) StateT[tn].cumAH++;
						}
						else	StateT[tn].cumACS++;
					}
				}
			}
		}
	}

	//// Giant compound if statement. If doing delays by admin unit, then window of case isolation dependent on admin unit-specific duration. This if statement ensures that this timepoint within window, regardless of how window defined.
	if ((P.DoInterventionDelaysByAdUnit && 
		(t >= AdUnits[Mcells[a->mcell].adunit].CaseIsolationTimeStart && (t < AdUnits[Mcells[a->mcell].adunit].CaseIsolationTimeStart + AdUnits[Mcells[a->mcell].adunit].CaseIsolationDuration)))	||
		(t >= AdUnits[Mcells[a->mcell].adunit].CaseIsolationTimeStart && (t < AdUnits[Mcells[a->mcell].adunit].CaseIsolationTimeStart + P.CaseIsolationPolicyDuration))								)
		if ((P.CaseIsolationProp == 1) || (ranf_mt(tn) < P.CaseIsolationProp))
		{
			Hosts[ai].isolation_start_time = ts; //// set isolation start time.
			if (HOST_ABSENT(ai))
			{
				if (a->absent_stop_time < ts + P.usCaseAbsenteeismDelay + P.usCaseIsolationDuration) //// ensure that absent_stop_time is at least now + CaseIsolationDuraton
					a->absent_stop_time = ts + P.usCaseAbsenteeismDelay + P.usCaseIsolationDuration;
			}
			else if (P.DoRealSymptWithdrawal) /* This calculates adult absenteeism from work due to care of isolated children.  */
			{
				Hosts[ai].absent_start_time = ts + P.usCaseIsolationDelay;
				Hosts[ai].absent_stop_time	= ts + P.usCaseIsolationDelay + P.usCaseIsolationDuration;
				if (P.DoPlaces)
				{
					if ((!HOST_QUARANTINED(ai)) && (Hosts[ai].PlaceLinks[NUM_PLACE_TYPES_NOAIR - 1] >= 0) && (HOST_AGE_YEAR(ai) >= P.CaseAbsentChildAgeCutoff))
						StateT[tn].cumAC++;
				}
				if ((P.DoHouseholds) && (P.DoPlaces) && (HOST_AGE_YEAR(ai) < P.CaseAbsentChildAgeCutoff))
				{
					if (!HOST_QUARANTINED(ai)) StateT[tn].cumACS++;
					if ((P.CaseAbsentChildPropAdultCarers == 1) || (ranf_mt(tn) < P.CaseAbsentChildPropAdultCarers))
					{
						j1 = Households[Hosts[ai].hh].FirstPerson; j2 = j1 + Households[Hosts[ai].hh].nh;
						f = 0;

						/*loop through household members. 
						Ask whether a) household member alive; AND 
									b) household member is not a child requiring adult to stay home; AND 
									c) household member doesn't have links to an office. 
						If any household member satisfies all of these conditions, don't do following code block.
						*/
						for (j = j1; (j < j2) && (!f); j++) 
							f = ((abs(Hosts[j].inf) != InfStat_Dead) && (HOST_AGE_YEAR(j) >= P.CaseAbsentChildAgeCutoff) && (Hosts[j].PlaceLinks[NUM_PLACE_TYPES_NOAIR - 1] < 0));
						if (!f) //// so if either a) a household member is dead; b) a household member is a child requiring adult to stay home; c) a household member has links to office.
						{
							for (j = j1; (j < j2) && (!f); j++) //// loop through people again, now checking whether household member already absent or quarantined. If they are, don't do following block.
								f = ((abs(Hosts[j].inf) != InfStat_Dead) && (HOST_AGE_YEAR(j) >= P.CaseAbsentChildAgeCutoff) && ((HOST_ABSENT(j)) || (HOST_QUARANTINED(j))));
							if (!f) //// so if no household members dead, or if one hh member is a child requiring adult supervision, and neither absent nor quarantined.
							{
								k = -1; /// don't think this is necessary...
								for (j = j1; (j < j2) & (!f); j++) /// loop again, checking whehter household members not children needing supervision and are alive.
									if ((HOST_AGE_YEAR(j) >= P.CaseAbsentChildAgeCutoff) && (abs(Hosts[j].inf) != InfStat_Dead)) { k = j; f = 1; }
								if (f) //// so finally, if at least one member of household is alive and does not need supervision by an adult, amend absent start and stop times
								{
									Hosts[k].absent_start_time = ts + P.usCaseIsolationDelay;
									Hosts[k].absent_stop_time = ts + P.usCaseIsolationDelay + P.usCaseIsolationDuration;
									StateT[tn].cumAA++;
								}
							}
						}
					}
				}
			}
		}

	if (P.DoDigitalContactTracing)  
		if ((t >= AdUnits[Mcells[a->mcell].adunit].DigitalContactTracingTimeStart) && (t < AdUnits[Mcells[a->mcell].adunit].DigitalContactTracingTimeStart + P.DigitalContactTracingPolicyDuration))
			if (Hosts[ai].digitalContactTracingUser) //if host is a digital contact tracing app user, need to add their contacts to a list of people to be notified
				for (i = 0; i < Hosts[ai].ncontacts; i++)
					if ((P.ProportionDigitalContactsIsolate == 1) || (ranf_mt(tn) < P.ProportionDigitalContactsIsolate)) //add contacts to queue of digital contacts, assuming a certain proportion self-isolate. Note: need to add contacts to the queue related to their own admin unit, not admin unit of infector
						StateT[tn].dct_queue[Mcells[Hosts[Hosts[ai].contacts[i]].mcell].adunit][StateT[tn].ndct_queue[Mcells[Hosts[Hosts[ai].contacts[i]].mcell].adunit]] = Hosts[ai].contacts[i];
}

void DoCase(int ai, double t, unsigned short int ts, int tn) //// makes an infectious (but asymptomatic) person symptomatic. Called in IncubRecoverySweep (and DoInfect if P.DoOneGen)
{
	int j, k, f, j1, j2;
	person* a;
	int age;

	age = HOST_AGE_GROUP(ai);
	if (age >= NUM_AGE_GROUPS) age = NUM_AGE_GROUPS - 1;
	a = Hosts + ai;
	if (a->inf == InfStat_InfectiousAlmostSymptomatic) //// if person latent/asymptomatically infected, but infectious
	{
		a->inf = InfStat_Case; //// make person symptomatic and infectious (i.e. a case)
		if (HOST_ABSENT(ai))
		{
			if (a->absent_stop_time < ts + P.usCaseAbsenteeismDelay + P.usCaseAbsenteeismDuration)
				a->absent_stop_time = ts + P.usCaseAbsenteeismDelay + P.usCaseAbsenteeismDuration;
		}
		else
		{
			a->absent_start_time = USHRT_MAX - 1;
			if (P.DoPlaces)
				for (j = 0; j < P.PlaceTypeNum; j++)
					if ((a->PlaceLinks[j] >= 0) && (j != HOTEL_PLACE_TYPE) && (!HOST_ABSENT(ai)) && (P.SymptPlaceTypeWithdrawalProp[j] > 0))
					{
						if ((P.SymptPlaceTypeWithdrawalProp[j] == 1) || (ranf_mt(tn) < P.SymptPlaceTypeWithdrawalProp[j]))
						{
							a->absent_start_time = ts + P.usCaseAbsenteeismDelay;
							a->absent_stop_time = ts + P.usCaseAbsenteeismDelay + P.usCaseAbsenteeismDuration;
#ifdef ABSENTEEISM_PLACE_CLOSURE
							if ((t >= P.PlaceCloseTimeStart) && (!P.DoAdminTriggers) && (!P.DoGlobalTriggers))
								for (j = 0; j < P.PlaceTypeNum; j++)
									if ((j != HOTEL_PLACE_TYPE) && (a->PlaceLinks[j] >= 0))
											DoPlaceClose(j, a->PlaceLinks[j], ts, tn, 0);
#endif
							if ((!HOST_QUARANTINED(ai)) && (Hosts[ai].PlaceLinks[NUM_PLACE_TYPES_NOAIR - 1] >= 0) && (HOST_AGE_YEAR(ai) >= P.CaseAbsentChildAgeCutoff))
								StateT[tn].cumAC++;
							/* This calculates adult absenteeism from work due to care of sick children. Note, children not at school not counted (really this should
							be fixed in population setup by having adult at home all the time for such kids. */
							if ((P.DoHouseholds) && (HOST_AGE_YEAR(ai) < P.CaseAbsentChildAgeCutoff))
							{
								if (!HOST_QUARANTINED(ai)) StateT[tn].cumACS++;
								if ((P.CaseAbsentChildPropAdultCarers == 1) || (ranf_mt(tn) < P.CaseAbsentChildPropAdultCarers))
								{
									j1 = Households[Hosts[ai].hh].FirstPerson; j2 = j1 + Households[Hosts[ai].hh].nh;
									f = 0;
									for (int j = j1; (j < j2) && (!f); j++)
										f = ((abs(Hosts[j].inf) != InfStat_Dead) && (HOST_AGE_YEAR(j) >= P.CaseAbsentChildAgeCutoff) && (Hosts[j].PlaceLinks[NUM_PLACE_TYPES_NOAIR - 1] < 0));
									if (!f)
									{
										for (int j = j1; (j < j2) && (!f); j++)
											f = ((abs(Hosts[j].inf) != InfStat_Dead) && (HOST_AGE_YEAR(j) >= P.CaseAbsentChildAgeCutoff) && ((HOST_ABSENT(j)) || (HOST_QUARANTINED(j))));
										if (!f)
										{
											k = -1;
											for (int j = j1; (j < j2) && (!f); j++)
												if ((HOST_AGE_YEAR(j) >= P.CaseAbsentChildAgeCutoff) && (abs(Hosts[j].inf) != InfStat_Dead)) { k = j; f = 1; }
											if (f)
											{
												Hosts[k].absent_start_time = ts + P.usCaseIsolationDelay;
												Hosts[k].absent_stop_time = ts + P.usCaseIsolationDelay + P.usCaseIsolationDuration;
												StateT[tn].cumAA++;
											}
										}
									}
								}
							}
						}
					}
		}

		//added some case detection code here: ggilani - 03/02/15
		if (Hosts[ai].detected == 1)
			//if ((P.ControlPropCasesId == 1) || (ranf_mt(tn) < P.ControlPropCasesId))
		{
			StateT[tn].cumDC++;
			StateT[tn].cumDC_adunit[Mcells[a->mcell].adunit]++;
			DoDetectedCase(ai, t, ts, tn);
		}

		if (HOST_TREATED(ai)) Cells[Hosts[ai].pcell].cumTC++;
		StateT[tn].cumC++;
		StateT[tn].cumCa[age]++;
		StateT[tn].cumC_country[Mcells[Hosts[ai].mcell].country]++; //add to cumulative count of cases in that country: ggilani - 12/11/14
		StateT[tn].cumC_keyworker[a->keyworker]++;


		if (P.DoSeverity)
		{
			if (a->Severity_Final == Severity_Mild)
				DoMild(ai, tn);
			else
				DoILI(ai, tn); //// symptomatic cases either mild or ILI at symptom onset. 
		}
		if (P.DoAdUnits) StateT[tn].cumC_adunit[Mcells[a->mcell].adunit]++;
	}
}

void DoFalseCase(int ai, double t, unsigned short int ts, int tn)
{
	person* a;

	/* Arguably adult absenteeism to take care of sick kids could be included here, but then output absenteeism would not be 'excess' absenteeism */
	a = Hosts + ai;
	if ((P.ControlPropCasesId == 1) || (ranf_mt(tn) < P.ControlPropCasesId))
	{
		if ((!P.DoEarlyCaseDiagnosis) || (State.cumDC >= P.PreControlClusterIdCaseThreshold)) StateT[tn].cumDC++;
		DoDetectedCase(ai, t, ts, tn);
	}
	StateT[tn].cumFC++;
}

void DoRecover(int ai, int run, int tn)
{
	int i, j, x, y;
	person* a;

	a = Hosts + ai;
	if (a->inf == InfStat_InfectiousAsymptomaticNotCase || a->inf == InfStat_Case)
	{
		i = a->listpos;
		Cells[a->pcell].I--;
		if (P.DoSIS)
		{
			if (Cells[a->pcell].I > 0)
			{
				Cells[a->pcell].susceptible[i] = Cells[a->pcell].infected[0];
				Hosts[Cells[a->pcell].susceptible[i]].listpos = i;
				if (Cells[a->pcell].L > 0)
				{
					Cells[a->pcell].latent[Cells[a->pcell].L] = Cells[a->pcell].latent[0];
					Hosts[Cells[a->pcell].latent[Cells[a->pcell].L]].listpos = Cells[a->pcell].S + Cells[a->pcell].L;
				}
			}
			else if (Cells[a->pcell].L > 0)
			{
				Cells[a->pcell].susceptible[i] = Cells[a->pcell].latent[0];
				Hosts[Cells[a->pcell].susceptible[i]].listpos = i;
			}
			Cells[a->pcell].susceptible[Cells[a->pcell].S] = ai;
			a->listpos = Cells[a->pcell].S;
			Cells[a->pcell].S++;
			Cells[a->pcell].latent++;
			Cells[a->pcell].infected++;
			a->susc = (float)(a->susc * P.SuscReductionFactorPerInfection);
			a->inf = InfStat_Susceptible;
			a->infector = -1;
			if (P.OutputBitmap)
			{
				if ((P.OutputBitmapDetected == 0) || ((P.OutputBitmapDetected == 1) && (Hosts[ai].detected == 1)))
				{

					x = ((int)(Households[a->hh].loc_x * P.scalex)) - P.bminx;
					y = ((int)(Households[a->hh].loc_y * P.scaley)) - P.bminy;
					if ((x >= 0) && (x < P.bwidth) && (y >= 0) && (y < P.bheight))
					{
						unsigned j = y * bmh->width + x;
						if (j < bmh->imagesize)
						{
#pragma omp atomic
							bmi2[j]--;
						}
					}
				}
			}

		}
		else
		{
			Cells[a->pcell].R++;
			j = Cells[a->pcell].S + Cells[a->pcell].L + Cells[a->pcell].I;
			if (i < Cells[a->pcell].S + Cells[a->pcell].L + Cells[a->pcell].I)
			{
				Cells[a->pcell].susceptible[i] = Cells[a->pcell].susceptible[j];
				Hosts[Cells[a->pcell].susceptible[i]].listpos = i;
				a->listpos = j;
				Cells[a->pcell].susceptible[j] = ai;
			}
			a->inf = InfStat_Recovered * a->inf / abs(a->inf);

			if (P.OutputBitmap)
			{
				if ((P.OutputBitmapDetected == 0) || ((P.OutputBitmapDetected == 1) && (Hosts[ai].detected == 1)))
				{
					x = ((int)(Households[a->hh].loc_x * P.scalex)) - P.bminx;
					y = ((int)(Households[a->hh].loc_y * P.scaley)) - P.bminy;
					if ((x >= 0) && (x < P.bwidth) && (y >= 0) && (y < P.bheight))
					{
						unsigned j = y * bmh->width + x;
						if (j < bmh->imagesize)
						{
#pragma omp atomic
							bmi3[j]++;
#pragma omp atomic
							bmi2[j]--;
						}
					}
				}
			}
		}
	}
	//else
	//fprintf(stderr, "\n ### %i %i  \n", ai, a->inf);
}

void DoDeath(int ai, int tn, int run)
{
	int i, x, y;
	person* a = Hosts + ai;

	if ((a->inf == InfStat_InfectiousAsymptomaticNotCase || a->inf == InfStat_Case)) //added infection status of 6 as well, as this also needs to deal with 'dead' infectious individuals: ggilani 25/10/14
	{
		a->inf = InfStat_Dead * a->inf / abs(a->inf);
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
		if (P.DoAdUnits) StateT[tn].cumD_adunit[Mcells[a->mcell].adunit]++;
		if (P.OutputBitmap)
		{
			if ((P.OutputBitmapDetected == 0) || ((P.OutputBitmapDetected == 1) && (Hosts[ai].detected == 1)))
			{
				x = ((int)(Households[a->hh].loc_x * P.scalex)) - P.bminx;
				y = ((int)(Households[a->hh].loc_y * P.scaley)) - P.bminy;
				if ((x >= 0) && (x < P.bwidth) && (y >= 0) && (y < P.bheight))
				{
					unsigned j = y * bmh->width + x;
					if (j < bmh->imagesize)
					{
#pragma omp atomic
						bmi3[j]++;
#pragma omp atomic
						bmi2[j]--;
					}
				}
			}
		}
	}
}

void DoTreatCase(int ai, unsigned short int ts, int tn)
{
	int x, y;

	if (State.cumT < P.TreatMaxCourses)
	{
#ifdef NO_TREAT_PROPH_CASES
		if (!HOST_TO_BE_TREATED(ai))
#endif
		{
			Hosts[ai].treat_start_time = ts + ((unsigned short int) (P.TimeStepsPerDay * P.TreatDelayMean));
			Hosts[ai].treat_stop_time = ts + ((unsigned short int) (P.TimeStepsPerDay * (P.TreatDelayMean + P.TreatCaseCourseLength)));
			StateT[tn].cumT++;
			if ((abs(Hosts[ai].inf) > InfStat_Susceptible) && (Hosts[ai].inf != InfStat_Dead_WasAsymp)) Cells[Hosts[ai].pcell].cumTC++;
			StateT[tn].cumT_keyworker[Hosts[ai].keyworker]++;
			if ((++Hosts[ai].num_treats) < 2) StateT[tn].cumUT++;
			Cells[Hosts[ai].pcell].tot_treat++;
			if (P.DoAdUnits) StateT[tn].cumT_adunit[Mcells[Hosts[ai].mcell].adunit]++;
			if (P.OutputBitmap)
			{
				x = ((int)(Households[Hosts[ai].hh].loc_x * P.scalex)) - P.bminx;
				y = ((int)(Households[Hosts[ai].hh].loc_y * P.scaley)) - P.bminy;
				if ((x >= 0) && (x < P.bwidth) && (y >= 0) && (y < P.bheight))
				{
					unsigned j = y * bmh->width + x;
					if (j < bmh->imagesize)
					{
#pragma omp atomic
						bmi4[j]++;
					}
				}
			}
		}
	}
}

void DoProph(int ai, unsigned short int ts, int tn)
{
	//// almost identical to DoProphNoDelay, except unsurprisingly this function includes delay between timestep and start of treatment. Also increments StateT[tn].cumT_keyworker by 1 every time. 
	int x, y;

	if (State.cumT < P.TreatMaxCourses)
	{
		Hosts[ai].treat_start_time = ts + ((unsigned short int) (P.TimeStepsPerDay * P.TreatDelayMean));
		Hosts[ai].treat_stop_time = ts + ((unsigned short int) (P.TimeStepsPerDay * (P.TreatDelayMean + P.TreatProphCourseLength)));
		StateT[tn].cumT++;
		StateT[tn].cumT_keyworker[Hosts[ai].keyworker]++;
		if ((++Hosts[ai].num_treats) < 2) StateT[tn].cumUT++;
		if (P.DoAdUnits)	StateT[tn].cumT_adunit[Mcells[Hosts[ai].mcell].adunit]++;
#pragma omp critical (tot_treat)
		Cells[Hosts[ai].pcell].tot_treat++;
		if (P.OutputBitmap)
		{
			x = ((int)(Households[Hosts[ai].hh].loc_x * P.scalex)) - P.bminx;
			y = ((int)(Households[Hosts[ai].hh].loc_y * P.scaley)) - P.bminy;
			if ((x >= 0) && (x < P.bwidth) && (y >= 0) && (y < P.bheight))
			{
				unsigned j = y * bmh->width + x;
				if (j < bmh->imagesize)
				{
#pragma omp atomic
					bmi4[j]++;
				}
			}
		}
	}
}

void DoProphNoDelay(int ai, unsigned short int ts, int tn, int nc)
{
	int x, y;

	if (State.cumT < P.TreatMaxCourses)
	{
		Hosts[ai].treat_start_time = ts;
		Hosts[ai].treat_stop_time = ts + ((unsigned short int) (P.TimeStepsPerDay * P.TreatProphCourseLength * nc));
		StateT[tn].cumT += nc;
		StateT[tn].cumT_keyworker[Hosts[ai].keyworker] += nc;
		if ((++Hosts[ai].num_treats) < 2) StateT[tn].cumUT++;
		if (P.DoAdUnits) StateT[tn].cumT_adunit[Mcells[Hosts[ai].mcell].adunit] += nc;
#pragma omp critical (tot_treat)
		Cells[Hosts[ai].pcell].tot_treat++;
		if (P.OutputBitmap)
		{
			x = ((int)(Households[Hosts[ai].hh].loc_x * P.scalex)) - P.bminx;
			y = ((int)(Households[Hosts[ai].hh].loc_y * P.scaley)) - P.bminy;
			if ((x >= 0) && (x < P.bwidth) && (y >= 0) && (y < P.bheight))
			{
				unsigned j = y * bmh->width + x;
				if (j < bmh->imagesize)
				{
#pragma omp atomic
					bmi4[j]++;
				}
			}
		}
	}
}

void DoPlaceClose(int i, int j, unsigned short int ts, int tn, int DoAnyway)
{
	int k, ai, j1, j2, l, f, m, f2;
	unsigned short trig;
	unsigned short int t_start, t_stop;

#ifdef ABSENTEEISM_PLACE_CLOSURE
	unsigned short int t_old, t_new;
#endif


	f2 = 0;
	/*	if((j<0)||(j>=P.Nplace[i]))
			fprintf(stderr,"** %i %i *\n",i,j);
		else
	*/
#ifdef ABSENTEEISM_PLACE_CLOSURE
	t_new = ts / P.TimeStepsPerDay;
#endif
	trig = 0;
	t_start = ts + ((unsigned short int) (P.TimeStepsPerDay * P.PlaceCloseDelayMean));
	if (P.DoInterventionDelaysByAdUnit)
	{
		k = Mcells[Places[i][j].mcell].adunit;
		t_stop = ts + ((unsigned short int) (P.TimeStepsPerDay * (P.PlaceCloseDelayMean + AdUnits[k].PlaceCloseDuration)));
	}
	else
	{
		t_stop = ts + ((unsigned short int) (P.TimeStepsPerDay * (P.PlaceCloseDelayMean + P.PlaceCloseDuration)));
	}
#pragma omp critical (closeplace)
	{
		if (Places[i][j].close_end_time < t_stop)
		{
			if ((!DoAnyway) && (Places[i][j].control_trig < USHRT_MAX - 2))
			{
#ifdef ABSENTEEISM_PLACE_CLOSURE
				t_old = Places[i][j].AbsentLastUpdateTime;
				if (t_new >= t_old + MAX_ABSENT_TIME)
					for (l = 0; l < MAX_ABSENT_TIME; l++) Places[i][j].Absent[l] = 0;
				else
					for (l = t_old; l < t_new; l++) Places[i][j].Absent[l % MAX_ABSENT_TIME] = 0;
				for (l = t_new; l < t_new + P.usCaseAbsenteeismDuration / P.TimeStepsPerDay; l++) Places[i][j].Absent[l % MAX_ABSENT_TIME]++;
				trig = Places[i][j].Absent[t_new % MAX_ABSENT_TIME];
				Places[i][j].AbsentLastUpdateTime = t_new;
				if ((P.PlaceCloseByAdminUnit) && (P.PlaceCloseAdunitPlaceTypes[i] > 0)
					&& (((double)trig) / ((double)Places[i][j].n) > P.PlaceCloseCasePropThresh))
				{
					//fprintf(stderr,"** %i %i %i %i %lg ## ",i,j,(int) Places[i][j].control_trig, (int) Places[i][j].n,P.PlaceCloseCasePropThresh);
					k = Mcells[Places[i][j].mcell].adunit;
					if (AdUnits[k].place_close_trig < USHRT_MAX - 1) AdUnits[k].place_close_trig++;
				}
#else
				trig = Places[i][j].control_trig;
				if ((P.PlaceCloseByAdminUnit) && (P.PlaceCloseAdunitPlaceTypes[i] > 0)
					&& (((double)Places[i][j].control_trig) / ((double)Places[i][j].n) > P.PlaceCloseCasePropThresh))
				{
					//fprintf(stderr,"** %i %i %i %i %lg ## ",i,j,(int) Places[i][j].control_trig, (int) Places[i][j].n,P.PlaceCloseCasePropThresh);
					k = Mcells[Places[i][j].mcell].adunit;
					if (AdUnits[k].place_close_trig < USHRT_MAX - 1) AdUnits[k].place_close_trig++;
				}
#endif
			}
			if (Places[i][j].control_trig < USHRT_MAX - 1)
			{
				if (P.PlaceCloseFracIncTrig > 0)
					k = (((double)trig) / ((double)Places[i][j].n) > P.PlaceCloseFracIncTrig);
				else
					k = (((int)trig) >= P.PlaceCloseIncTrig);
				if (((!P.PlaceCloseByAdminUnit) && (k)) || (DoAnyway))
				{
					if (P.DoPlaceCloseOnceOnly)
						Places[i][j].control_trig = USHRT_MAX - 1;  // Places only close once
					else
						Places[i][j].control_trig = 0;
					if ((P.PlaceCloseEffect[i] == 0) || (ranf_mt(tn) >= P.PlaceCloseEffect[i]))
					{
						if (Places[i][j].close_start_time > t_start) Places[i][j].close_start_time = t_start;
						Places[i][j].close_end_time = t_stop;
						f2 = 1;
					}
					else
						Places[i][j].close_start_time = Places[i][j].close_end_time = t_stop;
				}
			}
		}
	}
	if (f2)
	{
		if (P.DoRealSymptWithdrawal)
			for (k = 0; k < Places[i][j].n; k++)
			{
				ai = Places[i][j].members[k];
				if ((HOST_AGE_YEAR(ai) < P.CaseAbsentChildAgeCutoff) && (!HOST_ABSENT(ai)) && (!HOST_QUARANTINED(ai)))
				{
					StateT[tn].cumAPCS++;
					if ((P.CaseAbsentChildPropAdultCarers == 1) || (ranf_mt(tn) < P.CaseAbsentChildPropAdultCarers))
					{
						j1 = Households[Hosts[ai].hh].FirstPerson; j2 = j1 + Households[Hosts[ai].hh].nh;
						if ((j1 < 0) || (j2 > P.N)) fprintf(stderr, "++ %i %i %i (%i %i %i)##  ", ai, j1, j2, i, j, k);
						f = 0;
						for (l = j1; (l < j2) && (!f); l++)
							f = ((abs(Hosts[l].inf) != InfStat_Dead) && (HOST_AGE_YEAR(l) >= P.CaseAbsentChildAgeCutoff) &&
							((Hosts[l].PlaceLinks[NUM_PLACE_TYPES_NOAIR - 1] < 0) || (HOST_ABSENT(l)) || (HOST_QUARANTINED(l))));
						if (!f)
						{
							for (l = j1; (l < j2) && (!f); l++)
								if ((HOST_AGE_YEAR(l) >= P.CaseAbsentChildAgeCutoff) && (abs(Hosts[l].inf) != InfStat_Dead)) { m = l; f = 1; }
							if (f)
							{
								if (Hosts[m].absent_start_time > Places[i][j].close_start_time) Hosts[m].absent_start_time = Places[i][j].close_start_time;
								if (Hosts[m].absent_stop_time < Places[i][j].close_end_time) Hosts[m].absent_stop_time = Places[i][j].close_end_time;
								StateT[tn].cumAPA++;
							}
						}
					}
				}
				if (Hosts[ai].absent_start_time > Places[i][j].close_start_time) Hosts[ai].absent_start_time = Places[i][j].close_start_time;
				if (Hosts[ai].absent_stop_time < Places[i][j].close_end_time) Hosts[ai].absent_stop_time = Places[i][j].close_end_time;
				if ((HOST_AGE_YEAR(ai) >= P.CaseAbsentChildAgeCutoff) && (Hosts[ai].PlaceLinks[NUM_PLACE_TYPES_NOAIR - 1] >= 0)) StateT[tn].cumAPC++;
			}
	}
}

void DoPlaceOpen(int i, int j, unsigned short int ts, int tn)
{
	int k, ai, j1, j2, l, f, m;

#ifdef ABSENTEEISM_PLACE_CLOSURE
	unsigned short int t_old, t_new;
	t_new = ts / P.TimeStepsPerDay;
#endif
#pragma omp critical (closeplace)
	{

		if (ts < Places[i][j].close_end_time)
		{

			if (P.DoRealSymptWithdrawal)
				for (k = 0; k < Places[i][j].n; k++)
				{
					ai = Places[i][j].members[k];
					if ((HOST_AGE_YEAR(ai) < P.CaseAbsentChildAgeCutoff) && (!HOST_ABSENT(ai)) && (!HOST_QUARANTINED(ai)))
					{
						StateT[tn].cumAPCS++;
						j1 = Households[Hosts[ai].hh].FirstPerson; j2 = j1 + Households[Hosts[ai].hh].nh;
						f = 0;
						for (l = j1; (l < j2) && (!f); l++)
							f = ((abs(Hosts[l].inf) != InfStat_Dead) && (HOST_AGE_YEAR(l) >= P.CaseAbsentChildAgeCutoff) &&
							((Hosts[l].PlaceLinks[NUM_PLACE_TYPES_NOAIR - 1] < 0) || (HOST_ABSENT(l)) || (HOST_QUARANTINED(l))));
						if (!f)
						{
							for (l = j1; (l < j2) && (!f); l++)
								if ((HOST_AGE_YEAR(l) >= P.CaseAbsentChildAgeCutoff) && (abs(Hosts[l].inf) != InfStat_Dead)) { m = l; f = 1; }
							if (f)
							{
								if (Hosts[m].absent_stop_time == Places[i][j].close_end_time) Hosts[m].absent_stop_time = ts;
							}
						}
					}
					if (Hosts[ai].absent_stop_time == Places[i][j].close_end_time) Hosts[ai].absent_stop_time = ts;
				}
			Places[i][j].close_end_time = ts;
		}
	}
}


int DoVacc(int ai, unsigned short int ts)
{
	int x, y;

	if ((!P.DoDistributionVaccination) && (State.cumV >= P.VaccMaxCourses))
		return 2;
	else if ((HOST_TO_BE_VACCED(ai)) || (Hosts[ai].inf < InfStat_InfectiousAlmostSymptomatic) || (Hosts[ai].inf >= InfStat_Dead_WasAsymp))
		return 1;
	else
	{
		Hosts[ai].vacc_start_time = ts + ((unsigned short int) (P.TimeStepsPerDay * P.VaccDelayMean));

#pragma omp critical (state_cumV)
		State.cumV++;
#pragma omp critical (state_cumV_daily)
		if (P.VaccDosePerDay >= 0)
		{
			State.cumV_daily++;
		}
#pragma omp critical (tot_vacc)
		Cells[Hosts[ai].pcell].tot_vacc++;
		if (P.OutputBitmap)
		{
			x = ((int)(Households[Hosts[ai].hh].loc_x * P.scalex)) - P.bminx;
			y = ((int)(Households[Hosts[ai].hh].loc_y * P.scaley)) - P.bminy;
			if ((x >= 0) && (x < P.bwidth) && (y >= 0) && (y < P.bheight))
			{
				unsigned j = y * bmh->width + x;
				if (j < bmh->imagesize)
				{
#pragma omp atomic
					bmi4[j]++;
				}
			}
		}
	}
	return 0;
}

void DoVaccNoDelay(int ai, unsigned short int ts)
{
	int x, y;

	if ((State.cumVG < P.VaccMaxCourses) && (!HOST_TO_BE_VACCED(ai)) && (Hosts[ai].inf >= InfStat_InfectiousAlmostSymptomatic) && (Hosts[ai].inf < InfStat_Dead_WasAsymp)) 
	{
		Hosts[ai].vacc_start_time = ts;
#pragma omp critical (state_cumVG) //changed to VG
		State.cumVG++; //changed to VG
#pragma omp critical (state_cumV_daily)
		if (P.VaccDosePerDay >= 0)
		{
			State.cumVG_daily++;
		}
#pragma omp critical (tot_vacc)
		Cells[Hosts[ai].pcell].tot_vacc++;
		if (P.OutputBitmap)
		{
			x = ((int)(Households[Hosts[ai].hh].loc_x * P.scalex)) - P.bminx;
			y = ((int)(Households[Hosts[ai].hh].loc_y * P.scaley)) - P.bminy;
			if ((x >= 0) && (x < P.bwidth) && (y >= 0) && (y < P.bheight))
			{
				unsigned j = y * bmh->width + x;
				if (j < bmh->imagesize)
				{
#pragma omp atomic
					bmi4[j]++;
				}
			}
		}
	}
}

///// Change person status functions (e.g. change person from susceptible to latently infected). 
int ChooseFinalDiseaseSeverity(int AgeGroup, int tn)
{
	int DiseaseSeverity;
	double x;

	// assume normalised props

	x = ranf_mt(tn);
	if (x < P.Prop_ILI_ByAge[AgeGroup]) DiseaseSeverity = Severity_ILI;
	else if (x < P.Prop_ILI_ByAge[AgeGroup] + P.Prop_SARI_ByAge[AgeGroup]) DiseaseSeverity = Severity_SARI;
	else if (x < P.Prop_ILI_ByAge[AgeGroup] + P.Prop_SARI_ByAge[AgeGroup] + P.Prop_Critical_ByAge[AgeGroup]) DiseaseSeverity = Severity_Critical;
	else DiseaseSeverity = Severity_Mild;
	return DiseaseSeverity;
}

unsigned short int ChooseFromICDF(double *ICDF, double Mean, int tn)
{
	unsigned short int Value;
	int i;
	double q, ti;

	i = (int)floor(q = ranf_mt(tn) * CDF_RES); //// note q defined here as well as i. 
	q -= ((double)i); //// remainder
	ti = -Mean * log(q * ICDF[i + 1] + (1.0 - q) * ICDF[i]); //// weighted average (sort of) between quartile values from CDF_RES. logged as it was previously exponentiated in ReadParams. Minus as exp(-cdf) was done in ReadParaams. Sort of
	Value = (unsigned short int) floor(0.5 + (ti * P.TimeStepsPerDay));

	return Value;
}
