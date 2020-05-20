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
Severity ChooseFinalDiseaseSeverity(int, int);

void DoImmune(int ai)
{
	// This transfers a person straight from susceptible to immune. Used to start a run with a partially immune population.
	Person* a;
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
					bmRecovered[j]++;
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
	double q, x, y; //// q radius squared, x and y coords. q later changed to be quantile of inverse CDF to choose latent period.
	Person* a;

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
				if ((ix >= 0) && (ix < P.bwidth) && (iy >= 0) && (iy < P.bheight))
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
			DoRecover(ai, tn, run);
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
		Person* a = Hosts + ai;
		if (a->Severity_Current == Severity_Asymptomatic)
		{
			a->Severity_Current = Severity_Mild;
			StateT[tn].Mild++;
			StateT[tn].cumMild++;
			StateT[tn].Mild_age[HOST_AGE_GROUP(ai)]++;
			StateT[tn].cumMild_age[HOST_AGE_GROUP(ai)]++;
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
		Person* a = Hosts + ai;
		if (a->Severity_Current == Severity_Asymptomatic)
		{
			a->Severity_Current = Severity_ILI;
			StateT[tn].ILI++;
			StateT[tn].cumILI++;
			StateT[tn].ILI_age[HOST_AGE_GROUP(ai)]++;
			StateT[tn].cumILI_age[HOST_AGE_GROUP(ai)]++;
			if (P.DoAdUnits)
			{
				StateT[tn].ILI_adunit	[Mcells[a->mcell].adunit]++;
				StateT[tn].cumILI_adunit[Mcells[a->mcell].adunit]++;
			}
		}
	}
}
void DoSARI(int ai, int tn)
{
	if (P.DoSeverity) //// shouldn't need this but best be careful.
	{
		Person* a = Hosts + ai;
		if (a->Severity_Current == Severity_ILI)
		{
			a->Severity_Current = Severity_SARI;
			StateT[tn].ILI--;
			StateT[tn].ILI_age[HOST_AGE_GROUP(ai)]--;
			StateT[tn].SARI++;
			StateT[tn].cumSARI++;
			StateT[tn].SARI_age[HOST_AGE_GROUP(ai)]++;
			StateT[tn].cumSARI_age[HOST_AGE_GROUP(ai)]++;
			if (P.DoAdUnits)
			{
				StateT[tn].ILI_adunit		[Mcells[a->mcell].adunit]--;
				StateT[tn].SARI_adunit		[Mcells[a->mcell].adunit]++;
				StateT[tn].cumSARI_adunit	[Mcells[a->mcell].adunit]++;
			}
		}
	}
}
void DoCritical(int ai, int tn)
{
	if (P.DoSeverity) //// shouldn't need this but best be careful.
	{
		Person* a = Hosts + ai;
		if (a->Severity_Current == Severity_SARI)
		{
			a->Severity_Current = Severity_Critical;
			StateT[tn].SARI--;
			StateT[tn].SARI_age[HOST_AGE_GROUP(ai)]--;
			StateT[tn].Critical++;
			StateT[tn].cumCritical++;
			StateT[tn].Critical_age[HOST_AGE_GROUP(ai)]++;
			StateT[tn].cumCritical_age[HOST_AGE_GROUP(ai)]++;
			if (P.DoAdUnits)
			{
				StateT[tn].SARI_adunit			[Mcells[a->mcell].adunit]--;
				StateT[tn].Critical_adunit		[Mcells[a->mcell].adunit]++;
				StateT[tn].cumCritical_adunit	[Mcells[a->mcell].adunit]++;
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
		Person* a = Hosts + ai;
		if (a->Severity_Current == Severity_Critical && (!a->to_die)) //// second condition should be unnecessary but leave in for now.
		{
			a->Severity_Current = Severity_RecoveringFromCritical;
			StateT[tn].Critical--;
			StateT[tn].Critical_age[HOST_AGE_GROUP(ai)]--;
			StateT[tn].CritRecov++;
			StateT[tn].cumCritRecov++;
			StateT[tn].CritRecov_age[HOST_AGE_GROUP(ai)]++;
			StateT[tn].cumCritRecov_age[HOST_AGE_GROUP(ai)]++;
			if (P.DoAdUnits)
			{
				StateT[tn].Critical_adunit[Mcells[a->mcell].adunit]--;
				StateT[tn].CritRecov_adunit[Mcells[a->mcell].adunit]++;
				StateT[tn].cumCritRecov_adunit[Mcells[a->mcell].adunit]++;
			}
		}
	}
}
void DoDeath_FromCriticalorSARIorILI(int ai, int tn)
{
	Person* a = Hosts + ai;
	if (P.DoSeverity)
	{
		if (a->Severity_Current == Severity_Critical)
		{
			StateT[tn].Critical--;
			StateT[tn].Critical_age[HOST_AGE_GROUP(ai)]--;
			StateT[tn].cumDeath_Critical++;
			StateT[tn].cumDeath_Critical_age[HOST_AGE_GROUP(ai)]++;
			if (P.DoAdUnits)
			{
				StateT[tn].Critical_adunit			[Mcells[a->mcell].adunit]--;
				StateT[tn].cumDeath_Critical_adunit	[Mcells[a->mcell].adunit]++;
			}
			//// change current status (so that flags work if function called again for same person). Don't move this outside of this if statement, even though it looks like it can be moved safely. It can't.
			a->Severity_Current = Severity_Dead;
		}
		else if (a->Severity_Current == Severity_SARI)
		{
			StateT[tn].SARI--;
			StateT[tn].SARI_age[HOST_AGE_GROUP(ai)]--;
			StateT[tn].cumDeath_SARI++;
			StateT[tn].cumDeath_SARI_age[HOST_AGE_GROUP(ai)]++;
			if (P.DoAdUnits)
			{
				StateT[tn].SARI_adunit			[Mcells[a->mcell].adunit]--;
				StateT[tn].cumDeath_SARI_adunit	[Mcells[a->mcell].adunit]++;
			}
			//// change current status (so that flags work if function called again for same person). Don't move this outside of this if statement, even though it looks like it can be moved safely. It can't.
			a->Severity_Current = Severity_Dead;
		}
		else if (a->Severity_Current == Severity_ILI)
		{
			StateT[tn].ILI--;
			StateT[tn].ILI_age[HOST_AGE_GROUP(ai)]--;
			StateT[tn].cumDeath_ILI++;
			StateT[tn].cumDeath_ILI_age[HOST_AGE_GROUP(ai)]++;
			if (P.DoAdUnits)
			{
				StateT[tn].ILI_adunit			[Mcells[a->mcell].adunit]--;
				StateT[tn].cumDeath_ILI_adunit	[Mcells[a->mcell].adunit]++;
			}
			//// change current status (so that flags work if function called again for same person). Don't move this outside of this if statement, even though it looks like it can be moved safely. It can't.
			a->Severity_Current = Severity_Dead;
		}
	}
}
void DoRecover_FromSeverity(int ai, int tn)
{
	//// note function different from DoRecoveringFromCritical.
	//// DoRecover_FromSeverity assigns people to state Recovered (and bookkeeps accordingly).
	//// DoRecoveringFromCritical assigns people to intermediate state "recovering from critical condition" (and bookkeeps accordingly).

	//// moved this from DoRecover
	Person* a = Hosts + ai;

	if (P.DoSeverity)
		if (a->inf == InfStat_InfectiousAsymptomaticNotCase || a->inf == InfStat_Case) ///// i.e same condition in DoRecover (make sure you don't recover people twice).
		{
			if (a->Severity_Current == Severity_Mild)
			{
				StateT[tn].Mild--;
				StateT[tn].Mild_age[HOST_AGE_GROUP(ai)]--;
				if (P.DoAdUnits) StateT[tn].Mild_adunit[Mcells[a->mcell].adunit]--;
				//// change current status (so that flags work if function called again for same person). Don't move this outside of this if statement, even though it looks like it can be moved safely. It can't.
				a->Severity_Current = Severity_Recovered;
			}
			else if (a->Severity_Current == Severity_ILI)
			{
				StateT[tn].ILI--;
				StateT[tn].ILI_age[HOST_AGE_GROUP(ai)]--;
				if (P.DoAdUnits) StateT[tn].ILI_adunit[Mcells[a->mcell].adunit]--;
				//// change current status (so that flags work if function called again for same person). Don't move this outside of this if statement, even though it looks like it can be moved safely. It can't.
				a->Severity_Current = Severity_Recovered;
			}
			else if (a->Severity_Current == Severity_SARI)
			{
				StateT[tn].SARI--;
				StateT[tn].SARI_age[HOST_AGE_GROUP(ai)]--;
				if (P.DoAdUnits) StateT[tn].SARI_adunit[Mcells[a->mcell].adunit]--;
				//// change current status (so that flags work if function called again for same person). Don't move this outside of this if statement, even though it looks like it can be moved safely. It can't.
				a->Severity_Current = Severity_Recovered;
			}
			else if (a->Severity_Current == Severity_RecoveringFromCritical)
			{
				StateT[tn].CritRecov--; //// decrement CritRecov, not critical.
				StateT[tn].CritRecov_age[HOST_AGE_GROUP(ai)]--;
				if (P.DoAdUnits) StateT[tn].CritRecov_adunit[Mcells[a->mcell].adunit]--;
				//// change current status (so that flags work if function called again for same person). Don't move this outside of this if statement, even though it looks like it can be moved safely. It can't.
				a->Severity_Current = Severity_Recovered;
			}
		}
}

void DoIncub(int ai, unsigned short int ts, int tn, int run)
{
	Person* a;
	double q;
	int age;

	age = HOST_AGE_GROUP(ai);
	if (age >= NUM_AGE_GROUPS) age = NUM_AGE_GROUPS - 1;

	a = Hosts + ai;
	if (a->inf == InfStat_Latent)
	{
		a->infectiousness = (float)P.AgeInfectiousness[age];
		if (P.InfectiousnessSD > 0) a->infectiousness *= (float) gen_gamma_mt(1 / (P.InfectiousnessSD * P.InfectiousnessSD), 1 / (P.InfectiousnessSD * P.InfectiousnessSD), tn);
		q = P.ProportionSymptomatic[age]
			* (HOST_TREATED(ai) ? (1 - P.TreatSympDrop) : 1)
			* (HOST_VACCED(ai) ? (1 - P.VaccSympDrop) : 1);

		if (ranf_mt(tn) < q)
		{
			a->inf = InfStat_InfectiousAlmostSymptomatic;
			a->infectiousness *= (float)(-P.SymptInfectiousness);
		}
		else
			a->inf = InfStat_InfectiousAsymptomaticNotCase;

		if (!P.DoSeverity || a->inf == InfStat_InfectiousAsymptomaticNotCase) //// if not doing severity or if person asymptomatic.
		{
			if (P.DoInfectiousnessProfile)	a->recovery_or_death_time = a->latent_time + (unsigned short int) (P.InfectiousPeriod * P.TimeStepsPerDay);
			else							a->recovery_or_death_time = a->latent_time + ChooseFromICDF(P.infectious_icdf, P.InfectiousPeriod, tn);
		}
		else
		{
			int CaseTime = a->latent_time + ((int)(P.LatentToSymptDelay / P.TimeStep)); //// base severity times on CaseTime, not latent time. Otherwise there are edge cases where recovery time is zero days after latent_time and therefore before DoCase called in IncubRecoverySweep (i.e. people can recover before they've become a case!).

			//// choose final disease severity (either mild, ILI, SARI, Critical, not asymptomatic as covered above) by age
			a->Severity_Final = ChooseFinalDiseaseSeverity(age, tn);

			/// choose outcome recovery or death
			if (	((a->Severity_Final == Severity_Critical)	&& (ranf_mt(tn) < P.CFR_Critical_ByAge	[age]))		||
					((a->Severity_Final == Severity_SARI	)	&& (ranf_mt(tn) < P.CFR_SARI_ByAge		[age]))		||
					((a->Severity_Final == Severity_ILI		)	&& (ranf_mt(tn) < P.CFR_ILI_ByAge		[age]))		)
				a->to_die = 1;

			//// choose events and event times
			if (a->Severity_Final == Severity_Mild)
				a->recovery_or_death_time = CaseTime + ChooseFromICDF(P.MildToRecovery_icdf, P.Mean_MildToRecovery[age], tn);
			else if (a->Severity_Final == Severity_Critical)
			{
				a->SARI_time		= CaseTime		+ ChooseFromICDF(P.ILIToSARI_icdf		, P.Mean_ILIToSARI[age], tn);
				a->Critical_time	= a->SARI_time	+ ChooseFromICDF(P.SARIToCritical_icdf	, P.Mean_SARIToCritical[age], tn);
				if (a->to_die)
					a->recovery_or_death_time = a->Critical_time					+ ChooseFromICDF(P.CriticalToDeath_icdf		, P.Mean_CriticalToDeath[age], tn);
				else
				{
					a->RecoveringFromCritical_time	= a->Critical_time					+ ChooseFromICDF(P.CriticalToCritRecov_icdf	, P.Mean_CriticalToCritRecov[age], tn);
					a->recovery_or_death_time		= a->RecoveringFromCritical_time	+ ChooseFromICDF(P.CritRecovToRecov_icdf	, P.Mean_CritRecovToRecov[age], tn);
				}
			}
			else if (a->Severity_Final == Severity_SARI)
			{
				a->SARI_time = CaseTime + ChooseFromICDF(P.ILIToSARI_icdf, P.Mean_ILIToSARI[age], tn);
				if (a->to_die)
					a->recovery_or_death_time = a->SARI_time + ChooseFromICDF(P.SARIToDeath_icdf	, P.Mean_SARIToDeath[age], tn);
				else
					a->recovery_or_death_time = a->SARI_time + ChooseFromICDF(P.SARIToRecovery_icdf	, P.Mean_SARIToRecovery[age], tn);
			}
			else /*i.e. if Severity_Final == Severity_ILI*/
			{
				if (a->to_die)
					a->recovery_or_death_time = CaseTime + ChooseFromICDF(P.ILIToDeath_icdf		, P.Mean_ILIToDeath[age], tn);
				else
					a->recovery_or_death_time = CaseTime + ChooseFromICDF(P.ILIToRecovery_icdf	, P.Mean_ILIToRecovery[age], tn);
			}
		}

		if ((a->inf== InfStat_InfectiousAlmostSymptomatic) && ((P.ControlPropCasesId == 1) || (ranf_mt(tn) < P.ControlPropCasesId)))
		{
			Hosts[ai].detected = 1;
			Hosts[ai].detected_time = ts + (unsigned short int)(P.LatentToSymptDelay * P.TimeStepsPerDay);


			if ((P.DoDigitalContactTracing) && (Hosts[ai].detected_time >= (unsigned short int)(AdUnits[Mcells[Hosts[ai].mcell].adunit].DigitalContactTracingTimeStart * P.TimeStepsPerDay)) && (Hosts[ai].detected_time < (unsigned short int)((AdUnits[Mcells[Hosts[ai].mcell].adunit].DigitalContactTracingTimeStart + P.DigitalContactTracingPolicyDuration)*P.TimeStepsPerDay)) && (Hosts[ai].digitalContactTracingUser))
			{
				//set dct_trigger_time for index case
			if (P.DoDigitalContactTracing)	//set dct_trigger_time for index case
				if (Hosts[ai].dct_trigger_time == (USHRT_MAX - 1)) //if this hasn't been set in DigitalContactTracingSweep due to detection of contact of contacts, set it here
					Hosts[ai].dct_trigger_time = Hosts[ai].detected_time + (unsigned short int) (P.DelayFromIndexCaseDetectionToDCTIsolation * P.TimeStepsPerDay);
			}
		}

		//// update pointers
		Cells[a->pcell].L--;		//// one fewer person latently infected.
		Cells[a->pcell].infected--; //// first infected person is now one index earlier in array.
		Cells[a->pcell].I++;		//// one more infectious person.
		if (Cells[a->pcell].L > 0)
		{
			Cells[a->pcell].susceptible[a->listpos] = Cells[a->pcell].latent[Cells[a->pcell].L]; //// reset pointers.
			Hosts[Cells[a->pcell].susceptible[a->listpos]].listpos = a->listpos;
			a->listpos = Cells[a->pcell].S + Cells[a->pcell].L; //// change person a's listpos, which will now refer to their position among infectious people, not latent.
			Cells[a->pcell].infected[0] = ai; //// this person is now first infectious person in the array. Pointer was moved back one so now that memory address refers to person ai. Alternative would be to move everyone back one which would take longer.
		}
	}
}

void DoDetectedCase(int ai, double t, unsigned short int ts, int tn)
{
	//// Function DoDetectedCase does many things associated with various interventions.
	//// Enacts Household quarantine, case isolation, place closure.
	//// and therefore changes lots of quantities (e.g. quar_comply and isolation_start_time) associated with model macros e.g. HOST_ABSENT / HOST_ISOLATED

	int j, k, f, j1, j2, ad; // m, h, ad;
	Person* a = Hosts + ai;

	//// Increment triggers (Based on numbers of detected cases) for interventions. Used in TreatSweep function when not doing Global or Admin triggers. And not when doing ICU triggers.
	if (Mcells[a->mcell].treat_trig				< USHRT_MAX - 1) Mcells[a->mcell].treat_trig++;
	if (Mcells[a->mcell].vacc_trig				< USHRT_MAX - 1) Mcells[a->mcell].vacc_trig++;
	if (Mcells[a->mcell].move_trig				< USHRT_MAX - 1) Mcells[a->mcell].move_trig++;
	if (Mcells[a->mcell].socdist_trig			< USHRT_MAX - 1) Mcells[a->mcell].socdist_trig++;
	if (Mcells[a->mcell].keyworkerproph_trig	< USHRT_MAX - 1) Mcells[a->mcell].keyworkerproph_trig++;

	if (!P.AbsenteeismPlaceClosure)
	{
		if ((P.PlaceCloseRoundHousehold)&& (Mcells[a->mcell].place_trig < USHRT_MAX - 1)) Mcells[a->mcell].place_trig++;
		if ((t >= P.PlaceCloseTimeStart) && (!P.DoAdminTriggers) && (!((P.DoGlobalTriggers)&&(P.PlaceCloseCellIncThresh<1000000000))))
			for (j = 0; j < P.PlaceTypeNum; j++)
				if ((j != P.HotelPlaceType) && (a->PlaceLinks[j] >= 0))
				{
					DoPlaceClose(j, a->PlaceLinks[j], ts, tn, 0);
					if (!P.PlaceCloseRoundHousehold)
					{
						if (Mcells[Places[j][a->PlaceLinks[j]].mcell].place_trig < USHRT_MAX - 1)
						{
#pragma omp critical (place_trig)
							Mcells[Places[j][a->PlaceLinks[j]].mcell].place_trig++;
						}
					}
				}
	}

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
				Hosts[j1].quar_start_time = ts + ((unsigned short int) (P.TimeStepsPerDay * P.HQuarantineDelay));
				k = (ranf_mt(tn) < P.HQuarantinePropHouseCompliant) ? 1 : 0; //// Is household compliant? True or false
				if (k) StateT[tn].cumHQ++; ////  if compliant, increment cumulative numbers of households under quarantine.
				//// if household not compliant then neither is first person. Otheswise ask whether first person is compliant?
				///// cycle through remaining household members and repeat the above steps
				for (j = j1; j < j2; j++)
				{
					if(j>j1) Hosts[j].quar_start_time = Hosts[j1].quar_start_time;
					Hosts[j].quar_comply = ((k == 0) ? 0 : ((ranf_mt(tn) < P.HQuarantinePropIndivCompliant) ? 1 : 0));
					if ((Hosts[j].quar_comply) && (!HOST_ABSENT(j)))
					{
						if (HOST_AGE_YEAR(j) >= P.CaseAbsentChildAgeCutoff)
						{
							if (Hosts[j].PlaceLinks[P.PlaceTypeNoAirNum - 1] >= 0) StateT[tn].cumAH++;
						}
						else	StateT[tn].cumACS++;
					}
				}
			}
		}
	}

	//// Giant compound if statement. If doing delays by admin unit, then window of case isolation dependent on admin unit-specific duration. This if statement ensures that this timepoint within window, regardless of how window defined.
	if ((P.DoInterventionDelaysByAdUnit &&
		(t >= AdUnits[Mcells[a->mcell].adunit].CaseIsolationTimeStart && (t < AdUnits[Mcells[a->mcell].adunit].CaseIsolationTimeStart + AdUnits[Mcells[a->mcell].adunit].CaseIsolationPolicyDuration)))	||
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
					if ((!HOST_QUARANTINED(ai)) && (Hosts[ai].PlaceLinks[P.PlaceTypeNoAirNum - 1] >= 0) && (HOST_AGE_YEAR(ai) >= P.CaseAbsentChildAgeCutoff))
						StateT[tn].cumAC++;
				}
				if ((P.DoHouseholds) && (P.DoPlaces) && (HOST_AGE_YEAR(ai) < P.CaseAbsentChildAgeCutoff)) //// if host is a child who requires adult to stay at home.
				{
					if (!HOST_QUARANTINED(ai)) StateT[tn].cumACS++;
					if (Hosts[ai].ProbCare < P.CaseAbsentChildPropAdultCarers) //// if adult needs to stay at home (i.e. if Proportion of children at home for whom one adult also stays at home = 1 or coinflip satisfied.)
					{
						j1 = Households[Hosts[ai].hh].FirstPerson; j2 = j1 + Households[Hosts[ai].hh].nh;
						f = 0;

						//// in loop below, f true if any household member a) alive AND b) not a child AND c) has no links to workplace (or is absent from work or quarantined).
						for (j = j1; (j < j2) && (!f); j++)
							f = ((abs(Hosts[j].inf) != InfStat_Dead) && (HOST_AGE_YEAR(j) >= P.CaseAbsentChildAgeCutoff) && ((Hosts[j].PlaceLinks[P.PlaceTypeNoAirNum - 1] < 0) || (HOST_ABSENT(j)) || (HOST_QUARANTINED(j))));

						//// so !f true if any household member EITHER: a) dead; b) a child; c) has a link to an office and not currently absent or quarantined.
						if (!f) //// so if either a) a household member is dead; b) a household member is a child requiring adult to stay home; c) a household member has links to office.
						{
							for (j = j1; (j < j2) & (!f); j++) /// loop again, checking whether household members not children needing supervision and are alive.
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

	//add contacts to digital contact tracing, but only if considering contact tracing, we are within the window of the policy and the detected case is a user
	if ((P.DoDigitalContactTracing) && (t >= AdUnits[Mcells[Hosts[ai].mcell].adunit].DigitalContactTracingTimeStart) && (t < AdUnits[Mcells[Hosts[ai].mcell].adunit].DigitalContactTracingTimeStart + P.DigitalContactTracingPolicyDuration) && (Hosts[ai].digitalContactTracingUser))
	{

		// allow for DCT to isolate index cases
		if ((P.DCTIsolateIndexCases) && (Hosts[ai].index_case_dct==0))//(Hosts[ai].digitalContactTraced == 0)&& - currently removed this condition as it would mean that someone already under isolation wouldn't have their isolation extended
		{
			ad = Mcells[Hosts[ai].mcell].adunit;
			//if (AdUnits[j].ndct < AdUnits[j].n)
			if(StateT[tn].ndct_queue[ad] < AdUnits[ad].n)
			{
				//if we are isolating an index case, we set their infector as -1 in order to get the timings consistent.
				StateT[tn].dct_queue[ad][StateT[tn].ndct_queue[ad]++] = { ai,-1,ts };
			}
			else
			{
				fprintf(stderr, "No more space in queue! AdUnit: %i, ndct=%i, max queue length: %i\n", ad, AdUnits[j].ndct, AdUnits[ad].n);
				fprintf(stderr, "Error!\n");
			}
		}
		//currently commenting this out as household members will likely be picked up by household quarantine.
		//can add back in if needed, but would need to re-add a couple more local variables.

		//if(P.IncludeHouseholdDigitalContactTracing)
		//{
		//	//Then we want to find all their household and place group contacts to add to the contact tracing queue
		//	//Start with household contacts
		//	j1 = Households[Hosts[ai].hh].FirstPerson; j2 = j1 + Households[Hosts[ai].hh].nh;
		//	for (j = j1; j < j2; j++)
		//	{
		//		//if host is dead or the detected case, no need to add them to the list. They also need to be a user themselves
		//		if ((abs(Hosts[j].inf) != 5) && (j != ai) && (Hosts[j].digitalContactTracingUser) && (ranf_mt(tn)<P.ProportionDigitalContactsIsolate))
		//		{
		//			//add contact and detected infectious host to lists
		//			ad = Mcells[Hosts[j].mcell].adunit;
		//			if ((StateT[tn].ndct_queue[ad] < P.InfQueuePeakLength))
		//			{
		//				StateT[tn].dct_queue[ad][StateT[tn].ndct_queue[ad]++] = j;
		//				StateT[tn].contacts[ad][StateT[tn].ncontacts[ad]++] = ai;
		//			}
		//			else
		//			{
		//				fprintf(stderr, "No space left in queue! Thread: %i, AdUnit: %i\n", tn, ad);
		//			}
		//		}
		//	}
		//}
		//if(P.IncludePlaceGroupDigitalContactTracing)
		//{
		//	//then loop over place group contacts as well
		//	for (int i = 0; i < P.PlaceTypeNum; i++)
		//	{
		//		k = Hosts[ai].PlaceLinks[i];
		//		if (k >= 0)
		//		{
		//			//Find place group link
		//			m = Hosts[ai].PlaceGroupLinks[i];
		//			j1 = Places[i][k].group_start[m]; j2 = j1 + Places[i][k].group_size[m];
		//			for (j = j1; j < j2; j++)
		//			{
		//				h = Places[i][k].members[j];
		//				ad = Mcells[Hosts[h].mcell].adunit;
		//				//if host is dead or the detected case, no need to add them to the list. They also need to be a user themselves
		//				if ((abs(Hosts[h].inf) != 5) && (h != ai) && (Hosts[h].digitalContactTracingUser))// && (ranf_mt(tn)<P.ProportionDigitalContactsIsolate))
		//				{
		//					ad = Mcells[Hosts[h].mcell].adunit;
		//					if ((StateT[tn].ndct_queue[ad] < P.InfQueuePeakLength))
		//					{
		//						//PLEASE CHECK ALL THIS LOGIC CAREFULLY!
		//
		//						StateT[tn].dct_queue[ad][StateT[tn].ndct_queue[ad]++] = h;
		//						StateT[tn].contacts[ad][StateT[tn].ncontacts[ad]++] = ai; //keep a record of who the detected case was
		//					}
		//					else
		//					{
		//						fprintf(stderr, "No space left in queue! Thread: %i, AdUnit: %i\n", tn, ad);
		//					}
		//
		//				}
		//			}
		//		}
		//	}
		//}

	}

}

void DoCase(int ai, double t, unsigned short int ts, int tn) //// makes an infectious (but asymptomatic) person symptomatic. Called in IncubRecoverySweep (and DoInfect if P.DoOneGen)
{
	int j, k, f, j1, j2;
	Person* a;
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
		else if((P.DoRealSymptWithdrawal)&&(P.DoPlaces))
		{
			a->absent_start_time = USHRT_MAX - 1;
			for (j = 0; j < P.PlaceTypeNum; j++)
				if ((a->PlaceLinks[j] >= 0) && (j != P.HotelPlaceType) && (!HOST_ABSENT(ai)) && (P.SymptPlaceTypeWithdrawalProp[j] > 0))
				{
					if ((P.SymptPlaceTypeWithdrawalProp[j] == 1) || (ranf_mt(tn) < P.SymptPlaceTypeWithdrawalProp[j]))
					{
						a->absent_start_time = ts + P.usCaseAbsenteeismDelay;
						a->absent_stop_time = ts + P.usCaseAbsenteeismDelay + P.usCaseAbsenteeismDuration;
						if (P.AbsenteeismPlaceClosure)
						{
							if ((t >= P.PlaceCloseTimeStart) && (!P.DoAdminTriggers) && (!P.DoGlobalTriggers))
								for (j = 0; j < P.PlaceTypeNum; j++)
									if ((j != P.HotelPlaceType) && (a->PlaceLinks[j] >= 0))
											DoPlaceClose(j, a->PlaceLinks[j], ts, tn, 0);
						}
						if ((!HOST_QUARANTINED(ai)) && (Hosts[ai].PlaceLinks[P.PlaceTypeNoAirNum - 1] >= 0) && (HOST_AGE_YEAR(ai) >= P.CaseAbsentChildAgeCutoff))
							StateT[tn].cumAC++;
						/* This calculates adult absenteeism from work due to care of sick children. Note, children not at school not counted (really this should
						be fixed in population setup by having adult at home all the time for such kids. */
						if ((P.DoHouseholds) && (HOST_AGE_YEAR(ai) < P.CaseAbsentChildAgeCutoff))
						{
							if (!HOST_QUARANTINED(ai)) StateT[tn].cumACS++;
							if (Hosts[ai].ProbCare < P.CaseAbsentChildPropAdultCarers)
							{
								j1 = Households[Hosts[ai].hh].FirstPerson; j2 = j1 + Households[Hosts[ai].hh].nh;
								f = 0;
								for (int j = j1; (j < j2) && (!f); j++)
									f = ((abs(Hosts[j].inf) != InfStat_Dead) && (HOST_AGE_YEAR(j) >= P.CaseAbsentChildAgeCutoff) && ((Hosts[j].PlaceLinks[P.PlaceTypeNoAirNum - 1] < 0)|| (HOST_ABSENT(j)) || (HOST_QUARANTINED(j))));
								if (!f)
								{
									for (int j = j1; (j < j2) && (!f); j++)
										if ((HOST_AGE_YEAR(j) >= P.CaseAbsentChildAgeCutoff) && (abs(Hosts[j].inf) != InfStat_Dead)) { k = j; f = 1; }
									if (f)
									{
										if (!HOST_ABSENT(k)) Hosts[k].absent_start_time = ts + P.usCaseIsolationDelay;
										Hosts[k].absent_stop_time = ts + P.usCaseIsolationDelay + P.usCaseIsolationDuration;
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
			//if ((P.ControlPropCasesId == 1) || (ranf_mt(tn) < P.ControlPropCasesId))
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


		if (P.DoSeverity)
		{
			if (a->Severity_Final == Severity_Mild)
				DoMild(ai, tn);
			else
				DoILI(ai, tn); //// symptomatic cases either mild or ILI at symptom onset. SARI and Critical cases still onset with ILI.
		}
		if (P.DoAdUnits) StateT[tn].cumC_adunit[Mcells[a->mcell].adunit]++;
	}
}

void DoFalseCase(int ai, double t, unsigned short int ts, int tn)
{
	/* Arguably adult absenteeism to take care of sick kids could be included here, but then output absenteeism would not be 'excess' absenteeism */
	if ((P.ControlPropCasesId == 1) || (ranf_mt(tn) < P.ControlPropCasesId))
	{
		if ((!P.DoEarlyCaseDiagnosis) || (State.cumDC >= P.PreControlClusterIdCaseThreshold)) StateT[tn].cumDC++;
		DoDetectedCase(ai, t, ts, tn);
	}
	StateT[tn].cumFC++;
}

void DoRecover(int ai, int tn, int run)
{
	int i, j, x, y;
	Person* a;

	a = Hosts + ai;
	if (a->inf == InfStat_InfectiousAsymptomaticNotCase || a->inf == InfStat_Case)
	{
		i = a->listpos;
		Cells[a->pcell].I--; //// one less infectious person
		Cells[a->pcell].R++; //// one more recovered person
		j = Cells[a->pcell].S + Cells[a->pcell].L + Cells[a->pcell].I;
		if (i < Cells[a->pcell].S + Cells[a->pcell].L + Cells[a->pcell].I)
		{
			Cells[a->pcell].susceptible[i] = Cells[a->pcell].susceptible[j];
			Hosts[Cells[a->pcell].susceptible[i]].listpos = i;
			a->listpos = j;
			Cells[a->pcell].susceptible[j] = ai;
		}
		a->inf = (InfStat)(InfStat_Recovered * a->inf / abs(a->inf));

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
						bmRecovered[j]++;
#pragma omp atomic
						bmInfected[j]--;
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
						bmRecovered[j]++;
#pragma omp atomic
						bmInfected[j]--;
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
						bmTreated[j]++;
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
					bmTreated[j]++;
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
					bmTreated[j]++;
				}
			}
		}
	}
}

void DoPlaceClose(int i, int j, unsigned short int ts, int tn, int DoAnyway)
{
	//// DoPlaceClose function called in TreatSweep (with arg DoAnyway = 1) and DoDetectedCase (with arg DoAnyway = 0).
	//// Basic pupose of this function is to change Places[i][j].close_start_time and Places[i][j].close_end_time, so that macro PLACE_CLOSED will return true.
	//// This will then scale peoples household, place, and spatial infectiousness and susceptibilities in function InfectSweep (but not in functions ini CalcInfSusc.cpp)

	int k, ai, j1, j2, l, f, m, f2;
	unsigned short trig;
	unsigned short int t_start, t_stop;
	unsigned short int t_old, t_new;

	f2 = 0;
	/*	if((j<0)||(j>=P.Nplace[i]))
			fprintf(stderr,"** %i %i *\n",i,j);
		else
	*/
	t_new = (unsigned short)(((double)ts) / P.TimeStepsPerDay);
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
		//// close_start_time initialized to USHRT_MAX - 1.
		//// close_end_time initialized to zero in InitModel (so will pass this check on at least first call of this function).

		if (Places[i][j].close_end_time < t_stop)
		{
			if ((!DoAnyway) && (Places[i][j].control_trig < USHRT_MAX - 2))
			{
				Places[i][j].control_trig++;
				if (P.AbsenteeismPlaceClosure)
				{
					t_old = Places[i][j].AbsentLastUpdateTime;
					if (t_new >= t_old + P.MaxAbsentTime)
						for (l = 0; l < P.MaxAbsentTime; l++) Places[i][j].Absent[l] = 0;
					else
						for (l = t_old; l < t_new; l++) Places[i][j].Absent[l % P.MaxAbsentTime] = 0;
					for (l = t_new; l < t_new + P.usCaseAbsenteeismDuration / P.TimeStepsPerDay; l++) Places[i][j].Absent[l % P.MaxAbsentTime]++;
					trig = Places[i][j].Absent[t_new % P.MaxAbsentTime];
					Places[i][j].AbsentLastUpdateTime = t_new;
					if ((P.PlaceCloseByAdminUnit) && (P.PlaceCloseAdunitPlaceTypes[i] > 0)
						&& (((double)trig) / ((double)Places[i][j].n) > P.PlaceCloseCasePropThresh))
					{
						//fprintf(stderr,"** %i %i %i %i %lg ## ",i,j,(int) Places[i][j].control_trig, (int) Places[i][j].n,P.PlaceCloseCasePropThresh);
						k = Mcells[Places[i][j].mcell].adunit;
						if (AdUnits[k].place_close_trig < USHRT_MAX - 1) AdUnits[k].place_close_trig++;
					}
				}
				else
				{
					trig = Places[i][j].control_trig;
					if ((P.PlaceCloseByAdminUnit) && (P.PlaceCloseAdunitPlaceTypes[i] > 0)
						&& (((double)Places[i][j].control_trig) / ((double)Places[i][j].n) > P.PlaceCloseCasePropThresh))
					{
						//fprintf(stderr,"** %i %i %i %i %lg ## ",i,j,(int) Places[i][j].control_trig, (int) Places[i][j].n,P.PlaceCloseCasePropThresh);
						k = Mcells[Places[i][j].mcell].adunit;
						if (AdUnits[k].place_close_trig < USHRT_MAX - 1) AdUnits[k].place_close_trig++;
					}
				}
			}
			if (Places[i][j].control_trig < USHRT_MAX - 1) //// control_trig initialized to zero so this check will pass at least once
			{
				if (P.PlaceCloseFracIncTrig > 0)
					k = (((double)trig) / ((double)Places[i][j].n) > P.PlaceCloseFracIncTrig);
				else
					k = (((int)trig) >= P.PlaceCloseIncTrig);
				if (((!P.PlaceCloseByAdminUnit) && (k)) || (DoAnyway))
				{
					if (P.DoPlaceCloseOnceOnly)
						Places[i][j].control_trig = USHRT_MAX - 1;  //// Places only close once, and so this code block would not be entered again.
					else
						Places[i][j].control_trig = 0;				//// otherwise reset the trigger.

					//// set close_start_time and close_end_time

					if (Places[i][j].ProbClose >= P.PlaceCloseEffect[i]) //// if proportion of places of type i remaining open is 0 or if place is closed with prob 1 - PlaceCloseEffect[i]...
					{
						if (Places[i][j].close_start_time > t_start) Places[i][j].close_start_time = t_start;
						Places[i][j].close_end_time = t_stop;
						f2 = 1; /// /set flag to true so next part of function used.
					}
					else
						Places[i][j].close_start_time = Places[i][j].close_end_time = t_stop; //// ... otherwise set start and end of closure to be the same, which will cause macro PLACE_CLOSED to always return false.
				}
			}
		}
	}

	if (f2)
	{
		if (P.DoRealSymptWithdrawal)
			for (k = 0; k < Places[i][j].n; k++) //// loop over all people in place.
			{
				ai = Places[i][j].members[k];
				if (((P.PlaceClosePropAttending[i] == 0) || (Hosts[ai].ProbAbsent >= P.PlaceClosePropAttending[i])))
				{
					if ((!HOST_ABSENT(ai)) && (!HOST_QUARANTINED(ai)) && (HOST_AGE_YEAR(ai) < P.CaseAbsentChildAgeCutoff)) //// if person is a child and neither absent nor quarantined
					{
						StateT[tn].cumAPCS++;
						if (Hosts[ai].ProbCare < P.CaseAbsentChildPropAdultCarers) //// if child needs adult supervision
						{
							j1 = Households[Hosts[ai].hh].FirstPerson; j2 = j1 + Households[Hosts[ai].hh].nh;
							if ((j1 < 0) || (j2 > P.PopSize)) fprintf(stderr, "++ %i %i %i (%i %i %i)##  ", ai, j1, j2, i, j, k);
							f = 0;

							//// in loop below, f true if any household member a) alive AND b) not a child AND c) has no links to workplace (or is absent from work or quarantined).
							for (l = j1; (l < j2) && (!f); l++)
								f = ((abs(Hosts[l].inf) != InfStat_Dead) && (HOST_AGE_YEAR(l) >= P.CaseAbsentChildAgeCutoff) && ((Hosts[l].PlaceLinks[P.PlaceTypeNoAirNum - 1] < 0) || (HOST_QUARANTINED(l))));
							if (!f) //// so !f true if there's no living adult household member who is not quarantined already or isn't a home-worker.
							{
								for (l = j1; (l < j2) && (!f); l++) //// loop over all household members of child this place: find the adults and ensure they're not dead...
									if ((HOST_AGE_YEAR(l) >= P.CaseAbsentChildAgeCutoff) && (abs(Hosts[l].inf) != InfStat_Dead))
									{
										if (Hosts[l].absent_start_time > t_start) Hosts[l].absent_start_time = t_start;
										if (Hosts[l].absent_stop_time < t_stop) Hosts[l].absent_stop_time = t_stop;
										StateT[tn].cumAPA++;
										f = 1;
									}
							}
						}
					}
					//#pragma omp critical (closeplace3)
					{
						///// finally amend absent start and stop times if they contradict place start and stop times.
						if (Hosts[ai].absent_start_time > t_start) Hosts[ai].absent_start_time = t_start;
						if (Hosts[ai].absent_stop_time < t_stop) Hosts[ai].absent_stop_time = t_stop;
					}
					if ((HOST_AGE_YEAR(ai) >= P.CaseAbsentChildAgeCutoff) && (Hosts[ai].PlaceLinks[P.PlaceTypeNoAirNum - 1] >= 0)) StateT[tn].cumAPC++;
				}
			}
	}
}


void DoPlaceOpen(int i, int j, unsigned short int ts, int tn)
{
	int k, ai, j1, j2, l, f;

#pragma omp critical (openplace)
	{
		if (ts < Places[i][j].close_end_time)
		{
			if (P.DoRealSymptWithdrawal)
				for (k = 0; k < Places[i][j].n; k++)
				{
					ai = Places[i][j].members[k];
					if (Hosts[ai].absent_stop_time == Places[i][j].close_end_time) Hosts[ai].absent_stop_time = ts;
					if (Hosts[ai].ProbCare < P.CaseAbsentChildPropAdultCarers) //// if child needs adult supervision
					{
						if ((HOST_AGE_YEAR(ai) < P.CaseAbsentChildAgeCutoff) && (!HOST_QUARANTINED(ai)))
						{
							j1 = Households[Hosts[ai].hh].FirstPerson; j2 = j1 + Households[Hosts[ai].hh].nh;
							f = 0;
							for (l = j1; (l < j2) && (!f); l++)
								f = ((abs(Hosts[l].inf) != InfStat_Dead) && (HOST_AGE_YEAR(l) >= P.CaseAbsentChildAgeCutoff) && ((Hosts[l].PlaceLinks[P.PlaceTypeNoAirNum - 1] < 0) || (HOST_QUARANTINED(l))));
							if (!f)
							{
								for (l = j1; (l < j2) && (!f); l++)
									if ((HOST_AGE_YEAR(l) >= P.CaseAbsentChildAgeCutoff) && (abs(Hosts[l].inf) != InfStat_Dead) && (HOST_ABSENT(l)))
									{
										if (Hosts[l].absent_stop_time == Places[i][j].close_end_time) Hosts[l].absent_stop_time = ts;
									}
							}
						}
					}
				}
			Places[i][j].close_end_time = ts;
		}
	}
}

int DoVacc(int ai, unsigned short int ts)
{
	int x, y;

	if (State.cumV >= P.VaccMaxCourses)
		return 2;
	else if ((HOST_TO_BE_VACCED(ai)) || (Hosts[ai].inf < InfStat_InfectiousAlmostSymptomatic) || (Hosts[ai].inf >= InfStat_Dead_WasAsymp))
		return 1;
	else
	{
		Hosts[ai].vacc_start_time = ts + ((unsigned short int) (P.TimeStepsPerDay * P.VaccDelayMean));

#pragma omp critical (state_cumV)
		State.cumV++;
		if (P.VaccDosePerDay >= 0)
		{
#pragma omp critical (state_cumV_daily)
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
					bmTreated[j]++;
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
		if (P.VaccDosePerDay >= 0)
		{
#pragma omp critical (state_cumV_daily)
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
					bmTreated[j]++;
				}
			}
		}
	}
}

///// Change person status functions (e.g. change person from susceptible to latently infected).
Severity ChooseFinalDiseaseSeverity(int AgeGroup, int tn)
{
	Severity DiseaseSeverity;
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
