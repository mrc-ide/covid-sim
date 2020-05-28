#include "Latent.h"
#include "Model.h"
#include "ModelMacros.h"
#include "Rand.h"
#include "InfectiousAlmostSymptomatic.h"
#include "InfectiousAsymptomaticNotCase.h"
#include "Update.h"

void Latent::GetsWorse(int ai, double t, int tn, int run)
{
	
	int inf_stat = DoIncub(ai, t, tn, run); // TODO

	if (inf_stat == InfStat_InfectiousAlmostSymptomatic)
		Hosts->infectionState = new InfectiousAlmostSymptomatic();

	if (inf_stat == InfStat_InfectiousAsymptomaticNotCase)
		Hosts->infectionState = new InfectiousAsymptomaticNotCase();

}

void Latent::GetsBetter(int ai, double t, int tn, int run)
{

}

int Latent::DoIncub(int ai, unsigned short int ts, int tn, int run)
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
		if (P.InfectiousnessSD > 0) a->infectiousness *= (float)gen_gamma_mt(1 / (P.InfectiousnessSD * P.InfectiousnessSD), 1 / (P.InfectiousnessSD * P.InfectiousnessSD), tn);
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
			else							a->recovery_or_death_time = a->latent_time + P.infectious_icdf.choose(P.InfectiousPeriod, tn, P.TimeStepsPerDay);
		}
		else
		{
			int CaseTime = a->latent_time + ((int)(P.LatentToSymptDelay / P.TimeStep)); //// base severity times on CaseTime, not latent time. Otherwise there are edge cases where recovery time is zero days after latent_time and therefore before DoCase called in IncubRecoverySweep (i.e. people can recover before they've become a case!).

			//// choose final disease severity (either mild, ILI, SARI, Critical, not asymptomatic as covered above) by age
		// TODO	a->Severity_Final = ChooseFinalDiseaseSeverity(age, tn);

			/// choose outcome recovery or death
			if (((a->Severity_Final == Severity_Critical) && (ranf_mt(tn) < P.CFR_Critical_ByAge[age])) ||
				((a->Severity_Final == Severity_SARI) && (ranf_mt(tn) < P.CFR_SARI_ByAge[age])) ||
				((a->Severity_Final == Severity_ILI) && (ranf_mt(tn) < P.CFR_ILI_ByAge[age])))
				a->to_die = 1;

			//// choose events and event times
			if (a->Severity_Final == Severity_Mild)
				a->recovery_or_death_time = CaseTime + P.MildToRecovery_icdf.choose(P.Mean_MildToRecovery[age], tn, P.TimeStepsPerDay);
			else if (a->Severity_Final == Severity_Critical)
			{
				a->SARI_time = CaseTime + P.ILIToSARI_icdf.choose(P.Mean_ILIToSARI[age], tn, P.TimeStepsPerDay);
				a->Critical_time = a->SARI_time + P.SARIToCritical_icdf.choose(P.Mean_SARIToCritical[age], tn, P.TimeStepsPerDay);
				if (a->to_die)
					a->recovery_or_death_time = a->Critical_time + P.CriticalToDeath_icdf.choose(P.Mean_CriticalToDeath[age], tn, P.TimeStepsPerDay);
				else
				{
					a->RecoveringFromCritical_time = a->Critical_time + P.CriticalToCritRecov_icdf.choose(P.Mean_CriticalToCritRecov[age], tn, P.TimeStepsPerDay);
					a->recovery_or_death_time = a->RecoveringFromCritical_time + P.CritRecovToRecov_icdf.choose(P.Mean_CritRecovToRecov[age], tn, P.TimeStepsPerDay);
				}
			}
			else if (a->Severity_Final == Severity_SARI)
			{
				a->SARI_time = CaseTime + P.ILIToSARI_icdf.choose(P.Mean_ILIToSARI[age], tn, P.TimeStepsPerDay);
				if (a->to_die)
					a->recovery_or_death_time = a->SARI_time + P.SARIToDeath_icdf.choose(P.Mean_SARIToDeath[age], tn, P.TimeStepsPerDay);
				else
					a->recovery_or_death_time = a->SARI_time + P.SARIToRecovery_icdf.choose(P.Mean_SARIToRecovery[age], tn, P.TimeStepsPerDay);
			}
			else /*i.e. if Severity_Final == Severity_ILI*/
			{
				if (a->to_die)
					a->recovery_or_death_time = CaseTime + P.ILIToDeath_icdf.choose(P.Mean_ILIToDeath[age], tn, P.TimeStepsPerDay);
				else
					a->recovery_or_death_time = CaseTime + P.ILIToRecovery_icdf.choose(P.Mean_ILIToRecovery[age], tn, P.TimeStepsPerDay);
			}
		}

		if ((a->inf == InfStat_InfectiousAlmostSymptomatic) && ((P.ControlPropCasesId == 1) || (ranf_mt(tn) < P.ControlPropCasesId)))
		{
			Hosts[ai].detected = 1;
			Hosts[ai].detected_time = ts + (unsigned short int)(P.LatentToSymptDelay * P.TimeStepsPerDay);


			if ((P.DoDigitalContactTracing) && (Hosts[ai].detected_time >= (unsigned short int)(AdUnits[Mcells[Hosts[ai].mcell].adunit].DigitalContactTracingTimeStart * P.TimeStepsPerDay)) && (Hosts[ai].detected_time < (unsigned short int)((AdUnits[Mcells[Hosts[ai].mcell].adunit].DigitalContactTracingTimeStart + P.DigitalContactTracingPolicyDuration) * P.TimeStepsPerDay)) && (Hosts[ai].digitalContactTracingUser))
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

	return a->inf;
}
