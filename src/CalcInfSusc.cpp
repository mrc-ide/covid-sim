#include <cmath>

#include "CalcInfSusc.h"
#include "Constants.h"
#include "InfStat.h"
#include "Model.h"
#include "ModelMacros.h"
#include "Param.h"

//// Infectiousness functions (House, Place, Spatial, Person). Idea is that in addition to a person's personal infectiousness, they have separate "infectiousnesses" for their house, place and on other cells (spatial).
//// These functions consider one person only. A person has an infectiousness that is independent of other people. Slightly different therefore than susceptibility functions.
double CalcHouseInf(int person, unsigned short int TimeStepNow)
{
	return	((HOST_ISOLATED(person) && (Hosts[person].digitalContactTraced != 1)) ? P.Efficacies[CaseIsolation][House] : 1.0)
		*	((Hosts[person].digitalContactTraced==1) ? P.Efficacies[DigContactTracing][House] : 1.0)
		*	((HOST_QUARANTINED(person) && (Hosts[person].digitalContactTraced != 1) && (!(HOST_ISOLATED(person)))) ? P.Efficacies[HomeQuarantine][House] : 1.0)
		*	P.HouseholdDenomLookup[Households[Hosts[person].hh].nhr - 1]
		*   ((Hosts[person].care_home_resident) ? P.CareHomeResidentHouseholdScaling : 1.0)
		*   (HOST_TREATED(person) ? P.TreatInfDrop : 1.0)
		*   (HOST_VACCED(person) ? P.VaccInfDrop : 1.0)
		*   ((P.NoInfectiousnessSDinHH)? ((Hosts[person].infectiousness < 0) ? P.SymptInfectiousness : P.AsymptInfectiousness):fabs(Hosts[person].infectiousness))  // removed call to CalcPersonInf to allow infectiousness to be const in hh
		*   P.infectiousness[TimeStepNow - Hosts[person].latent_time - 1];
}

double CalcPlaceInf(int person, int PlaceType, unsigned short int TimeStepNow)
{
	return	((HOST_ISOLATED(person) && (Hosts[person].digitalContactTraced != 1)) ? P.CaseIsolationEffectiveness : 1.0)
		*	((Hosts[person].digitalContactTraced==1) ? P.DCTCaseIsolationEffectiveness : 1.0)
		*	((HOST_QUARANTINED(person) && (!Hosts[person].care_home_resident) && (Hosts[person].digitalContactTraced != 1) && (!(HOST_ISOLATED(person)))) ? P.HQuarantinePlaceEffect[PlaceType] : 1.0)
		*	((Hosts[person].is_case() && (!Hosts[person].care_home_resident)) ? P.SymptPlaceTypeContactRate[PlaceType] : 1.0)
		*	P.PlaceTypeTrans[PlaceType] / P.PlaceTypeGroupSizeParam1[PlaceType] * CalcPersonInf(person, TimeStepNow);
}

double CalcSpatialInf(int person, unsigned short int TimeStepNow)
{
	return	((HOST_ISOLATED(person) && (Hosts[person].digitalContactTraced != 1)) ? P.CaseIsolationEffectiveness : 1.0)
		*	((Hosts[person].digitalContactTraced==1) ? P.DCTCaseIsolationEffectiveness : 1.0)
		*   ((HOST_QUARANTINED(person) && (!Hosts[person].care_home_resident) && (Hosts[person].digitalContactTraced != 1) && (!(HOST_ISOLATED(person)))) ? P.HQuarantineSpatialEffect : 1.0)
		*	(Hosts[person].is_case() ? P.SymptSpatialContactRate : 1.0)
		*	P.RelativeSpatialContact[HOST_AGE_GROUP(person)]
		*	CalcPersonInf(person, TimeStepNow); 		/*	*Hosts[person].spatial_norm */
}

double CalcPersonInf(int person, unsigned short int TimeStepNow)
{
	return	(HOST_TREATED(person) ? P.TreatInfDrop : 1.0)
		*	(HOST_VACCED(person) ? P.VaccInfDrop : 1.0)
		*	fabs(Hosts[person].infectiousness)
		*	P.infectiousness[TimeStepNow - Hosts[person].latent_time - 1];
}

//// Susceptibility functions (House, Place, Spatial, Person). Similarly, idea is that in addition to a person's personal susceptibility, they have separate "susceptibilities" for their house, place and on other cells (spatial)
//// These functions consider two people. A person has a susceptibility TO ANOTHER PERSON/infector. Slightly different therefore than infectiousness functions.
double CalcHouseSusc(int ai, unsigned short int TimeStepNow, int infector)
{
	return CalcPersonSusc(ai, TimeStepNow, infector)
		* ((Mcells[Hosts[ai].mcell].socdist == TreatStat::Treated) ? ((Hosts[ai].esocdist_comply) ? P.EnhancedSocDistHouseholdEffectCurrent : P.SocDistHouseholdEffectCurrent) : 1.0)
		* ((Hosts[ai].digitalContactTraced==1) ? P.DCTCaseIsolationHouseEffectiveness : 1.0)
		* ((Hosts[ai].care_home_resident) ? P.CareHomeResidentHouseholdScaling : 1.0);
}
double CalcPlaceSusc(int ai, int PlaceType, unsigned short int TimeStepNow)
{
	return		((HOST_QUARANTINED(ai) && (!Hosts[ai].care_home_resident) && (Hosts[ai].digitalContactTraced != 1)) ? P.HQuarantinePlaceEffect[PlaceType] : 1.0)
		* ((Mcells[Hosts[ai].mcell].socdist == TreatStat::Treated) ? ((Hosts[ai].esocdist_comply) ? P.EnhancedSocDistPlaceEffectCurrent[PlaceType] : P.SocDistPlaceEffectCurrent[PlaceType]) : 1.0)
		* ((Hosts[ai].digitalContactTraced==1) ? P.DCTCaseIsolationEffectiveness : 1.0);
}
double CalcSpatialSusc(int ai, unsigned short int TimeStepNow)
{
	return	 ((HOST_QUARANTINED(ai) && (!Hosts[ai].care_home_resident) && (Hosts[ai].digitalContactTraced != 1)) ? P.HQuarantineSpatialEffect : 1.0)
		* ((Mcells[Hosts[ai].mcell].socdist == TreatStat::Treated) ? ((Hosts[ai].esocdist_comply) ? P.EnhancedSocDistSpatialEffectCurrent : P.SocDistSpatialEffectCurrent) : 1.0)
		* ((Hosts[ai].digitalContactTraced == 1) ? P.DCTCaseIsolationEffectiveness : 1.0)
		* P.RelativeSpatialContactSusc[HOST_AGE_GROUP(ai)];
}
double CalcPersonSusc(int ai, unsigned short int TimeStepNow, int infector)
{
	return		P.WAIFW_Matrix[HOST_AGE_GROUP(ai)][HOST_AGE_GROUP(infector)]
		* P.AgeSusceptibility[HOST_AGE_GROUP(ai)] * Hosts[ai].susc
		*	(HOST_TREATED(ai) ? P.TreatSuscDrop : 1.0)
		*	(HOST_VACCED(ai) ? (HOST_VACCED_SWITCH(ai) ? P.VaccSuscDrop2 : P.VaccSuscDrop) : 1.0);
}
