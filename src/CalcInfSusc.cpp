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
	return	((HOST_ISOLATED(person) && (Hosts[person].digitalContactTraced != 1)) ? P.Efficacies[CaseIsolation][PlaceType] : 1.0)
		*	((Hosts[person].digitalContactTraced == 1) ? P.Efficacies[DigContactTracing][PlaceType] : 1.0)
		*	((HOST_QUARANTINED(person) && (!Hosts[person].care_home_resident) && (Hosts[person].digitalContactTraced != 1) && (!(HOST_ISOLATED(person)))) ? P.Efficacies[HomeQuarantine][PlaceType] : 1.0)
		*	((Hosts[person].is_case() && (!Hosts[person].care_home_resident)) ? P.SymptPlaceTypeContactRate[PlaceType] : 1.0)
		*	P.PlaceTypeTrans[PlaceType] / P.PlaceTypeGroupSizeParam1[PlaceType] * CalcPersonInf(person, TimeStepNow);
}

double CalcSpatialInf(int person, unsigned short int TimeStepNow)
{
	return	((HOST_ISOLATED(person) && (Hosts[person].digitalContactTraced != 1)) ? P.Efficacies[CaseIsolation][Spatial] : 1.0)
		*	((Hosts[person].digitalContactTraced==1) ? P.Efficacies[DigContactTracing][Spatial] : 1.0)
		*   ((HOST_QUARANTINED(person) && (!Hosts[person].care_home_resident) && (Hosts[person].digitalContactTraced != 1) && (!(HOST_ISOLATED(person)))) ? P.Efficacies[HomeQuarantine][Spatial] : 1.0)
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
double CalcHouseSusc(int person, unsigned short int TimeStepNow, int infector)
{
	return CalcPersonSusc(person, TimeStepNow, infector)
		* ((Mcells[Hosts[person].mcell].socdist == TreatStat::Treated) ? ((Hosts[person].esocdist_comply) ? P.Efficacies[EnhancedSocialDistancing][House] : P.Efficacies[SocialDistancing][House]) : 1.0)
		* ((Hosts[person].digitalContactTraced == 1) ? P.Efficacies[DigContactTracing][House] : 1.0)
		* ((Hosts[person].care_home_resident) ? P.CareHomeResidentHouseholdScaling : 1.0);
}
double CalcPlaceSusc(int person, int PlaceType, unsigned short int TimeStepNow)
{
	return		((HOST_QUARANTINED(person) && (!Hosts[person].care_home_resident) && (Hosts[person].digitalContactTraced != 1)) ? P.Efficacies[HomeQuarantine][PlaceType] : 1.0)
		* ((Mcells[Hosts[person].mcell].socdist == TreatStat::Treated) ? ((Hosts[person].esocdist_comply) ? P.Efficacies[EnhancedSocialDistancing][PlaceType] : P.Efficacies[SocialDistancing][PlaceType]) : 1.0)
		* ((Hosts[person].digitalContactTraced == 1) ? P.Efficacies[DigContactTracing][PlaceType] : 1.0);
}
double CalcSpatialSusc(int person, unsigned short int TimeStepNow)
{
	return	 ((HOST_QUARANTINED(person) && (!Hosts[person].care_home_resident) && (Hosts[person].digitalContactTraced != 1)) ? P.Efficacies[HomeQuarantine][Spatial] : 1.0)
		* ((Mcells[Hosts[person].mcell].socdist == TreatStat::Treated) ? ((Hosts[person].esocdist_comply) ? P.Efficacies[EnhancedSocialDistancing][Spatial] : P.Efficacies[SocialDistancing][Spatial]) : 1.0)
		* ((Hosts[person].digitalContactTraced == 1) ? P.Efficacies[DigContactTracing][Spatial] : 1.0)
		* P.RelativeSpatialContactSusc[HOST_AGE_GROUP(person)];
}
double CalcPersonSusc(int person, unsigned short int TimeStepNow, int infector)
{
	return		P.WAIFW_Matrix[HOST_AGE_GROUP(person)][HOST_AGE_GROUP(infector)]
		* P.AgeSusceptibility[HOST_AGE_GROUP(person)] * Hosts[person].susc
		*	(HOST_TREATED(person) ? P.TreatSuscDrop : 1.0)
		*	(HOST_VACCED(person) ? (HOST_VACCED_SWITCH(person) ? P.VaccSuscDrop2 : P.VaccSuscDrop) : 1.0);
}
