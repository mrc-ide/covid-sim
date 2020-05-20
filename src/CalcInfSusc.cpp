#include <math.h>

#include "CalcInfSusc.h"
#include "Constants.h"
#include "InfStat.h"
#include "Model.h"
#include "ModelMacros.h"
#include "Param.h"

//// Infectiousness functions (House, Place, Spatial, Person). Idea is that in addition to a person's personal infectiousness, they have separate "infectiousnesses" for their house, place and on other cells (spatial).
//// These functions consider one person only. A person has an infectiousness that is independent of other people. Slightly different therefore than susceptibility functions.
double CalcHouseInf(int j, unsigned short int ts)
{
	return	((HOST_ISOLATED(j) && (Hosts[j].digitalContactTraced != 1)) ? P.CaseIsolationHouseEffectiveness : 1.0)
		*	((Hosts[j].digitalContactTraced==1) ? P.DCTCaseIsolationHouseEffectiveness : 1.0)
		*	((HOST_QUARANTINED(j) && (Hosts[j].digitalContactTraced != 1) && (!(HOST_ISOLATED(j))))? P.HQuarantineHouseEffect : 1.0)
		*	P.HouseholdDenomLookup[Households[Hosts[j].hh].nhr - 1] * CalcPersonInf(j, ts);
}
double CalcPlaceInf(int j, int k, unsigned short int ts)
{
	return	((HOST_ISOLATED(j) && (Hosts[j].digitalContactTraced != 1)) ? P.CaseIsolationEffectiveness : 1.0)
		*	((Hosts[j].digitalContactTraced==1) ? P.DCTCaseIsolationEffectiveness : 1.0)
		*	((HOST_QUARANTINED(j) && (Hosts[j].digitalContactTraced != 1) && (!(HOST_ISOLATED(j)))) ? P.HQuarantinePlaceEffect[k] : 1.0)
		*	((Hosts[j].inf == InfStat_Case) ? P.SymptPlaceTypeContactRate[k] : 1.0)
		*	P.PlaceTypeTrans[k] / P.PlaceTypeGroupSizeParam1[k] * CalcPersonInf(j, ts);
}
double CalcSpatialInf(int j, unsigned short int ts)
{
	return	((HOST_ISOLATED(j) && (Hosts[j].digitalContactTraced != 1)) ? P.CaseIsolationEffectiveness : 1.0)
		*	((Hosts[j].digitalContactTraced==1) ? P.DCTCaseIsolationEffectiveness : 1.0)
		*   ((HOST_QUARANTINED(j) && (Hosts[j].digitalContactTraced != 1) && (!(HOST_ISOLATED(j)))) ? P.HQuarantineSpatialEffect : 1.0)
		*	((Hosts[j].inf == InfStat_Case) ? P.SymptSpatialContactRate : 1.0)
		*	P.RelativeSpatialContact[HOST_AGE_GROUP(j)]
		*	CalcPersonInf(j, ts); 		/*	*Hosts[j].spatial_norm */
}
double CalcPersonInf(int j, unsigned short int ts)
{
	return	(HOST_TREATED(j) ? P.TreatInfDrop : 1.0)
		*	(HOST_VACCED(j) ? P.VaccInfDrop : 1.0)
		*	fabs(Hosts[j].infectiousness)
		*	P.infectiousness[ts - Hosts[j].latent_time - 1];
}

//// Susceptibility functions (House, Place, Spatial, Person). Similarly, idea is that in addition to a person's personal susceptibility, they have separate "susceptibilities" for their house, place and on other cells (spatial)
//// These functions consider two people. A person has a susceptibility TO ANOTHER PERSON/infector. Slightly different therefore than infectiousness functions.
double CalcHouseSusc(int ai, unsigned short int ts, int infector, int tn)
{
	return CalcPersonSusc(ai, ts, infector, tn)
		* ((Mcells[Hosts[ai].mcell].socdist == 2) ? ((Hosts[ai].esocdist_comply) ? P.EnhancedSocDistHouseholdEffectCurrent : P.SocDistHouseholdEffectCurrent) : 1.0)
		* (Hosts[ai].digitalContactTraced==1 ? P.DCTCaseIsolationHouseEffectiveness : 1.0);
}
double CalcPlaceSusc(int ai, int k, unsigned short int ts, int infector, int tn)
{
	return		((HOST_QUARANTINED(ai) && (Hosts[ai].digitalContactTraced != 1)) ? P.HQuarantinePlaceEffect[k] : 1.0)
		* ((Mcells[Hosts[ai].mcell].socdist == 2) ? ((Hosts[ai].esocdist_comply) ? P.EnhancedSocDistPlaceEffectCurrent[k] : P.SocDistPlaceEffectCurrent[k]) : 1.0)
		*	(Hosts[ai].digitalContactTraced==1 ? P.DCTCaseIsolationEffectiveness : 1.0);
}
double CalcSpatialSusc(int ai, unsigned short int ts, int infector, int tn)
{
	return	 ((HOST_QUARANTINED(ai) && (Hosts[ai].digitalContactTraced != 1)) ? P.HQuarantineSpatialEffect : 1.0)
		* ((Mcells[Hosts[ai].mcell].socdist == 2) ? ((Hosts[ai].esocdist_comply) ? P.EnhancedSocDistSpatialEffectCurrent : P.SocDistSpatialEffectCurrent) : 1.0)
		*	(Hosts[ai].digitalContactTraced==1 ? P.DCTCaseIsolationEffectiveness : 1.0);
}
double CalcPersonSusc(int ai, unsigned short int ts, int infector, int tn)
{
	return		P.WAIFW_Matrix[HOST_AGE_GROUP(ai)][HOST_AGE_GROUP(infector)]
		* P.AgeSusceptibility[HOST_AGE_GROUP(ai)] * Hosts[ai].susc
		*	(HOST_TREATED(ai) ? P.TreatSuscDrop : 1.0)
		*	(HOST_VACCED(ai) ? (HOST_VACCED_SWITCH(ai) ? P.VaccSuscDrop2 : P.VaccSuscDrop) : 1.0);
}
