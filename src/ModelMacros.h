#ifndef COVIDSIM_MODELMACROS_H_INCLUDED_
#define COVIDSIM_MODELMACROS_H_INCLUDED_

#include <climits>

#define HOST_TREATED(x)				((Hosts[x].treat_stop_time > TimeStepNow) && (Hosts[x].treat_start_time <= TimeStepNow))
#define HOST_TO_BE_TREATED(x)		(Hosts[x].treat_stop_time > TimeStepNow)
#define PLACE_TREATED(x, y)			(Places[x][y].treat_end_time > TimeStepNow)
#define PLACE_CLOSED(x, y)			((Places[x][y].close_start_time <= TimeStepNow) && (Places[x][y].close_end_time > TimeStepNow))
#define HOST_TO_BE_VACCED(x)		(Hosts[x].vacc_start_time < USHRT_MAX - 1)
#define HOST_VACCED(x)				(Hosts[x].vacc_start_time + P.usVaccTimeToEfficacy <= TimeStepNow)
#define HOST_VACCED_SWITCH(x)		(Hosts[x].vacc_start_time >= P.usVaccTimeEfficacySwitch)
#define HOST_QUARANTINED(x)			((HostsQuarantine[x].comply == 1) && (HostsQuarantine[x].start_time + P.usHQuarantineHouseDuration > TimeStepNow) && (HostsQuarantine[x].start_time <= TimeStepNow))
#define HOST_TO_BE_QUARANTINED(x)	((HostsQuarantine[x].start_time + P.usHQuarantineHouseDuration > TimeStepNow) && (HostsQuarantine[x].comply < 2))
#define HOST_ISOLATED(x)			((Hosts[x].isolation_start_time + P.usCaseIsolationDelay <= TimeStepNow) && (Hosts[x].isolation_start_time + P.usCaseIsolationDelay + P.usCaseIsolationDuration > TimeStepNow))
#define HOST_ABSENT(x)				((Hosts[x].absent_start_time <= TimeStepNow) && (Hosts[x].absent_stop_time > TimeStepNow))

/*
  #define NO_TREAT_PROPH_CASES
*/

#define HOST_AGE_YEAR(x)			(Hosts[x].age)
#define HOST_AGE_GROUP(x)			(Hosts[x].age / AGE_GROUP_WIDTH)

#endif // COVIDSIM_MODELMACROS_H_INCLUDED_
