#pragma once

#ifndef SPATIALSIM_MODELMACROS_H_INCLUDED_
#define SPATIALSIM_MODELMACROS_H_INCLUDED_

#include <limits.h>

#define HOST_TREATED(x) ((Hosts[x].treat_stop_time > ts) && (Hosts[x].treat_start_time <= ts))
#define HOST_TO_BE_TREATED(x) (Hosts[x].treat_stop_time > ts)
#define PLACE_TREATED(x, y) (Places[x][y].treat_end_time > ts)
#define PLACE_CLOSED(x, y) ((Places[x][y].close_start_time <= ts) && (Places[x][y].close_end_time > ts))
#define HOST_TO_BE_VACCED(x) (Hosts[x].vacc_start_time < USHRT_MAX - 1)
#define HOST_VACCED(x) (Hosts[x].vacc_start_time+P.usVaccTimeToEfficacy<=ts)
#define HOST_VACCED_SWITCH(x) (Hosts[x].vacc_start_time >= P.usVaccTimeEfficacySwitch)
#define HOST_QUARANTINED(x) ((Hosts[x].quar_comply == 1) && (Hosts[x].quar_start_time + P.usHQuarantineHouseDuration > ts) && (Hosts[x].quar_start_time <= ts))
#define HOST_TO_BE_QUARANTINED(x) ((Hosts[x].quar_start_time + P.usHQuarantineHouseDuration > ts) && (Hosts[x].quar_comply < 2))
#define HOST_ISOLATED(x) ((Hosts[x].isolation_start_time + P.usCaseIsolationDelay <= ts) && (Hosts[x].isolation_start_time + P.usCaseIsolationDelay + P.usCaseIsolationDuration > ts))
#define HOST_ABSENT(x) ((Hosts[x].absent_start_time <= ts) && (Hosts[x].absent_stop_time > ts))
//add macro to see if host is isolated due to digital contact tracing
#define HOST_ISOLATED_DCT(x) ((Hosts[x].dct_start_time <= ts) && (Hosts[x].dct_end_time > ts))

/*
  #define NO_TREAT_PROPH_CASES
*/

#define HOST_AGE_YEAR(x) (Hosts[x].age)
#define HOST_AGE_GROUP(x) (Hosts[x].age / AGE_GROUP_WIDTH)

#endif // SPATIALSIM_MODELMACROS_H_INCLUDED_
