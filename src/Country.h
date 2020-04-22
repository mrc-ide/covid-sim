#ifndef COVIDSIM_COUNTRY_H_INCLUDED_
#define COVIDSIM_COUNTRY_H_INCLUDED_


#define MAX_HOUSEHOLD_SIZE 10
#define MAX_INTERVENTION_TYPES 1
#define MAX_INTERVENTIONS_PER_ADUNIT 10

const int ADUNIT_LOOKUP_SIZE = 1000000;
const int MAX_COUNTRIES = 100;
const int NUM_PLACE_TYPES = 4;
const int MAX_ADUNITS = 3200;

// Maximal absent time - used for array sizing
const int MAX_ABSENT_TIME = 8;

#ifdef COUNTRY_US
#define ABSENTEEISM_PLACE_CLOSURE
#else
/*
#define ABSENTEEISM_PLACE_CLOSURE
*/
#endif

#ifdef COUNTRY_US
#define MAX_DIST 2000
#define NK_HR 2500
#define NKR 4000000
#else
#define MAX_DIST 2000
#define NK_HR 400
#define NKR 2000000
#endif

#endif // COVIDSIM_COUNTRY_H_INCLUDED_
