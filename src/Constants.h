#ifndef COVIDSIM_CONSTANTS_H_INCLUDED_
#define COVIDSIM_CONSTANTS_H_INCLUDED_

/**
 * Math constant defined as the ratio of a circle's circumference to its diameter.
 *
 * TODO: since all calculations using this constant are being automatically
 * type-casted to double, should the precision be extended for more accuracy in
 * the simulations?
 *
 * Eventually could be replaced with C++20's std::numbers::pi.
 * https://en.cppreference.com/w/cpp/header/numbers
 */
constexpr double PI = 3.1415926535; // full double precision: 3.14159265358979323846

const int OUTPUT_DIST_SCALE = 1000;
const int MAX_PLACE_SIZE = 20000;
const int MAX_NUM_SEED_LOCATIONS = 10000;

const int MAX_CLP_COPIES = 50;

const int CDF_RES = 20;
const int INFPROF_RES = 56;

const int NUM_AGE_GROUPS = 17;
const int AGE_GROUP_WIDTH = 5;
const int DAYS_PER_YEAR = 364;
const int INFECT_TYPE_MASK = 16;
const int MAX_GEN_REC = 20;
const int MAX_SEC_REC = 500;
const int INF_QUEUE_SCALE = 5;
const int MAX_TRAVEL_TIME = 14;

const int MAX_INFECTIOUS_STEPS = 2550;

const int MAX_NUM_THREADS = 96;
const int CACHE_LINE_SIZE = 64;

// define maximum number of contacts
const int MAX_CONTACTS = 500;

const int MAX_DUR_IMPORT_PROFILE = 10245;

const int MAX_AIRPORTS = 5000;
const int NNA = 10;
// Need to use define for MAX_DIST_AIRPORT_TO_HOTEL to avoid differences between GCC and clang in requirements to share const doubles in OpenMP default(none) pragmas
#define MAX_DIST_AIRPORT_TO_HOTEL 200000.0
const int MIN_HOTELS_PER_AIRPORT = 20;
const int HOTELS_PER_1000PASSENGER = 10;

const int MAX_NUM_INTERVENTION_CHANGE_TIMES = 100; // may want to make this intervention-specifc, but keep simple for now.

// MS generates a lot of C26451 overflow warnings. Below is shorthand
// to increase the cast size and clean them up.

#define _I64(x) static_cast<int64_t>(x)

#endif // COVIDSIM_CONSTANTS_H_INCLUDED_
