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

/**
 * An arc minute of latitude along any line of longitude in meters.
 *
 * Also known as the International Nautical Mile.
 *
 * @see https://en.wikipedia.org/wiki/Nautical_mile
 */
constexpr int NMI = 1852;

/**
 * The number of arc minutes in one degree.
 *
 * @see https://en.wikipedia.org/wiki/Minute_and_second_of_arc
 */
constexpr int ARCMINUTES_PER_DEGREE = 60;

/**
 * The number of degrees in a complete rotation.
 *
 * @see https://en.wikipedia.org/wiki/Turn_(angle)
 */
constexpr int DEGREES_PER_TURN = 360;

/**
 * The earth's circumference in meters.
 *
 * The units of cancellation:
 *    meters/minute * minutes/degree * degrees = meters
 */
constexpr int EARTH_CIRCUMFERENCE = NMI * ARCMINUTES_PER_DEGREE * DEGREES_PER_TURN;

/**
 * The earth's diameter in meters.
 */
constexpr double EARTH_DIAMETER = EARTH_CIRCUMFERENCE / PI;

/**
 * The Earth's radius in meters.
 *
 * The previous hardcoded value used 6366707 which was derived from the
 * following formula:
 *
 *     Earth's radius (m) = Earth's circumference / 2 * Pi
 *
 * where Earth's circumference can be derived with the following formula:
 *
 *     Earth's circumference (m) = NMI * ARCMINUTES_PER_DEGREE * DEGREES_PER_TURN
 */
constexpr double EARTHRADIUS = EARTH_DIAMETER / 2;

const int OUTPUT_DIST_SCALE = 1000;
const int MAX_PLACE_SIZE = 20000;
const int MAX_NUM_SEED_LOCATIONS = 10000;

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

// define maximum number of contacts
const int MAX_CONTACTS = 500;

const int MAX_DUR_IMPORT_PROFILE = 10245;

const int MAX_AIRPORTS = 5000;
const int NNA = 10;
const double MAX_DIST_AIRPORT_TO_HOTEL = 200000.0;
const int MIN_HOTELS_PER_AIRPORT = 20;
const int HOTELS_PER_1000PASSENGER = 10;

const int MAX_NUM_INTERVENTION_CHANGE_TIMES = 100; // may want to make this intervention-specifc, but keep simple for now.

#endif // COVIDSIM_CONSTANTS_H_INCLUDED_
