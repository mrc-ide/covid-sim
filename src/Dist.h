#ifndef COVIDSIM_DIST_H_INCLUDED_
#define COVIDSIM_DIST_H_INCLUDED_

#include "Models/Person.h"
#include "Models/Cell.h"
#include "Models/Microcell.h"
#include "Constants.h"

double dist2UTM(double, double, double, double);
double dist2(Person*, Person*);
double dist2_cc(Cell*, Cell*);
double dist2_cc_min(Cell*, Cell*);
double dist2_mm(Microcell*, Microcell*);
double dist2_raw(double, double, double, double);
double periodic_xy(double x, double y);

#endif // COVIDSIM_DIST_H_INCLUDED_
