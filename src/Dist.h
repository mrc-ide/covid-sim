#ifndef COVIDSIM_DIST_H_INCLUDED_
#define COVIDSIM_DIST_H_INCLUDED_

#include "Model.h"

#include "geometry/Geometry.h"

extern double sinx[361], cosx[361], asin2sqx[1001];
double dist2UTM(double, double, double, double);
double dist2(Person*, Person*);
double dist2_cc(Cell*, Cell*);
double dist2_cc_min(Cell*, Cell*);
double dist2_mm(Microcell*, Microcell*);
double dist2_raw(double, double, double, double);

double dist2UTM(const CovidSim::Geometry::Vector2f& p1, const CovidSim::Geometry::Vector2f& p2);
double dist2_raw(const CovidSim::Geometry::Vector2f& p1, const CovidSim::Geometry::Vector2f& p2);

#endif // COVIDSIM_DIST_H_INCLUDED_
