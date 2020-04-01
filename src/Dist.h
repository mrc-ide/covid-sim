#pragma once

#ifndef SPATIALSIM_DIST_H_INCLUDED_
#define SPATIALSIM_DIST_H_INCLUDED_

#include "Model.h"
extern double sinx[361], cosx[361], asin2sqx[1001];
double dist2UTM(double, double, double, double);
double dist2(person*, person*);
double dist2_cc(cell*, cell*);
double dist2_cc_min(cell*, cell*);
double dist2_mm(microcell*, microcell*);
double dist2_raw(double, double, double, double);

#endif // SPATIALSIM_DIST_H_INCLUDED_
