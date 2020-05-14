#pragma once

#include "Model.h"
extern double sinx[361], cosx[361], asin2sqx[1001];
double dist2UTM(double, double, double, double);
double dist2(Person*, Person*);
double dist2_cc(Cell*, Cell*);
double dist2_cc_min(Cell*, Cell*);
double dist2_mm(Microcell*, Microcell*);
double dist2_raw(double, double, double, double);
