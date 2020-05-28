#ifndef COVIDSIM_DIST_H_INCLUDED_
#define COVIDSIM_DIST_H_INCLUDED_

extern double sinx[361], cosx[361], asin2sqx[1001];
double dist2UTM(double, double, double, double);
double dist2_raw(double, double, double, double);

#endif // COVIDSIM_DIST_H_INCLUDED_
