#ifndef COVIDSIM_CALCINFSUSC_H_
#define COVIDSIM_CALCINFSUSC_H_

double CalcHouseInf(int, unsigned short int);
double CalcPlaceInf(int, int, unsigned short int);
double CalcSpatialInf(int, unsigned short int);
double CalcPersonInf(int, unsigned short int);
double CalcHouseSusc(int, unsigned short int, int, int);
double CalcPlaceSusc(int, int, unsigned short int, int, int);
double CalcSpatialSusc(int, unsigned short int, int, int);
double CalcPersonSusc(int, unsigned short int, int, int);

#endif // COVIDSIM_CALCINFSUSC_
